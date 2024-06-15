import numpy as np
import uproot
import pandas as pd
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import train_test_split
from tensorflow.keras.preprocessing.sequence import pad_sequences
import tensorflow as tf
import tf2onnx
import onnx

# train_signal.py
from model2 import VAE
from model2 import KLWeightCallback
from model2 import Generator

epochs = 100
batch_size = 500
model_path = 'model_tpx4_new2'
data_grid_size = 6

condition_columns = ['x', 'y', 'px', 'py']
nconditions = len(condition_columns)

# Load data from the ROOT file
file_path = 'output/Out_Convert_tpx4-6.root'

#vae = create_model()
vae = VAE(latent_dim=2,nconditions=nconditions,weight=1.0,grid_size=data_grid_size)

vae.compile(optimizer=Adam())

# Assuming the ROOT file structure: MCParticles and PixelHits trees
with uproot.open(file_path) as file:
    tree = file['events']

    # Extracting data from the ROOT file
    df = tree.arrays(['x', 'y', 'px', 'py', 'pixel_x', 'pixel_y', 'charge', 'time'], library='pd')

    #limit the number of rows
    #df = df.head(10000)
    

    # Normalize the 'x', 'y', 'px', and 'py' columns
    #df['x']  = (df['x'] - df['x'].min()) / (df['x'].max() - df['x'].min())
    #df['y']  = (df['y'] - df['y'].min()) / (df['y'].max() - df['y'].min())
    #df['px'] = (df['px'] - df['px'].min()) / (df['px'].max() - df['px'].min())
    #df['py'] = (df['py'] - df['py'].min()) / (df['py'].max() - df['py'].min())

    # Define a function to create a sparse tensor from a row
    def row_to_sparse_tensor(row):
        charge_indices = np.column_stack([row['pixel_x'], row['pixel_y'], np.zeros(len(row['pixel_x']))])
        time_indices = np.column_stack([row['pixel_x'], row['pixel_y'], np.ones(len(row['pixel_x']))])
        indices = np.concatenate([charge_indices, time_indices])
        values = np.concatenate([row['charge'], row['time']])
        dense_shape = [data_grid_size, data_grid_size, 2]
        sparse_tensor = tf.sparse.reorder(tf.sparse.SparseTensor(indices, values, dense_shape))
        return tf.sparse.to_dense(sparse_tensor)

    # Apply the function to each row of the DataFrame
    #target_tensors = df.apply(row_to_sparse_tensor, axis=1)
    target_tensors = tf.stack(df.apply(row_to_sparse_tensor, axis=1).to_list())

    # Reshape the target_tensors so that other than the event dimension, the other dimensions are the flat
    target_tensors = tf.reshape(target_tensors, (target_tensors.shape[0], -1)).numpy()

    # Create input tensors
    conditions_tensors = df[condition_columns].to_numpy()
    #conditions_tensors = df[['x', 'y']].to_numpy()
    #conditions_tensors = df[[]].to_numpy()

    # Concatenate the conditions_tensors and target_tensors along final axis
    input_tensors = np.concatenate([conditions_tensors, target_tensors], axis=1)

    # Split the input and target tensors into training and validation sets
    input_train, input_val, target_train, target_val = train_test_split(input_tensors, target_tensors, test_size=0.2)
    
    callback = KLWeightCallback(vae, 0.02)

    # Now you can use these variables in the fit function
    vae.fit(input_train,target_train, validation_data=[input_val, target_val], epochs=epochs, batch_size=batch_size, callbacks=[callback])
    #vae.fit(input_train,target_train, validation_data=[input_val, target_val], epochs=epochs, batch_size=batch_size)
    
    model_name = model_path+'.keras'
    vae.save(model_name)
    
    decoder = Generator(vae)

    outTest = decoder(conditions_tensors[:1])
    #outTest = decoder.predict(conditions_tensors[:1])
    print(outTest)

    input_signature = [tf.TensorSpec([None,nconditions], tf.float32, name='x')]
    
    # Convert the model
    onnx_model, _ = tf2onnx.convert.from_keras(decoder,input_signature, opset=13)
    onnx.save(onnx_model, model_path+".onnx")
