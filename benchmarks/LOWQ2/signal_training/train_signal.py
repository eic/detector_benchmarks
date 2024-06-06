import numpy as np
import uproot
import pandas as pd
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import train_test_split
from tensorflow.keras.preprocessing.sequence import pad_sequences
import tensorflow as tf

# train_signal.py
from model import create_model

# Load data from the ROOT file
file_path = 'output/Out_Convert.root'

vae = create_model()

def mse_round_error(y_true, y_pred):
    mse = tf.reduce_mean(tf.square(y_true - y_pred))
    round_error = tf.reduce_mean(tf.square(y_pred - tf.round(y_pred)))
    return mse + round_error

# Compile the model
vae.compile(optimizer=Adam(), loss="mse")
#vae.compile(optimizer=Adam(), loss=mse_round_error)


# Assuming the ROOT file structure: MCParticles and PixelHits trees
with uproot.open(file_path) as file:
    tree = file['events']

    # Extracting data from the ROOT file
    df = tree.arrays(['x', 'y', 'pixel_x', 'pixel_y', 'charge', 'time'], library='pd')

    #limit the number of rows
    #df = df.head(10000)
    

    # Normalize the 'x' and 'y' columns
    df['x'] = (df['x'] - df['x'].min()) / (df['x'].max() - df['x'].min())
    df['y'] = (df['y'] - df['y'].min()) / (df['y'].max() - df['y'].min())

    # Define a function to create a sparse tensor from a row
    def row_to_sparse_tensor(row):
        charge_indices = np.column_stack([row['pixel_x'], row['pixel_y'], np.zeros(len(row['pixel_x']))])
        time_indices = np.column_stack([row['pixel_x'], row['pixel_y'], np.ones(len(row['pixel_x']))])
        indices = np.concatenate([charge_indices, time_indices])
        values = np.concatenate([row['charge'], row['time']])
        dense_shape = [10, 10, 2]
        sparse_tensor = tf.sparse.reorder(tf.sparse.SparseTensor(indices, values, dense_shape))
        return tf.sparse.to_dense(sparse_tensor)

    # Apply the function to each row of the DataFrame
    #target_tensors = df.apply(row_to_sparse_tensor, axis=1)
    target_tensors = tf.stack(df.apply(row_to_sparse_tensor, axis=1).to_list())

    # Create input tensors
    input_tensors = df[['x', 'y']].to_numpy()


    target_tensors_np = target_tensors.numpy()

    # Split the input and target tensors into training and validation sets
    input_train, input_val, target_train, target_val = train_test_split(input_tensors, target_tensors_np, test_size=0.2)
    
    # Now you can use these variables in the fit function
    vae.fit(input_train, target_train, validation_data=(input_val, target_val), epochs=100, batch_size=64)

    model_path = 'model.h5'
    vae.save(model_path)