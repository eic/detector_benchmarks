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
from model_ad import VAE
from model_ad import Generator
from model_ad import LatentSpace

epochs = 200
batch_size = 5000
model_path = 'model_digitization_genprop'
data_grid_size = 6

condition_columns = ['x', 'y', 'px', 'py']
#condition_columns = ['x', 'y']
nconditions = len(condition_columns)

nInput = nconditions + data_grid_size*data_grid_size*2

# Load data from the ROOT file
file_path = 'output/Out_Convert_genprop.root'

#vae = create_model()
vae = VAE(latent_dim=10,nconditions=nconditions,grid_size=data_grid_size)

#vae.compile(optimizer=Adam())
vae.compile(r_optimizer=Adam(),a_optimizer=Adam())

# Assuming the ROOT file structure: MCParticles and PixelHits trees
with uproot.open(file_path) as file:
    tree = file['events']

    # Extracting data from the ROOT file
    #df = tree.arrays(['x', 'y', 'px', 'py', 'pixel_x', 'pixel_y', 'charge', 'time'], library='pd')
    df = tree.arrays(['x', 'y', 'px', 'py', 'charge', 'time'], entry_stop=1000000)
    
    target_tensors = np.concatenate([df['charge'].to_numpy(), df['time'].to_numpy()], axis=1)
  
    # Create input tensors
    conditions_tensors = df[condition_columns].to_numpy()
    conditions_tensors = np.array([list(t) for t in conditions_tensors])
  
    # Concatenate the conditions_tensors and target_tensors along final axis
    input_tensors = np.concatenate([conditions_tensors, target_tensors], axis=1)

    # Split the input and target tensors into training and validation sets
    input_train, input_val, target_train, target_val = train_test_split(input_tensors, target_tensors, test_size=0.25)
    
    vae.adapt(input_train)

    # Now you can use these variables in the fit function
    #vae.fit(input_train,target_train, validation_data=[input_val, target_val], epochs=epochs, batch_size=batch_size, callbacks=[callback])
    vae.fit(input_train,target_train, validation_data=[input_val, target_val], epochs=epochs, batch_size=batch_size)
    
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

    latent_encoder  = LatentSpace(vae)
    outTest_latent = latent_encoder(input_tensors[:1])
    print(outTest_latent)


    input_signature_latent = [tf.TensorSpec([None,nInput], tf.float32, name='x')]
    
    onnx_model_latent, _ = tf2onnx.convert.from_keras(latent_encoder,input_signature_latent, opset=13)
    onnx.save(onnx_model_latent, model_path+"_latent.onnx")