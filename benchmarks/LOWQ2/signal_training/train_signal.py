import numpy as np
import uproot
import pandas as pd
import tensorflow as tf
from tensorflow.keras.layers import Input, LSTM, Dense, Lambda, RepeatVector, TimeDistributed
from tensorflow.keras.models import Model
from tensorflow.keras import backend as K
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

# Load data from the ROOT file
file_path = 'output/Out.realistic.root'

# Assuming the ROOT file structure: MCParticles and PixelHits trees
with uproot.open(file_path) as file:
    mc_particles_tree = file['MCParticle/mydetector']
    pixel_hits_tree = file['PixelHit/mydetector']

    print(f'MCParticles columns: {mc_particles_tree.keys()}')
    print(pixel_hits_tree)
    print(f'PixelHits columns: {pixel_hits_tree.keys()}')

    # Load data into DataFrames
    mc_particles_df = mc_particles_tree.arrays(library='np')
    print(mc_particles_df)
    print(mc_particles_df["mydetector"])
    print(dir(mc_particles_df["mydetector"]))
    
    print(mc_particles_df["mydetector"].flatten())
    for value in mc_particles_df["mydetector"]:
        print(value)
        print(dir(value))
        print("BLABABASADSDAS")
        for vec in value:
            print(vec)
            print(dir(vec))
            print(vec.member_names)
            print(vec.getLocalStartPoint())

    mc_particles_df = mc_particles_tree.arrays(library='pd')
    pixel_hits_df = pixel_hits_tree.arrays(library='pd')

# Select relevant columns from MCParticles and PixelHits
mc_particles_data = mc_particles_df[['momentum_x', 'momentum_y', 'momentum_z', 'position_x', 'position_y', 'position_z']].values

# Extract PixelHits and pad sequences to the same length
pixel_hits_data = pixel_hits_df[['charge', 'time', 'pixel_x', 'pixel_y']].groupby(level=0).apply(lambda x: x.values.tolist())
pixel_hits_data = tf.keras.preprocessing.sequence.pad_sequences(pixel_hits_data, dtype='float32', padding='post')

# Ensure consistent length
max_sequence_length = pixel_hits_data.shape[1]

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(mc_particles_data, pixel_hits_data, test_size=0.2, random_state=42)

# Define the VAE model
input_dim = X_train.shape[1]
timesteps = max_sequence_length
output_dim = pixel_hits_data.shape[2]
latent_dim = 16

# Encoder
inputs = Input(shape=(input_dim,))
h = Dense(128, activation='relu')(inputs)
h = Dense(64, activation='relu')(h)
h = Dense(32, activation='relu')(h)

z_mean = Dense(latent_dim)(h)
z_log_var = Dense(latent_dim)(h)

# Sampling function
def sampling(args):
    z_mean, z_log_var = args
    epsilon = K.random_normal(shape=(K.shape(z_mean)[0], latent_dim), mean=0., stddev=1.0)
    return z_mean + K.exp(z_log_var / 2) * epsilon

z = Lambda(sampling, output_shape=(latent_dim,))([z_mean, z_log_var])

# Decoder
decoder_h = Dense(32, activation='relu')
decoder_h2 = Dense(64, activation='relu')
decoder_h3 = Dense(128, activation='relu')
repeat_latent = RepeatVector(timesteps)
decoder_lstm = LSTM(128, return_sequences=True, activation='relu')
decoder_mean = TimeDistributed(Dense(output_dim, activation='linear'))

h_decoded = decoder_h(z)
h_decoded = decoder_h2(h_decoded)
h_decoded = decoder_h3(h_decoded)
h_decoded = repeat_latent(h_decoded)
h_decoded = decoder_lstm(h_decoded)
outputs = decoder_mean(h_decoded)

# Define the VAE model
vae = Model(inputs, outputs)

# Define VAE loss
def vae_loss(inputs, outputs):
    reconstruction_loss = tf.reduce_mean(tf.square(inputs - outputs))
    kl_loss = -0.5 * tf.reduce_mean(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var))
    return reconstruction_loss + kl_loss

vae.compile(optimizer='adam', loss=vae_loss)

# Train the VAE
history = vae.fit(X_train, y_train, epochs=50, batch_size=32, validation_split=0.2, verbose=1)

# Evaluate the VAE on the test set
loss = vae.evaluate(X_test, y_test, verbose=1)
print(f'Test loss: {loss}')

# Save the trained model
vae.save('vae_allpix2_model.h5')

# Predict on test data
predictions = vae.predict(X_test)

# Plot training history
plt.plot(history.history['loss'], label='Training loss')
plt.plot(history.history['val_loss'], label='Validation loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.show()
