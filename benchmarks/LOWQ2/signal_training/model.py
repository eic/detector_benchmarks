import tensorflow as tf
import tensorflow.keras as tfk
import tensorflow_probability as tfp
from tensorflow.keras.layers import Input, Dense, Concatenate, Lambda, Reshape, BatchNormalization
from tensorflow.keras.models import Model

tfkl = tf.keras.layers
tfpl = tfp.layers
tfd = tfp.distributions

import collections

collections.Mapping = collections.abc.Mapping

collections.Sequence = collections.abc.Sequence

# def create_model():
#     # Define the input
#     inputs = Input(shape=(2,))  # Static size input

#     # Process the input
#     x = Dense(64, activation='tanh')(inputs)
#     x = BatchNormalization()(x)  # Batch normalization to introduce stability


#     # Generate a latent variable
#     latent_dim = 20
#     # Encoder to map inputs to latent distribution parameters
#     z_mean = Dense(latent_dim)(x)
#     z_log_var = Dense(latent_dim)(x)
    
#     # Sample from the latent distribution
#     def sampling(args):
#         z_mean, z_log_var = args
#         batch = tf.shape(z_mean)[0]
#         dim = tf.shape(z_mean)[1]
#         epsilon = tf.random.normal(shape=(batch, dim))
#         return z_mean + tf.exp(0.5 * z_log_var) * epsilon

#     z = Lambda(sampling)([z_mean, z_log_var])

#     # Concatenate the processed input and the latent variable
#     x = Concatenate()([x, z])

#     x = Dense(64, activation='tanh')(x)
#     x = BatchNormalization()(x)  # Batch normalization to introduce stability

#     # Output a 10x10 grid of times and charges
#     outputs = Dense(10 * 10 * 2, activation='relu')(x)  # time and charge for each pixel
#     outputs = tf.keras.layers.Reshape((10, 10, 2))(outputs)

#     # Create the model
#     model = Model(inputs, outputs)

#     return model, z_mean, z_log_var

def create_encoder(input_shape=(2,),encoded_size=20):
    prior = tfd.Independent(tfd.Normal(loc=tf.zeros(encoded_size), scale=1),
                        reinterpreted_batch_ndims=1)
    encoder = tfk.Sequential([
        tfkl.InputLayer(input_shape=input_shape),
        tfkl.Dense(64, activation='tanh'),
        tfkl.BatchNormalization(),
        tfkl.Dense(encoded_size),
        tfkl.Dense(tfpl.MultivariateNormalTriL.params_size(encoded_size),
               activation=None),        
        tfpl.MultivariateNormalTriL(
            encoded_size,
            activity_regularizer=tfpl.KLDivergenceRegularizer(prior)),
    ])

    return encoder

def create_decoder(output_shape=(10,10,2),encoded_size=20):
    decoder = tfk.Sequential([
        tfkl.InputLayer(input_shape=(encoded_size,)),
        tfkl.Dense(64, activation='tanh'),
        tfkl.BatchNormalization(),
        tfkl.Dense(tf.reduce_prod(output_shape), activation='tanh'),
        tfkl.Reshape(output_shape),
    ])

    return decoder


def create_model():
    
    encoder = create_encoder()
    decoder = create_decoder()
    vae = tfk.Model(inputs=encoder.inputs,
                outputs=decoder(encoder.outputs[0]))
    
    return vae