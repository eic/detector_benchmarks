import tensorflow as tf
from tensorflow.keras.layers import Input, Dense, LSTM, Concatenate, Lambda
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Conv2D, Flatten, Reshape, BatchNormalization

'''
def sample_gumbel(shape, eps=1e-20): 
    """Sample from Gumbel(0, 1)"""
    U = tf.random.uniform(shape,minval=0,maxval=1)
    return -tf.math.log(-tf.math.log(U + eps) + eps)

def gumbel_softmax_sample(logits, temperature): 
    """ Draw a sample from the Gumbel-Softmax distribution"""
    y = logits + sample_gumbel(tf.shape(logits))
    return tf.nn.softmax( y / temperature)

def gumbel_softmax(logits, temperature, hard=False):
    """Sample from the Gumbel-Softmax distribution and optionally discretize."""
    y = gumbel_softmax_sample(logits, temperature)
    if hard:
        k = tf.shape(logits)[-1]
        y_hard = tf.cast(tf.one_hot(tf.argmax(y, -1), k), y.dtype)
        y = tf.stop_gradient(y_hard - y) + y
    return y
'''

def create_model():
    # Define the input
    inputs = Input(shape=(2,))  # Static size input

    # Process the input
    x = Dense(64, activation='tanh')(inputs)
    x = BatchNormalization()(x)  # Batch normalization to introduce stability


    # Generate a latent variable
    latent_dim = 20
    # Encoder to map inputs to latent distribution parameters
    z_mean = Dense(latent_dim)(x)
    z_log_var = Dense(latent_dim)(x)
    
    # Sample from the latent distribution
    def sampling(args):
        z_mean, z_log_var = args
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.random.normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon

    z = Lambda(sampling)([z_mean, z_log_var])

    # Concatenate the processed input and the latent variable
    x = Concatenate()([x, z])

    x = Dense(64, activation='tanh')(x)
    x = BatchNormalization()(x)  # Batch normalization to introduce stability

    # Output a 10x10 grid of times and charges
    outputs = Dense(10 * 10 * 2, activation='relu')(x)  # time and charge for each pixel
    outputs = tf.keras.layers.Reshape((10, 10, 2))(outputs)

    # Apply Gumbel-Softmax
    #outputs = Lambda(lambda x: gumbel_softmax(x, temperature=0.5, hard=True))(outputs)
        
    # Create the model
    model = Model(inputs, outputs)

    return model