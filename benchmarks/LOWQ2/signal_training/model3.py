import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions

import collections

collections.Mapping = collections.abc.Mapping

collections.Sequence = collections.abc.Sequence


def create_model(input_shape=(2,), encoded_size=20, output_shape=(10, 10, 2)):
    model = tf.keras.models.Sequential([
        tf.keras.layers.InputLayer(input_shape=input_shape),
        tf.keras.layers.Dense(64, activation='tanh'),
        tf.keras.layers.BatchNormalization(),
        tf.keras.layers.Dense(encoded_size),
        tf.keras.layers.Dense(encoded_size),
        tf.keras.layers.Lambda(lambda x: tfd.Independent(tfd.Normal(loc=x, scale=1))),
        tf.keras.layers.Lambda(lambda x: x.sample()),
        tf.keras.layers.Dense(64, activation='tanh'),
        tf.keras.layers.BatchNormalization(),
        tf.keras.layers.Dense(tf.reduce_prod(output_shape), activation='relu'),
        tf.keras.layers.Reshape(output_shape),
    ])
    return model
