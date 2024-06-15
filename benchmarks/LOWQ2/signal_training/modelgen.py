import tensorflow as tf
import tensorflow_probability as tfp

tfkl = tf.keras.layers
tfpl = tfp.layers
tfd = tfp.distributions

class KLWeightCallback(tf.keras.callbacks.Callback):
    def __init__(self, model, increase_per_epoch):
        super().__init__()
        self.model = model
        self.increase_per_epoch = increase_per_epoch

    def on_epoch_end(self, epoch, logs=None):
        self.model.kl_weight.assign_add(self.increase_per_epoch)
        self.model.kl_weight.assign(min(self.model.kl_weight, 1.0))
        


class VAE(tf.keras.Model):
    def __init__(self, latent_dim=400, input_shape=(4,), output_shape=(10, 10, 2), weight=0.0, increase_per_epoch=0.1 ):
        super(VAE, self).__init__()
        self.latent_dim = latent_dim        
        self.decoder = tf.keras.Sequential([
            tfkl.InputLayer(input_shape=(latent_dim,)),
            tfkl.Dense(64, activation='relu'),
            tfkl.Dense(256, activation='relu'),
            tfkl.Dense(tf.reduce_prod(output_shape), activation='relu'),
            tfkl.Reshape(output_shape),
        ])

    def on_epoch_end(self, epoch, logs=None):
        self.weight += self.increase_per_epoch
        self.weight = min(self.weight, 1.0)

    def encode(self, x):
        mean_logvar = self.encoder(x)
        N = mean_logvar.shape[-1] // 2
        mean, logvar = mean_logvar[:, :N], mean_logvar[:, N:]
        return mean, logvar
    
    def reparameterize(self, mean, logvar):
        eps_shape = tf.shape(mean)  # Get the dynamic shape at runtime
        eps = tf.random.normal(shape=eps_shape)
        return eps * tf.exp(logvar * 0.5) + mean
        # # Convert logvar to rate parameter of Poisson distribution
        # rate = tf.exp(-logvar)
        # # Sample from Uniform distribution
        # u = tf.random.uniform(shape=tf.shape(mean), minval=0, maxval=1)
        # # Inverse transform sampling to generate Exponential random variable
        # return -tf.math.log(u) / rate + mean
       
    def decode(self, z):
        return self.decoder(z)


    def call(self, inputs, kl_weight=0.0):   
        
        mean   = tf.Variable(tf.random.normal(shape=(self.latent_dim,)))
        logvar = tf.Variable(tf.random.normal(shape=(self.latent_dim,)))

        z = tf.random.normal(shape=(self.latent_dim,))
        
        reconstruction = self.decode(z)
        # Reconstruction loss
        #reconstruction_loss = tf.reduce_mean(tf.reduce_sum(tf.math.squared_difference(y, reconstruction)))
        # KL divergence
        kl_loss = -0.5 * tf.reduce_sum(1 + logvar - tf.square(mean) - tf.exp(logvar), axis=1)
        #kl_loss = self.kl_divergence_exponential(mean, logvar)
        #kl_loss = lambda: -0.5 * tf.reduce_mean(1 + logvar - tf.square(mean) - tf.exp(logvar), axis=1)
        #total_loss = tf.reduce_mean(reconstruction_loss + kl_loss)
        self.add_loss(kl_weight*tf.reduce_mean(kl_loss))
        return reconstruction