import tensorflow as tf
#import tensorflow_probability as tfp

tfkl = tf.keras.layers
#tfpl = tfp.layers
#tfd = tfp.distributions

class KLWeightCallback(tf.keras.callbacks.Callback):
    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        self._model = model

    def __init__(self, model, increase_per_epoch):
        super(KLWeightCallback, self).__init__()
        self.model = model
        self.increase_per_epoch = increase_per_epoch

    def on_epoch_end(self, epoch, logs=None):
        self.model.kl_weight.assign_add(self.increase_per_epoch)
        self.model.kl_weight.assign(min(self.model.kl_weight, 1.0))
        


class VAE(tf.keras.Model):
    def __init__(self, latent_dim=50, nconditions=4, grid_size=10, weight=0.0, increase_per_epoch=0.1 ):
        super(VAE, self).__init__()
        self.latent_dim = latent_dim
        self.nconditions = nconditions
        self.flat_shape = grid_size*grid_size*2
        self.kl_weight = tf.Variable(weight, trainable=False)

        self.conditions_encoder = tf.keras.Sequential([
            tfkl.InputLayer(shape=(self.nconditions,)),
            tfkl.Dense(64, activation='tanh'),
            tfkl.Dense(64, activation='tanh'),
            tfkl.Dense(64, activation='tanh'),
            tfkl.Dense(self.latent_dim * 2)
        ])

        self.images_encoder = tf.keras.Sequential([
            tfkl.InputLayer(shape=(self.flat_shape,)),
            tfkl.Dense(self.flat_shape, activation='tanh'),
            tfkl.Dense(64, activation='tanh'),
            tfkl.Dense(64, activation='tanh'),
            tfkl.Dense(64, activation='tanh'),
            tfkl.Dense(self.latent_dim)
        ])        
        
        self.decoder = tf.keras.Sequential([            
            tfkl.InputLayer(shape=(latent_dim,),name='input_layer'),
            tfkl.Dense(64, activation='relu'),
            tfkl.Dense(64, activation='relu'),
            tfkl.Dense(self.flat_shape, activation='relu'),
            tfkl.Dense(self.flat_shape, name='output_layer')
        ])

    def on_epoch_end(self, epoch, logs=None):
        self.weight += self.increase_per_epoch
        self.weight = min(self.weight, 1.0)

    
    def encode_conditions(self, conditions):
        mean, logvar = tf.split(self.conditions_encoder(conditions), num_or_size_splits=2, axis=1)
        return mean, logvar

    def encode_images(self, images):
        return self.images_encoder(images)
    
    def reparameterize(self, mean, logvar):
        eps = tf.random.normal(shape=tf.shape(mean))
        return eps * tf.exp(logvar * .5) + mean
    
    def decode(self, z):
        return self.decoder(z)

    def call(self, inputs, kl_weight=0.0):   
        conditions = inputs[:,0:self.nconditions]
        x = inputs[:,self.nconditions:]

        encoded_images = self.encode_images(x)
        mean, logvar = self.encode_conditions(conditions)
        z = self.reparameterize(mean, logvar)

        reconstruction = self.decode(encoded_images)
        reconstruction_z = self.decode(z)

        # # Reconstruction loss
        # reconstruction_loss = tf.reduce_mean(tf.reduce_sum(tf.math.squared_difference(x, reconstruction)))
        # # KL divergence
        # kl_loss = -0.5 * tf.reduce_sum(1 + logvar - tf.square(mean) - tf.exp(logvar), axis=1)*kl_weight
        # # Latent loss between predictions of condition and images
        # latent_loss = tf.reduce_mean(tf.reduce_sum(tf.math.squared_difference(mean, encoded_images)))
        
        # Square diff
        recon_diff = tf.math.squared_difference(x, reconstruction)
        recon_diff_z = tf.math.squared_difference(x, reconstruction_z)

        # Multiply diffs together
        recon_loss = tf.reduce_mean(recon_diff*recon_diff_z, axis=1)

        # Reconstruction loss
        # reconstruction_loss = tf.reduce_mean(tf.math.squared_difference(x, reconstruction), axis=1)
        # reconstruction_loss_z = tf.reduce_mean(tf.math.squared_difference(x, reconstruction_z), axis=1)*0.01
        # KL divergence
        kl_loss = -0.5 * tf.reduce_sum(1 + logvar - tf.square(mean) - tf.exp(logvar), axis=1)*kl_weight
        # Latent loss between predictions of condition and images
        #latent_loss = tf.reduce_mean(tf.math.squared_difference(z, encoded_images), axis=1) #Might still be good
        #tf.print('reconstruction_loss',tf.reduce_mean(reconstruction_loss))
        #tf.print('kl_loss',tf.reduce_mean(kl_loss))
        #tf.print('latent_loss',tf.reduce_mean(latent_loss))

        #total_loss = tf.reduce_mean(reconstruction_loss + kl_loss + reconstruction_loss_z)
        total_loss = tf.reduce_mean(recon_loss + kl_loss)
        # add loss
        self.add_loss(total_loss)
        return reconstruction
    

class Generator(tf.keras.Model):
    def __init__(self, original_model):
        super(Generator, self).__init__()
        self.conditions_encoder = original_model.conditions_encoder
        self.decoder    = original_model.decoder
        self.latent_dim = original_model.latent_dim
        self.nconditions = original_model.nconditions
        self.output_names = ['output']

        self.input_layer = tfkl.InputLayer(shape=(self.latent_dim,),name='input_layer')
        
    def reparameterize(self, mean, logvar):
        eps = tf.random.normal(shape=tf.shape(mean))
        return eps * tf.exp(logvar * .5) + mean
    

    def reconstruct(self, conditions):      
        mean, logvar = tf.split(self.conditions_encoder(conditions), num_or_size_splits=2, axis=1)
        z = self.reparameterize(mean, logvar)
        decoded = self.decoder(z)
        return decoded

    def call(self, conditions):
        reconstructed = self.reconstruct(conditions)
        return reconstructed