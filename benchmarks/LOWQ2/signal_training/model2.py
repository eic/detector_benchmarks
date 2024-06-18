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
        self.model.kl_weight.assign(min(self.model.kl_weight, 0.01))
        
class Encoder(tf.keras.Model):
    def __init__(self, nconditions, latent_dim):
        super(Encoder, self).__init__()

        self.image_encoder = tf.keras.Sequential([
            tfkl.InputLayer(shape=(self.flat_shape+self.nconditions,)),
            tfkl.Dense(self.flat_shape, activation='relu'),
            tfkl.Dense(64, activation='relu'),
            tfkl.Dense(2*latent_dim),
        ])

    def call(self, inputs):
        conditions = inputs[:,0:self.nconditions]
        x = inputs[:,self.nconditions:]
        mean, logvar = self.encode(inputs)
        return mean, logvar

class VAE(tf.keras.Model):
    def __init__(self, latent_dim=50, nconditions=4, grid_size=10, weight=0.0, increase_per_epoch=0.1 ):
        super(VAE, self).__init__()
        self.latent_dim = latent_dim
        self.nconditions = nconditions
        self.flat_shape = grid_size*grid_size*2
        self.kl_weight = tf.Variable(weight, trainable=False)
        
        self.image_encoder = tf.keras.Sequential([
            tfkl.InputLayer(shape=(self.flat_shape+self.nconditions,)),
            tfkl.Dense(self.flat_shape, activation='relu'),
            tfkl.Dense(64, activation='relu'),
            tfkl.Dense(2*latent_dim),
        ])
          
        
        self.decoder = tf.keras.Sequential([            
            tfkl.InputLayer(shape=(latent_dim+nconditions,),name='input_layer'),
            tfkl.Dense(64, activation='relu'),
            tfkl.Dense(self.flat_shape, activation='relu'),
            tfkl.Dense(self.flat_shape, name='output_layer')
        ])

    def on_epoch_end(self, epoch, logs=None):
        self.weight += self.increase_per_epoch
        self.weight = min(self.weight, 1.0)

    def encode(self, x):
        mean_logvar = self.image_encoder(x)
        N = mean_logvar.shape[-1] // 2
        mean, logvar = mean_logvar[:, :N], mean_logvar[:, N:]
        return mean, logvar
        
    def reparameterize(self, mean, logvar):
        eps = tf.random.normal(shape=tf.shape(mean))
        return eps * tf.exp(logvar * .5) + mean
    
    def decode(self, conditions, z):
        z = tf.concat([conditions,z], axis=-1)
        return self.decoder(z)
    

    def rbf_kernel(self, x, y, gamma=1.0):
        dist = tf.norm(tf.expand_dims(x, -2) - tf.expand_dims(y, -1), axis=2)
        kxy = tf.exp(-gamma * dist)
        return kxy

    def compute_mmd(self, z, condition, gamma=1.0): 
        z_expanded = tf.expand_dims(z, 1)
        condition_expanded = tf.expand_dims(condition, 0)
        z_tile = tf.tile(z_expanded, [1, tf.shape(condition)[0], 1])
        condition_tile = tf.tile(condition_expanded, [tf.shape(z)[0], 1, 1])

        print(z_tile.shape, condition_tile.shape)

        kzz = self.rbf_kernel(z_tile, z_tile, gamma)
        kcc = self.rbf_kernel(condition_tile, condition_tile, gamma)
        kzc = self.rbf_kernel(z_tile, condition_tile, gamma)
        mmd = tf.reduce_mean(kzz) + tf.reduce_mean(kcc) - 2 * tf.reduce_mean(kzc)
        return mmd

    def call(self, inputs, kl_weight=0.0):   
        conditions = inputs[:,0:self.nconditions]
        x = inputs[:,self.nconditions:]
        
        mean, logvar = self.encode(inputs)
        z = self.reparameterize(mean, logvar)
        reconstruction = self.decode(conditions,mean)

        # # Reconstruction loss
        # reconstruction_loss = tf.reduce_mean(tf.reduce_sum(tf.math.squared_difference(x, reconstruction)))
        # # KL divergence
        # kl_loss = -0.5 * tf.reduce_sum(1 + logvar - tf.square(mean) - tf.exp(logvar), axis=1)*kl_weight
        # # Latent loss between predictions of condition and images
        # latent_loss = tf.reduce_mean(tf.reduce_sum(tf.math.squared_difference(mean, encoded_images)))
        

        # Reconstruction loss
        reconstruction_loss = tf.reduce_mean(tf.reduce_sum(tf.math.squared_difference(x, reconstruction)))
        # reconstruction_loss_z = tf.reduce_mean(tf.math.squared_difference(x, reconstruction_z), axis=1)*0.01
        # KL divergence
        kl_loss = -0.5 * tf.reduce_sum(1 + logvar - tf.square(mean) - tf.exp(logvar), axis=1)
        # Latent loss between predictions of condition and images
        #latent_loss = tf.reduce_mean(tf.math.squared_difference(z, encoded_images), axis=1) #Might still be good
        
        #mmd_loss = self.compute_mmd(mean, conditions)
        
        #total_loss = tf.reduce_mean(reconstruction_loss + kl_loss) + mmd_loss
        total_loss = tf.reduce_mean(reconstruction_loss + kl_loss)
        #total_loss = tf.reduce_mean(reconstruction_loss + kl_loss)
        #total_loss = tf.reduce_mean(recon_loss + kl_loss)
        # add loss
        self.add_loss(total_loss)
        return reconstruction
    

class Generator(tf.keras.Model):
    def __init__(self, original_model):
        super(Generator, self).__init__()
        # self.conditions_encoder = original_model.conditions_encoder
        self.decoder    = original_model.decoder
        self.latent_dim = original_model.latent_dim
        self.nconditions = original_model.nconditions
        self.output_names = ['output']

    #     # Generate coordinates in latent space and concatenate with conditions 
    #     self.GenerateLatentSpace = tf.keras.Sequential([            
    #         tfkl.Concatenate(axis=-1)           
    #     ])

    #     #self.input_layer = tfkl.InputLayer(shape=(self.latent_dim,),name='input_layer')
    
    # def concatLatentSpace(self, conditions, z):
    #     #z = tf.random.normal(shape=(tf.shape(conditions)[0], self.latent_dim))
    #     condition_input_layer  = tfkl.InputLayer(shape=(self.nconditions,))     
    #     latent_input_layer     = tfkl.InputLayer(shape=(self.latent_dim,))     
    #     condition_input = condition_input_layer(conditions)
    #     latent_input = latent_input_layer(z)
    #     concat = tfkl.Concatenate(axis=-1)([condition_input, latent_input])


    #     return concat

    def call(self, conditions):
        # Concatenate conditions and z along the last axis
        z = tf.random.normal(shape=(tf.shape(conditions)[0], self.latent_dim))
        input = tf.concat([conditions,z], axis=-1)
        #concat = self.concatLatentSpace(conditions,z)
        reconstructed = self.decoder(input)
        return reconstructed
        
    # def reparameterize(self, mean, logvar):
    #     eps = tf.random.normal(shape=tf.shape(mean))
    #     return eps * tf.exp(logvar * .5) + mean
    

    # def reconstruct(self, conditions):      
    #     mean, logvar = tf.split(self.conditions_encoder(conditions), num_or_size_splits=2, axis=1)
    #     z = self.reparameterize(mean, logvar)
    #     decoded = self.decoder(z)
    #     return decoded

    # def call(self, conditions):
    #     reconstructed = self.reconstruct(conditions)
    #     return reconstructed