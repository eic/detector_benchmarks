import tensorflow as tf
#import tensorflow_probability as tfp

tfkl = tf.keras.layers
#tfpl = tfp.layers
#tfd = tfp.distributions
        
######################################################
# Define the Encoder
######################################################
class Encoder(tf.keras.Model):
    def __init__(self, nconditions, latent_dim, flat_shape=200):
        super(Encoder, self).__init__()
        self.normalizer = tfkl.Normalization(name='encode_normalizer')
        self.encode = tf.keras.Sequential([
            tfkl.InputLayer(shape=(flat_shape+nconditions,)),
            self.normalizer,
            tfkl.Dense(1024, activation='relu'),
            tfkl.Dense(256, activation='relu'),
            tfkl.Dense(2*latent_dim),
        ])

    def adapt(self, data):
        self.normalizer.adapt(data)

    def call(self, inputs):

        encoded = self.encode(inputs)
        mean, logvar = tf.split(encoded, num_or_size_splits=2, axis=1)
        return mean, logvar
    
######################################################
# Define the Decoder
######################################################
class Decoder(tf.keras.Model):
    def __init__(self, nconditions, latent_dim, flat_shape=200):
        super(Decoder, self).__init__()
        self.normalizer = tfkl.Normalization(name='decode_normalizer')
        self.decode = tf.keras.Sequential([            
            tfkl.InputLayer(shape=(latent_dim+nconditions,),name='input_layer'),
            tfkl.Dense(256, activation='relu'),
            tfkl.Dense(1024, activation='relu'),
            tfkl.Dense(flat_shape, name='output_layer')
        ])

    def adapt(self, data):
        self.normalizer.adapt(data)

    def call(self, conditions, z):
        normalized_conditions = self.normalizer(conditions)
        inputs = tf.concat([normalized_conditions, z], axis=-1)
        return self.decode(inputs)
    
######################################################
# Define the Discriminator
######################################################
class Adversarial(tf.keras.Model):
    def __init__(self, latent_dim, nconditions):
        super(Adversarial, self).__init__()
        self.reconstruct_conditions = tf.keras.Sequential([
            tfkl.InputLayer(shape=(latent_dim,)),
            tfkl.Dense(1024, activation='relu'),
            tfkl.Dense(256, activation='relu'),
            tfkl.Dense(nconditions, activation='linear')
        ])

    def call(self, inputs):
        x = self.reconstruct_conditions(inputs)
        return x

######################################################
# Define the Conditional, Adviserial Variational Autoencoder
######################################################
class VAE(tf.keras.Model):
    def __init__(self, latent_dim=50, nconditions=4, grid_size=6 ):
        super(VAE, self).__init__()
        self.flat_shape  = grid_size*grid_size*2
        self.latent_dim  = latent_dim
        self.nconditions = nconditions        

        self.encoder     = Encoder(nconditions, latent_dim, self.flat_shape)
        self.decoder     = Decoder(nconditions, latent_dim, self.flat_shape)
        self.adversarial = Adversarial(latent_dim, nconditions)
        
    def compile(self, r_optimizer, a_optimizer):
        super().compile()
        self.optimizer_encoder_decoder = r_optimizer
        self.optimizer_adversarial     = a_optimizer

    def reparameterize(self, mean, logvar):
        eps = tf.random.normal(shape=tf.shape(mean))
        return eps * tf.exp(logvar * .5) + mean
    
    def adapt(self, data):
        conditions = data[:,0:self.nconditions]
        self.encoder.adapt(data)
        self.decoder.adapt(conditions)
    
    def call(self, inputs):   
        conditions = inputs[:,0:self.nconditions]
        x = inputs[:,self.nconditions:]
        
        mean, logvar = self.encoder(inputs)
        #z = self.reparameterize(mean, logvar)
        reconstruction = self.decoder(conditions,mean)
        #advisarial_conditions = self.adversarial(mean)


        # Reconstruction loss
        #reconstruction_loss = tf.reduce_mean(tf.reduce_sum(tf.math.squared_difference(x, reconstruction)))
        #kl_loss = -0.5 * tf.reduce_sum(1 + logvar - tf.square(mean) - tf.exp(logvar), axis=1)
        #adversarial_loss = tf.reduce_mean(tf.reduce_sum(tf.math.squared_difference(conditions, advisarial_conditions)))
        
        #total_loss = tf.reduce_mean(reconstruction_loss + kl_loss - adversarial_loss)
        
        # add loss
        #self.add_loss(total_loss)
        return reconstruction
    
    def train_step(self, input):
        data, targets = input  # Unpack the inputs from the data tuple
        conditions = data[:,0:self.nconditions]
        x = data[:,self.nconditions:]

        with tf.GradientTape() as tape:
            # Forward pass
            mean, logvar = self.encoder(data)
            z = self.reparameterize(mean, logvar)
            reconstruction = self.decoder(conditions, z)            
            adversarial_conditions = self.adversarial(mean)
            # Compute losses
            reconstruction_loss = tf.reduce_sum(tf.math.squared_difference(x, reconstruction), axis=-1)
            kl_loss = -0.5 * tf.reduce_sum(1 + logvar - tf.square(mean) - tf.exp(logvar), axis=-1)
            adversarial_lossA = tf.reduce_sum(tf.math.squared_difference(conditions, adversarial_conditions), axis=-1)

            # Compute total loss
            total_loss = tf.reduce_mean(reconstruction_loss + kl_loss - adversarial_lossA)

        # Compute gradients with respect to the weights of the encoder and decoder
        grads_encoder_decoder = tape.gradient(total_loss, self.encoder.trainable_variables + self.decoder.trainable_variables)

        # Apply gradients
        self.optimizer_encoder_decoder.apply_gradients(zip(grads_encoder_decoder, self.encoder.trainable_variables + self.decoder.trainable_variables))

        with tf.GradientTape() as tape:
            # Compute adversarial loss
            adversarial_conditions = self.adversarial(mean)
            adversarial_loss = tf.reduce_mean(tf.reduce_sum(tf.math.squared_difference(conditions, adversarial_conditions), axis=-1))

        # Compute gradients with respect to the weights of the adversarial network
        grads_adversarial = tape.gradient(adversarial_loss, self.adversarial.trainable_variables)

        # Apply gradients
        self.optimizer_adversarial.apply_gradients(zip(grads_adversarial, self.adversarial.trainable_variables))

        return {"loss": total_loss, "reconstruction_loss": reconstruction_loss, "kl_loss": kl_loss, "adversarial_loss": adversarial_loss}
    
    def test_step(self, input):
        data, targets = input
        conditions = data[:,0:self.nconditions]
        x = data[:,self.nconditions:]
        
        mean, logvar = self.encoder(data)
        z = self.reparameterize(mean, logvar)
        reconstruction = self.decoder(conditions,z)
        advisarial_conditions = self.adversarial(mean)

        # Reconstruction loss
        reconstruction_loss = tf.reduce_sum(tf.math.squared_difference(x, reconstruction), axis=-1)
        kl_loss = -0.5 * tf.reduce_sum(1 + logvar - tf.square(mean) - tf.exp(logvar), axis=-1)
        adversarial_loss = tf.reduce_sum(tf.math.squared_difference(conditions, advisarial_conditions), axis=-1)
        
        total_loss = tf.reduce_mean(reconstruction_loss + kl_loss - adversarial_loss)
        
        # Return a dictionary containing the loss and the metrics
        return {"loss": total_loss, "reconstruction_loss": reconstruction_loss, "kl_loss": kl_loss, "adversarial_loss": adversarial_loss}
        #return reconstruction


######################################################
# Define the Generator
######################################################
class Generator(tf.keras.Model):
    def __init__(self, original_model):
        super(Generator, self).__init__()
        # self.conditions_encoder = original_model.conditions_encoder
        self.decoder    = original_model.decoder
        self.latent_dim = original_model.latent_dim
        self.nconditions = original_model.nconditions
        self.output_names = ['output']

    def call(self, conditions):
        # Concatenate conditions and z along the last axis
        z = tf.random.normal(shape=(tf.shape(conditions)[0], self.latent_dim))
        #input = tf.concat([conditions,z], axis=-1)
        #concat = self.concatLatentSpace(conditions,z)
        reconstructed = self.decoder(conditions, z)
        return reconstructed
    
######################################################
# Define latent space encoder
######################################################
class LatentSpace(tf.keras.Model):
    def __init__(self, original_model):
        super(LatentSpace, self).__init__()
        self.flat_shape  = original_model.flat_shape
        self.nconditions = original_model.nconditions
        self.input_layer = tfkl.InputLayer(input_shape=(self.flat_shape+self.nconditions,))
        self.encoder    = original_model.encoder
        self.output_names = ['output']

    def call(self, input):
        mean, logvar = self.encoder(input)
        return mean, logvar