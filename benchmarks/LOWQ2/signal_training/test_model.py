import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf

# Define the number of plots
num_plots = 3
model_path = 'model.h5'

model = tf.keras.models.load_model(model_path)

#static input tensor
#input_tensor = np.array([0.5, 0.5])
input_tensors = [np.array([0.5, 0.5]), np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.5, 1.0])]

# Generate and plot the outputs
for j, input_tensor in enumerate(input_tensors):
    
    # Initialize the figure
    fig, axs = plt.subplots(num_plots, num_plots * 2, figsize=(40, 20))

    # Flatten the axes
    axs = axs.flatten()
    for ax in axs:
        # Plot data with a color scale from 0 to 4
        #im = ax.imshow(data, vmin=0, vmax=4)
        
        # Set x and y axes limits
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        
        # Add color bar
        #fig.colorbar(im, ax=ax, orientation='vertical')

    for i in range(num_plots * num_plots):
        # Predict the output for the input tensor
        output = model.predict(np.expand_dims(input_tensor, axis=0))

        print(output)
        print(output.shape)

        # Plot the output grid for the first channel
        im_charge = axs[i * 2].imshow(output[0,:,:,0], cmap='viridis', extent=[0, 10, 0, 10], vmin=0, vmax=2)
        axs[i * 2].set_title('Charge')
        fig.colorbar(im_charge, ax=axs[i * 2], orientation='vertical')

        # Plot the output grid for the second channel
        im_time = axs[i * 2 + 1].imshow(output[0,:,:,1], cmap='viridis', extent=[0, 10, 0, 10], vmin=0, vmax=2)
        axs[i * 2 + 1].set_title('Time')
        fig.colorbar(im_time, ax=axs[i * 2 + 1], orientation='vertical')


    # Show the plot
    plt.show()
    fig.savefig(f'output_plot_{j}.png')
