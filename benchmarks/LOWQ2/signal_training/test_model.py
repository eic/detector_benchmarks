import matplotlib.pyplot as plt
import numpy as np
import onnxruntime as ort

output_dir = 'plots/'
model_base = "model_digitization"
model_name = model_base+".onnx"
# Load the ONNX model
sess = ort.InferenceSession(model_name)

# Define the number of plots
num_plots = 3
data_grid_size = 6

#static input tensor
input_tensors = np.array([[[0.5, 0.5, 0.0, 0.0]],[[0.25, 0.25, 0.0, 0.0]],[[0.0, 0.0,-0.05,-0.05]],[[0.5, 0.1,0.0,0.0]],[[0.0, 0.5,0.05,0.05]],[[0.25, 0.5,0.05,0.05]]], dtype=np.float32)
input_tags    = ['x', 'y', 'px', 'py']

input_name = sess.get_inputs()[0].name

# Generate and plot the outputs
for j, input_tensor in enumerate(input_tensors[:,:,0:4]):
    
    print(input_tensor) 

    # Initialize the figure
    fig, axs = plt.subplots(num_plots, num_plots * 2, figsize=(40, 20))

    # Flatten the axes
    axs = axs.flatten()
    for ax in axs:
        # Plot data with a color scale from 0 to 4
        #im = ax.imshow(data, vmin=0, vmax=4)
        
        # Set x and y axes limits
        ax.set_xlim(0, data_grid_size)
        ax.set_ylim(0, data_grid_size)
        
        # Add color bar
        #fig.colorbar(im, ax=ax, orientation='vertical')

    for i in range(num_plots * num_plots):
        # Predict the output for the input tensor
        output = sess.run(None, {input_name: input_tensor})
        output = output[0]
        output = output.reshape((1, 2, data_grid_size, data_grid_size))

        round_output = np.round(output)
        #round_output = output

        # Plot the output grid for the first channel
        im_charge = axs[i * 2].imshow(round_output[0,0,:,:], cmap='viridis', extent=[0, data_grid_size, 0, data_grid_size], vmin=0.1, vmax=3)
        axs[i * 2].set_title('Charge')
        fig.colorbar(im_charge, ax=axs[i * 2], orientation='vertical')

        # Plot the output grid for the second channel
        im_time = axs[i * 2 + 1].imshow(round_output[0,1,:,:], cmap='viridis', extent=[0, data_grid_size, 0, data_grid_size], vmin=0.1, vmax=10)
        axs[i * 2 + 1].set_title('Time')
        fig.colorbar(im_time, ax=axs[i * 2 + 1], orientation='vertical')

    # Set the output name based on the input tensor
    output_name = output_dir + model_base
    # loop over inner ddimension of input tensor
    for i in range(input_tensor.shape[-1]):
        output_name += '_' + input_tags[i] + '-' + str(input_tensor[0][i])
    #'_x-' + str(input_tensor[0][0]) + '_y-' + str(input_tensor[0][1]) + '_px-' + str(input_tensor[0][2]) + '_py-' + str(input_tensor[0][3]) + '.png'
    output_name += '.png'

    # Show the plot
    plt.show()
    fig.savefig(output_name)
