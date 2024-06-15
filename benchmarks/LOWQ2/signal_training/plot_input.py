import matplotlib.pyplot as plt
import numpy as np
import uproot
from tensorflow import sparse, stack

# Load data from the ROOT file
file_path = 'output/Out_Convert.root'

num_plots = 10

# Assuming the ROOT file structure: MCParticles and PixelHits trees
infile = uproot.open(file_path)
tree  = infile['events']

# Extracting data from the ROOT file
df = tree.arrays(['x', 'y', 'px', 'py', 'pixel_x', 'pixel_y', 'charge', 'time'], library='pd',entry_start=0, entry_stop=num_plots*num_plots)

# Define a function to create a sparse tensor from a row
def row_to_sparse_tensor(row):
    charge_indices = np.column_stack([row['pixel_x'], row['pixel_y'], np.zeros(len(row['pixel_x']))])
    time_indices = np.column_stack([row['pixel_x'], row['pixel_y'], np.ones(len(row['pixel_x']))])
    indices = np.concatenate([charge_indices, time_indices])
    values = np.concatenate([row['charge'], row['time']])
    dense_shape = [10, 10, 2]
    sparse_tensor = sparse.reorder(sparse.SparseTensor(indices, values, dense_shape))
    return sparse.to_dense(sparse_tensor)

# Apply the function to each row of the DataFrame
#target_tensors = df.apply(row_to_sparse_tensor, axis=1)
target_tensors = stack(df.apply(row_to_sparse_tensor, axis=1).to_list())

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

# Plot the data for the entries in the dataframe on a grid
for i, tensor in enumerate(target_tensors):

    # Create empty 2D grids for charge and time
    charge_grid = np.zeros((10, 10))
    time_grid = np.zeros((10, 10))


    #print(tensor[:,:,0])
    # Plot the charge data
    im_charge = axs[i * 2].imshow(tensor[:,:,0], cmap='viridis', extent=[0, 10, 0, 10], vmin=0, vmax=4)
    axs[i * 2].set_title('Charge')
    fig.colorbar(im_charge, ax=axs[i * 2], orientation='vertical')

    # Plot the time data
    im_time = axs[i * 2 + 1].imshow(tensor[:,:,1], cmap='viridis', extent=[0, 10, 0, 10], vmin=0, vmax=4)
    axs[i * 2 + 1].set_title('Time')
    fig.colorbar(im_time, ax=axs[i * 2 + 1], orientation='vertical')

# Show the plot
plt.show()

# Save the plot
fig.savefig('data_plot.png')