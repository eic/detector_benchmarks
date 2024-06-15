import matplotlib.pyplot as plt
import numpy as np
import uproot
from tensorflow import sparse, stack

# Load data from the ROOT file
file_path = 'output/Out_Convert_tpx4-6.root'
output_dir = 'plots/'

num_plots = 3
data_grid_size = 6

# Assuming the ROOT file structure: MCParticles and PixelHits trees
infile = uproot.open(file_path)
tree  = infile['events']

#static input tensor
input_tensors = np.array([[0.5, 0.5, 0.0, 0.0],[0.25, 0.25, 0.0, 0.0],[0.0, 0.0,-0.05,-0.05],[0.5, 0.1,0.0,0.0],[0.0, 0.5,0.05,0.05],[0.25, 0.5,0.05,0.05]], dtype=np.float32)

# Range of input values to accept
input_range = np.array([0.05,0.05,0.01,0.01])

# Extracting data from the ROOT file
df = tree.arrays(['x', 'y', 'px', 'py', 'pixel_x', 'pixel_y', 'charge', 'time'], library='pd')

# Plot x, y 2d distribution in a 1x1 range
plt.figure()
plt.hist2d(df['x'], df['y'], bins=(100, 100), range=[[0, 1], [0, 1]], cmap='viridis')
plt.xlabel('x [pixel pitch]')
plt.ylabel('y [pixel pitch]')
plt.savefig(output_dir + 'x_y_distribution.png')

# Plot px, py 2d distribution in a 1x1 range
plt.figure()
plt.hist2d(df['px'], df['py'], bins=(100, 100), range=[[-1, 1], [-1, 1]], cmap='viridis')
plt.xlabel('px')
plt.ylabel('py')
plt.savefig(output_dir + 'px_py_distribution.png')

# Plot the number of pixels with signal for each entry
plt.figure()
plt.hist(df['pixel_x'].apply(len), bins=12, range=(0, 12))
plt.xlabel('Number of hit pixels')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'num_hit.png')

# Plot the x_pixel distribution
plt.figure()
plt.hist(df['pixel_x'], bins=data_grid_size, range=(0, data_grid_size))
plt.xlabel('x pixel')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'x_pixel_distribution.png')

# Plot the y_pixel distribution
plt.figure()
plt.hist(df['pixel_y'], bins=data_grid_size, range=(0, data_grid_size))
plt.xlabel('y pixel')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'y_pixel_distribution.png')

# Plot the charge distribution
plt.figure()
plt.hist(df['charge'], bins=6, range=(0, 6))
plt.xlabel('Charge')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'charge_distribution.png')

# Plot the time distribution
plt.figure()
plt.hist(df['time'], bins=10, range=(0, 10))
plt.xlabel('Time')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'time_distribution.png')

for j, input_tensor in enumerate(input_tensors[:,0:4]):

    output_extension = 'x-' + str(input_tensor[0]) + '_y-' + str(input_tensor[1]) + '_px-' + str(input_tensor[2]) + '_py-' + str(input_tensor[3]) + '.png'

    print(input_tensor)
    print(len(df))
    # Filter the df by +/- the input range on x, y, px and py
    df2 = df[(df['x'] > input_tensor[0] - input_range[0]) & (df['x'] < input_tensor[0] + input_range[0])]
    print(len(df2))
    df2 = df2[(df2['y'] > input_tensor[1] - input_range[1]) & (df2['y'] < input_tensor[1] + input_range[1])]
    print(len(df2))
    df2 = df2[(df2['px'] > input_tensor[2] - input_range[2]) & (df2['px'] < input_tensor[2] + input_range[2])]
    print(len(df2))
    df2 = df2[(df2['py'] > input_tensor[3] - input_range[3]) & (df2['py'] < input_tensor[3] + input_range[3])]
    print(len(df2))
    print('')
    #print number of entries

    # Plot the length of pixel_x for each entry
    plt.figure()
    plt.hist(df2['pixel_x'].apply(len), bins=12, range=(0, 12))
    plt.xlabel('Number of hit pixels')
    plt.ylabel('Number of entries')
    plt.savefig(output_dir + 'num_hit_'+output_extension)

    df_sample = df2.head(n=num_plots*num_plots).reset_index()

    df_sample['pixel_x'] = df_sample['pixel_x'].apply(lambda x: np.array(x, dtype=object))
    df_sample['pixel_y'] = df_sample['pixel_y'].apply(lambda x: np.array(x, dtype=object))
    df_sample['charge'] = df_sample['charge'].apply(lambda x: np.array(x, dtype=object))
    df_sample['time'] = df_sample['time'].apply(lambda x: np.array(x, dtype=object))

    # Define a function to create a sparse tensor from a row
    def row_to_sparse_tensor(row):

        charge_indices = np.column_stack([row['pixel_x'], row['pixel_y'], np.zeros(len(row['pixel_x']))])
        time_indices = np.column_stack([row['pixel_x'], row['pixel_y'], np.ones(len(row['pixel_x']))])
        indices = np.concatenate([charge_indices, time_indices])
        charge_values = np.array([float(i) for i in row['charge']])
        time_values = np.array([float(i) for i in row['time']])
        values = np.concatenate([charge_values,time_values])
        dense_shape = [data_grid_size, data_grid_size, 2]
        sparse_tensor = sparse.reorder(sparse.SparseTensor(indices, values, dense_shape))
        return sparse.to_dense(sparse_tensor)

    # Apply the function to each row of the DataFrame
    #target_tensors = df.apply(row_to_sparse_tensor, axis=1)
    target_tensors = stack(df_sample.apply(row_to_sparse_tensor, axis=1))#.to_list())

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

    # Plot the data for the entries in the dataframe on a grid
    for i, tensor in enumerate(target_tensors):

        # Create empty 2D grids for charge and time
        charge_grid = np.zeros((data_grid_size, data_grid_size))
        time_grid = np.zeros((data_grid_size, data_grid_size))


        #print(tensor[:,:,0])
        # Plot the charge data
        im_charge = axs[i * 2].imshow(tensor[:,:,0], cmap='viridis', extent=[0, data_grid_size, 0, data_grid_size], vmin=0, vmax=3)
        axs[i * 2].set_title('Charge')
        fig.colorbar(im_charge, ax=axs[i * 2], orientation='vertical')

        # Plot the time data
        im_time = axs[i * 2 + 1].imshow(tensor[:,:,1], cmap='viridis', extent=[0, data_grid_size, 0, data_grid_size], vmin=0, vmax=10)
        axs[i * 2 + 1].set_title('Time')
        fig.colorbar(im_time, ax=axs[i * 2 + 1], orientation='vertical')

    # Show the plot
    plt.show()

    # Build output name based on input tensor
    output_name = output_dir + 'tpx4_input_' + output_extension

    # Save the plot
    fig.savefig(output_name)