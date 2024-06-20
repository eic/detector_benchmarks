import matplotlib.pyplot as plt
import numpy as np
import uproot
from tensorflow import sparse, stack

# Load data from the ROOT file
file_path = 'output/Out_Convert_Big.root'
output_dir = 'plots/'

num_plots = 3
data_grid_size = 6
sensor_thickness = 300.0/55.0 # Thickness in pixel dimensions

# Assuming the ROOT file structure: MCParticles and PixelHits trees
infile = uproot.open(file_path)
tree  = infile['events']

#static input tensor
input_tensors = np.array([[0.5, 0.5, 0.0, 0.0],[0.25, 0.25, 0.0, 0.0],[0.0, 0.0,-0.05,-0.05],[0.5, 0.1,0.0,0.0],[0.0, 0.5,0.05,0.05],[0.25, 0.5,0.05,0.05]], dtype=np.float32)

# Range of input values to accept
input_range = np.array([0.05,0.05,0.01,0.01])

# Extracting data from the ROOT file
df = tree.arrays(['x', 'y', 'px', 'py', 'charge', 'time'], library="pd")

# Plot x, y 2d distribution in a 1x1 range
plt.figure()
plt.hist2d(df['x'], df['y'], bins=(400, 400), range=[[-3, 3], [-3, 3]], cmap='viridis')
plt.xlabel('x [pixel pitch]')
plt.ylabel('y [pixel pitch]')
plt.grid(True)
plt.savefig(output_dir + 'x_y_distribution.png')

# Plot px, py 2d distribution in a 1x1 range
plt.figure()
plt.hist2d(df['px'], df['py'], bins=(100, 100), range=[[-1, 1], [-1, 1]], cmap='viridis')
plt.xlabel('px')
plt.ylabel('py')
plt.savefig(output_dir + 'px_py_distribution.png')


#nHitsA = np.sum(np.where(df['time']>0, df['time'], 0), axis=1)
timeValues   = np.stack(df['time'].to_numpy())
chargeValues = np.stack(df['charge'].to_numpy())

nHits = np.sum(timeValues>0,axis=1)
# Plot the number of pixels with signal for each entry
plt.figure()
plt.hist(nHits, bins=12, range=(0, 12))
plt.xlabel('Number of hit pixels')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'num_hit.png')

# # Plot the x_pixel distribution
# plt.figure()
# plt.hist(np.concatenate(df['pixel_x'].values), bins=data_grid_size, range=(0, data_grid_size))
# plt.xlabel('x pixel')
# plt.ylabel('Number of entries')
# plt.savefig(output_dir + 'x_pixel_distribution.png')

# # Plot the y_pixel distribution
# plt.figure()
# plt.hist(np.concatenate(df['pixel_y'].values), bins=data_grid_size, range=(0, data_grid_size))
# plt.xlabel('y pixel')
# plt.ylabel('Number of entries')
# plt.savefig(output_dir + 'y_pixel_distribution.png')

# Plot the charge distribution
plt.figure()
plt.hist(chargeValues[np.where(chargeValues>0)] , bins=6, range=(0, 6))
plt.xlabel('Charge')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'charge_distribution.png')

# Plot the time distribution
plt.figure()
plt.hist(timeValues[np.where(timeValues>0)], bins=30, range=(0, 30))
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

    if len(df2) == 0:
        continue
    timeValues2  = np.stack(df2['time'].to_numpy())
    nHits2       = np.sum(timeValues2>0,axis=1)
    # Plot the length of pixel_x for each entry
    plt.figure()
    plt.hist(nHits2, bins=12, range=(0, 12))
    plt.xlabel('Number of hit pixels')
    plt.ylabel('Number of entries')
    plt.savefig(output_dir + 'num_hit_'+output_extension)

    df_sample = df2.head(n=num_plots*num_plots).reset_index()

    target_tensors = np.stack([df_sample['charge'].to_numpy(), df_sample['time'].to_numpy()], axis=2)
    print(target_tensors.shape)

    # Initialize the figure
    fig, axs = plt.subplots(num_plots, num_plots * 2, figsize=(40, 20))

    # Flatten the axes
    axs = axs.flatten()
    for ax in axs:
        # Plot data with a color scale from 0 to 4
        #im = ax.imshow(data, vmin=0, vmax=4)
        
        # Set x and y axes limits
        ax.set_xlim(-data_grid_size/2, data_grid_size/2)
        ax.set_ylim(-data_grid_size/2, data_grid_size/2)

    # Plot the data for the entries in the dataframe on a grid
    for i, tensor in enumerate(target_tensors):

        # Create empty 2D grids for charge and time
        charge_grid = np.reshape(tensor[:,0].astype(float),(data_grid_size, data_grid_size))
        time_grid   = np.reshape(tensor[:,1].astype(float),(data_grid_size, data_grid_size))

        axs[i * 2].arrow(df_sample['x'][i], df_sample['y'][i], df_sample['px'][i]*sensor_thickness, df_sample['py'][i]*sensor_thickness, head_width=0.01, head_length=0.01, fc='red', ec='red')
        axs[i * 2 + 1].arrow(df_sample['x'][i], df_sample['y'][i], df_sample['px'][i]*sensor_thickness, df_sample['py'][i]*sensor_thickness, head_width=0.01, head_length=0.01, fc='red', ec='red')

        # Plot the charge data
        im_charge = axs[i * 2].imshow(charge_grid, cmap='viridis', extent=[-data_grid_size/2, data_grid_size/2, -data_grid_size/2, data_grid_size/2], vmin=0, vmax=3)
        axs[i * 2].set_title('Charge')
        fig.colorbar(im_charge, ax=axs[i * 2], orientation='vertical')

        # Plot the time data
        im_time = axs[i * 2 + 1].imshow(time_grid, cmap='viridis', extent=[-data_grid_size/2, data_grid_size/2, -data_grid_size/2, data_grid_size/2], vmin=0, vmax=10)
        axs[i * 2 + 1].set_title('Time')
        fig.colorbar(im_time, ax=axs[i * 2 + 1], orientation='vertical')

    # Show the plot
    plt.show()

    # Build output name based on input tensor
    output_name = output_dir + 'tpx4_input_' + output_extension

    # Save the plot
    fig.savefig(output_name)