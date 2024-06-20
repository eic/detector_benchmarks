import matplotlib.pyplot as plt
import numpy as np
import uproot
from tensorflow import sparse, stack
import onnxruntime as ort

output_dir = 'plots/'
model_base = "model_digitization"
model_name = model_base+".onnx"
# Load the ONNX model
sess = ort.InferenceSession(model_name)

input_name = sess.get_inputs()[0].name

# Load data from the ROOT file
file_path = 'output/Out_Convert_Big.root'
output_dir = 'plots/'

num_plots = 3
data_grid_size = 6

# Assuming the ROOT file structure: MCParticles and PixelHits trees
infile = uproot.open(file_path)
tree  = infile['events']

# Extracting data from the ROOT file
df = tree.arrays(['x', 'y', 'px', 'py', 'pixel_x', 'pixel_y', 'charge', 'time'], library='pd')

input_data = df[['x', 'y', 'px', 'py']].values.astype(np.float32)
#input_data = df[['x', 'y']].values.astype(np.float32)

# Predict the output for the input tensor
output = sess.run(None, {input_name: input_data})
output = output[0]
output = output.reshape((len(input_data), 2, data_grid_size, data_grid_size))

round_output = np.round(output)
round_output = np.transpose(round_output, (0, 2, 3, 1))

# Calculate the number of pixels with hits > 0 for each entry
output_nhits = np.sum(round_output[:,:,:,0] > 0, axis=(1,2))
# Plot the number of pixels with a hit for predicted entry
plt.figure()
plt.hist(output_nhits, bins=12, range=(0, 12))
plt.xlabel('Number of hit pixels')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'num_hit_pred.png')

# Plot the charge distribution
plt.figure()
plt.hist(round_output[round_output[:,:,:,0] > 0][...,0], bins=6, range=(0, 6))
plt.xlabel('Charge')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'charge_distribution_pred.png')

# Plot the time distribution
plt.figure()
plt.hist(round_output[round_output[:,:,:,0] > 0][...,1], bins=30, range=(0, 30))
plt.xlabel('Time')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'time_distribution_pred.png')

input_tensors = np.array([[0.5, 0.5, 0.0, 0.0],[0.25, 0.25, 0.0, 0.0],[0.0, 0.0,-0.05,-0.05],[0.5, 0.1,0.0,0.0],[0.0, 0.5,0.05,0.05],[0.25, 0.5,0.05,0.05]], dtype=np.float32)

input_range = np.array([0.05,0.05,0.01,0.01])


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

    input_indeces = df2.index.values

    # Plot the length of pixel_x for each entry
    plt.figure()
    plt.hist(np.sum(round_output[input_indeces][:,:,:,0] > 0, axis=(1,2)), bins=12, range=(0, 12))
    plt.xlabel('Number of hit pixels')
    plt.ylabel('Number of entries')
    plt.savefig(output_dir + 'num_hit_pred_'+output_extension)
    #print number of entries