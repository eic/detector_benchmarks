import matplotlib.pyplot as plt
import numpy as np
import uproot
from tensorflow import sparse, stack
import onnxruntime as ort

output_dir = 'plots/'
model_base = "model_tpx4_new2"
model_name = model_base+".onnx"
# Load the ONNX model
sess = ort.InferenceSession(model_name)

input_name = sess.get_inputs()[0].name

# Load data from the ROOT file
file_path = 'output/Out_Convert_tpx4-6.root'
output_dir = 'plots/'

num_plots = 3
data_grid_size = 6

# Assuming the ROOT file structure: MCParticles and PixelHits trees
infile = uproot.open(file_path)
tree  = infile['events']

# Extracting data from the ROOT file
df = tree.arrays(['x', 'y', 'px', 'py', 'pixel_x', 'pixel_y', 'charge', 'time'], library='pd')

input_data = df[['x', 'y', 'px', 'py']].values.astype(np.float32)

# Predict the output for the input tensor
output = sess.run(None, {input_name: input_data})
output = output[0]
output = output.reshape((len(input_data), data_grid_size, data_grid_size, 2))

round_output = np.round(output)

# Calculate the number of pixels with hits > 0 for each entry
output_nhits = np.sum(round_output[:,:,:,0] > 0, axis=(1,2))
# Plot the number of pixels with a hit for predicted entry
plt.figure()
plt.hist(output_nhits, bins=12, range=(0, 12))
plt.xlabel('Number of hit pixels')
plt.ylabel('Number of entries')
plt.savefig(output_dir + 'num_hit_pred.png')