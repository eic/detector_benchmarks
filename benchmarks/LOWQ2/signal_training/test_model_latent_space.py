import matplotlib.pyplot as plt
import numpy as np
import onnxruntime as ort
import uproot
import pandas as pd
import tensorflow as tf

output_dir = 'plots/'
model_base = "model_tpx4_new3"
model_name = model_base+"_latent.onnx"
# Load the ONNX model
sess = ort.InferenceSession(model_name)

condition_columns    = ['x', 'y', 'px', 'py']
nConditions = len(condition_columns)

input_name = sess.get_inputs()[0].name

# Load data from the ROOT file
file_path = 'output/Out_Convert_tpx4-6.root'
output_dir = 'plots/'

# Assuming the ROOT file structure: MCParticles and PixelHits trees
infile = uproot.open(file_path)
tree  = infile['events']

# Extracting data from the ROOT file
df = tree.arrays(['x', 'y', 'px', 'py', 'pixel_x', 'pixel_y', 'charge', 'time'], library='pd')

data_grid_size = 6

input_data = df[condition_columns].values.astype(np.float32)

# Define a function to create a sparse tensor from a row
def row_to_sparse_tensor(row):
    charge_indices = np.column_stack([row['pixel_x'], row['pixel_y'], np.zeros(len(row['pixel_x']))])
    time_indices   = np.column_stack([row['pixel_x'], row['pixel_y'], np.ones(len(row['pixel_x']))])
    indices        = np.concatenate([charge_indices, time_indices])
    values         = np.concatenate([row['charge'], row['time']])
    dense_shape    = [data_grid_size, data_grid_size, 2]
    sparse_tensor  = tf.sparse.reorder(tf.sparse.SparseTensor(indices, values, dense_shape))
    return tf.sparse.to_dense(sparse_tensor)

# Apply the function to each row of the DataFrame
#target_tensors = df.apply(row_to_sparse_tensor, axis=1)
target_tensors = tf.stack(df.apply(row_to_sparse_tensor, axis=1).to_list())

# Reshape the target_tensors so that other than the event dimension, the other dimensions are the flat
target_tensors = tf.reshape(target_tensors, (target_tensors.shape[0], -1)).numpy()

# Create input tensors
conditions_tensors = df[condition_columns].to_numpy()
#conditions_tensors = df[['x', 'y']].to_numpy()
#conditions_tensors = df[[]].to_numpy()

# Concatenate the conditions_tensors and target_tensors along final axis
input_tensors = np.concatenate([conditions_tensors, target_tensors], axis=1)

# Predict the output for the input tensor
output = sess.run(None, {input_name: input_tensors})

nOutputs = output[0].shape[-1]

# Plot 2D scatter plots showing the correlation between the first 2 latent dimensions
plt.figure()
plt.scatter(output[0][:,0], output[0][:,1])
plt.xlabel('Latent dimension 1')
plt.ylabel('Latent dimension 2')
plt.savefig(output_dir + 'latent_space.png')



# Plot 2D scatter plots showing the correlation between the conditions and output dimensions
# A matrix of images

fig, axs = plt.subplots(nConditions,nOutputs , figsize=(20, 10))
for i in range(nConditions):
    for j in range(nOutputs):
        axs[i,j].scatter(input_data[:,i], output[0][:,j])
        axs[i,j].set_xlabel(condition_columns[i])
        axs[i,j].set_ylabel('Output '+str(j))
plt.savefig(output_dir + 'scatter.png')

