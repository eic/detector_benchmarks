import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import onnxruntime as ort
import uproot
import pandas as pd
import tensorflow as tf

output_dir = 'plots/'
model_base = "model_digitization"
model_name = model_base+"_latent.onnx"
# Load the ONNX model
sess = ort.InferenceSession(model_name)

condition_columns    = ['x', 'y', 'px', 'py']
condition_ranges     = [[0, 0.5], [0, 0.5], [-0.2, 0.2], [-0.2, 0.2]]
nConditions = len(condition_columns)

input_name = sess.get_inputs()[0].name

# Load data from the ROOT file
file_path = 'output/Out_Convert_Big.root'
output_dir = 'plots/'

# Assuming the ROOT file structure: MCParticles and PixelHits trees
infile = uproot.open(file_path)
tree  = infile['events']

# Extracting data from the ROOT file
df = tree.arrays(['x', 'y', 'px', 'py', 'pixel_x', 'pixel_y', 'charge', 'time'], entry_stop=100000)

data_grid_size = 6

target_tensors = np.concatenate([df['charge'].to_numpy().astype(np.float32), df['time'].to_numpy().astype(np.float32)], axis=1)
       
conditions_tensors = df[condition_columns].to_numpy()
conditions_tensors = np.array([list(t) for t in conditions_tensors]).astype(np.float32)
 
input_tensors = np.concatenate([conditions_tensors, target_tensors], axis=1)

# Predict the output for the input tensor
output = sess.run(None, {input_name: input_tensors})

nOutputs = output[0].shape[-1]

# Plot 2D scatter plots showing the correlation between the first 2 latent dimensions
plt.figure()
plt.scatter(output[0][:,0], output[0][:,1])
plt.xlim([-5, 5])
plt.ylim([-5, 5])
plt.xlabel('Latent dimension 1')
plt.ylabel('Latent dimension 2')
plt.savefig(output_dir + model_base + '_latent_space-01.png')

#Plot a histogram of each latent dimension on a grid
fig, axs = plt.subplots(int(nOutputs/5),5, figsize=(20, 10))
for i in range(nOutputs):
    row=int(i//5)
    col=int(i%5)
    axs[row,col].hist(output[0][:,i], bins=400, range=(-5,5))
    #axs[row,col].set_yscale('log')
    axs[row,col].set_title('Latent dimension '+str(i))



plt.savefig(output_dir + model_base + '_latent_histograms.png')

# Plot 2D scatter plots showing the correlation between the conditions and output dimensions
# A matrix of images

for i in range(nConditions):
    out_tag = model_base + '_latent_hist_'+condition_columns[i]
    nRows = 5
    fig, axs = plt.subplots(nRows, nOutputs//nRows, figsize=(20, 10))
    for j in range(nOutputs):
        col = int(j//nRows)
        row = int(j%nRows)
        axs[row,col].hist2d(conditions_tensors[:,i], output[0][:,j], bins=(200, 200), cmap=plt.cm.jet, range=[condition_ranges[i], [-5, 5]])#, norm=colors.LogNorm())
        axs[row,col].set_xlabel(condition_columns[i])
        axs[row,col].set_ylabel('Output '+str(j))
    plt.savefig(output_dir + out_tag + '.png')
