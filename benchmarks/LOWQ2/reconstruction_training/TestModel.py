import onnxruntime as ort
import argparse
import numpy as np
from ProcessData import create_arrays
import matplotlib.pyplot as plt

# Parse arguments
parser = argparse.ArgumentParser(description='Train a regression model for the Tagger.')
parser.add_argument('dataFiles', type=str, nargs='+', help='Path to the data files')
parser.add_argument('beamEnergy', type=float, help='Electron beam energy')
args = parser.parse_args()
dataFiles = args.dataFiles
beamEnergy = args.beamEnergy

input_data, target_data = create_arrays(dataFiles, beamEnergy)

target_data = np.array(target_data)

# Load the ONNX model
onnx_file_path = "regression_model.onnx"
session = ort.InferenceSession(onnx_file_path)

# Run the model on the input data
input_name = session.get_inputs()[0].name
output_name = session.get_outputs()[0].name
input_data = np.array(input_data,dtype=np.float32)
output = session.run([output_name], {input_name: input_data})
output = np.array(output[0])#*beamEnergy
print(output)

diff = (target_data - output)/target_data
diffrange = [[-2,2],[-2,2],[-2,2]]
datarange = [[-0.02,0.02],[-0.02,0.02],[-1,0]]


# Creates histograms to compare the target and output data
fig, axs = plt.subplots(3, 3, figsize=(6, 12))
for i in range(3):
    # 2D histograms showing trends in the data
    axs[0,i].hist2d(target_data[:,i], output[:,i], bins=100, range=[datarange[i],datarange[i]], cmap="viridis", label="Output vs Target")
    axs[0,i].set_xlabel(f"Variable {i} Target")
    axs[0,i].set_ylabel(f"Variable {i} Output")

    axs[1,i].hist(diff[:,i], bins=100, alpha=0.5, range=diffrange[i], label="Difference")
    axs[1,i].set_xlabel(f"Variable {i} Difference")
    axs[1,i].set_ylabel("Counts")

    axs[2,i].hist2d(target_data[:,i], diff[:,i], bins=100, range=[datarange[i],diffrange[i]], cmap="viridis", label="Difference vs Target")
    axs[2,i].set_xlabel(f"Variable {i} Target")
    axs[2,i].set_ylabel(f"Variable {i} Difference")


plt.show()
plt.savefig("output_vs_target.png")
