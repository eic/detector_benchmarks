import onnxruntime as ort
import argparse
import numpy as np
from ProcessData import create_arrays
import matplotlib.pyplot as plt

# Parse arguments
parser = argparse.ArgumentParser(description='Train a regression model for the Tagger.')
parser.add_argument('--modelFile', type=str, default="regression_model.onnx", help='Path to the ONNX model file')
parser.add_argument('--dataFiles', type=str, nargs='+', help='Path to the data files')
parser.add_argument('--beamEnergy', type=float, help='Electron beam energy')
parser.add_argument('--outGraphFile', type=str, default="output_vs_target.png", help='Output file for the graph')
parser.add_argument('--outGraphFile2', type=str, default="output_vs_target2.png", help='Output file for the graph')
args = parser.parse_args()
modelFile     = args.modelFile
dataFiles     = args.dataFiles
beamEnergy    = args.beamEnergy
outGraphFile  = args.outGraphFile
outGraphFile2 = args.outGraphFile2

input_data, target_data = create_arrays(dataFiles, beamEnergy)

target_data = np.array(target_data)

# Load the ONNX model
session = ort.InferenceSession(modelFile)

# Run the model on the input data
input_name = session.get_inputs()[0].name
output_name = session.get_outputs()[0].name
input_data = np.array(input_data,dtype=np.float32)
output = session.run([output_name], {input_name: input_data})
output = np.array(output[0])#*beamEnergy
#print(output)

out_theta = np.arctan2(np.sqrt(output[:,0]**2 + output[:,1]**2),output[:,2])
out_phi = np.arctan2(output[:,1],output[:,0])
out_mag = np.sqrt(output[:,0]**2 + output[:,1]**2 + output[:,2]**2)
in_theta = np.arctan2(np.sqrt(target_data[:,0]**2 + target_data[:,1]**2),target_data[:,2])
in_phi = np.arctan2(target_data[:,1],target_data[:,0])
in_mag = np.sqrt(target_data[:,0]**2 + target_data[:,1]**2 + target_data[:,2]**2)


thetadiff = out_theta - in_theta
phidiff = out_phi - in_phi
# Move phidiff to within -pi/2 and pi/2
phidiff = (phidiff + np.pi) % (2 * np.pi) - np.pi
magdiff = out_mag - in_mag

diff = (target_data - output)/target_data
diffrange = [[-5,5],[-5,5],[-2,2]]
datarange = [[-0.02,0.02],[-0.02,0.02],[-1,0]]


# Creates histograms to compare the target and output data
fig, axs = plt.subplots(3, 3, figsize=(12, 12))
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
plt.savefig(outGraphFile)

# Creates histograms to compare theta, phi and mag target and output data
fig2, axs2 = plt.subplots(3, 3, figsize=(12, 12))

thetarange = [np.pi-0.01,np.pi]
phirange = [-np.pi,np.pi]
magrange = [0,1]

thetadiffrange = [-0.02,0.02]
phidiffrange = [-np.pi,np.pi]
magdiffrange = [-0.2,0.2]

# 2D histograms showing trends in the data
axs2[0,0].hist2d(out_theta, in_theta, bins=100, range=[thetarange,thetarange], cmap="viridis", label="Output vs Target")
axs2[0,0].set_xlabel("Theta Target")
axs2[0,0].set_ylabel("Theta Output")

axs2[0,1].hist2d(out_phi, in_phi, bins=100, range=[phirange,phirange], cmap="viridis", label="Output vs Target")
axs2[0,1].set_xlabel("Phi Target")
axs2[0,1].set_ylabel("Phi Output")

axs2[0,2].hist2d(out_mag, in_mag, bins=100, range=[magrange,magrange], cmap="viridis", label="Output vs Target")
axs2[0,2].set_xlabel("Mag Target")
axs2[0,2].set_ylabel("Mag Output")

axs2[1,0].hist(thetadiff, bins=100, alpha=0.5, range=thetadiffrange, label="Difference")
axs2[1,0].set_xlabel("Theta Difference")
axs2[1,0].set_ylabel("Counts")

axs2[1,1].hist(phidiff, bins=100, alpha=0.5, range=phidiffrange, label="Difference")
axs2[1,1].set_xlabel("Phi Difference")
axs2[1,1].set_ylabel("Counts")

axs2[1,2].hist(magdiff, bins=100, alpha=0.5, range=magdiffrange, label="Difference")
axs2[1,2].set_xlabel("Mag Difference")
axs2[1,2].set_ylabel("Counts")

axs2[2,0].hist2d(in_theta, thetadiff, bins=100, range=[thetarange,thetadiffrange], cmap="viridis", label="Difference vs Target")
axs2[2,0].set_xlabel("Theta Target")
axs2[2,0].set_ylabel("Theta Difference")

axs2[2,1].hist2d(in_phi, phidiff, bins=100, range=[phirange,phidiffrange], cmap="viridis", label="Difference vs Target")
axs2[2,1].set_xlabel("Phi Target")
axs2[2,1].set_ylabel("Phi Difference")

axs2[2,2].hist2d(in_mag, magdiff, bins=100, range=[magrange,magdiffrange], cmap="viridis", label="Difference vs Target")
axs2[2,2].set_xlabel("Mag Target")
axs2[2,2].set_ylabel("Mag Difference")

plt.show()
plt.savefig(outGraphFile2)