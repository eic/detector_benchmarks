import onnxruntime as ort
import argparse
import numpy as np
from ProcessData import create_arrays
from RegressionModel import ProjectToX0Plane
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import norm
from scipy.optimize import curve_fit

# Parse arguments
parser = argparse.ArgumentParser(description='Train a regression model for the Tagger.')
parser.add_argument('--modelFile', type=str, default="regression_model.onnx", help='Path to the ONNX model file')
parser.add_argument('--dataFiles', type=str, nargs='+', help='Path to the data files')
parser.add_argument('--outDir', type=str, default=".", help='Output directory')
args = parser.parse_args()
modelFile     = args.modelFile
dataFiles     = args.dataFiles
outDir        = args.outDir
outGraphFile  = outDir + "/output_vs_target.png"
outGraphFile2 = outDir + "/transformed_output_vs_target.png"
outGraphFile3 = outDir + "/transformed_cut_output_vs_target.png"
outGraphFile4 = outDir + "/projected_output_vs_target.png"
correlationFile = outDir + "/correlations.png"
projectedCorrelationFile = outDir + "/projected_correlations.png"
differenceCorrelationFile = outDir + "/difference_correlations.png"

input_data, target_data = create_arrays(dataFiles)

target_data = np.array(target_data)

# Load the ONNX model
session = ort.InferenceSession(modelFile)

# Run the model on the input data
input_name = session.get_inputs()[0].name
output_name = session.get_outputs()[0].name
input_data = np.array(input_data,dtype=np.float32)
output = session.run([output_name], {input_name: input_data})
output = np.array(output[0])

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
diffrange = [[-5,5],[-5,5],[-0.5,0.5]]
datarange = [[-0.02,0.02],[-0.02,0.02],[-1,0]]

# Use the 'seismic' colormap
cmap = plt.get_cmap('seismic')

# Creates histograms to compare the target and output data
fig, axs = plt.subplots(3, 3, figsize=(12, 12))
for i in range(3):
    # 2D histograms showing trends in the data
    axs[0,i].hist2d(target_data[:,i], output[:,i], bins=100, range=[datarange[i],datarange[i]], cmap="seismic", norm=LogNorm(), label="Output vs Target")
    axs[0,i].set_xlabel(f"Variable {i} Target")
    axs[0,i].set_ylabel(f"Variable {i} Output")

    axs[1,i].hist(diff[:,i], bins=100, alpha=0.5, range=diffrange[i], label="Difference")
    axs[1,i].set_xlabel(f"Variable {i} Difference")
    axs[1,i].set_ylabel("Counts")

    axs[2,i].hist2d(target_data[:,i], diff[:,i], bins=100, range=[datarange[i],diffrange[i]], cmap="seismic", norm=LogNorm(), label="Difference vs Target")
    axs[2,i].set_xlabel(f"Variable {i} Target")
    axs[2,i].set_ylabel(f"Variable {i} Difference")

plt.show()
plt.savefig(outGraphFile)

# Creates histograms to compare theta, phi and mag target and output data
fig2, axs2 = plt.subplots(3, 3, figsize=(12, 12))

thetarange = [np.pi-0.01,np.pi]
phirange = [-np.pi,np.pi]
magrange = [0,1]

thetadiffrange = [-0.01,0.01]
phidiffrange = [-np.pi,np.pi]
magdiffrange = [-0.1,0.1]

# 2D histograms showing trends in the data
axs2[0,0].hist2d(out_theta, in_theta, bins=100, range=[thetarange,thetarange], cmap="seismic", norm=LogNorm(), label="Output vs Target")
axs2[0,0].set_xlabel("Theta Target")
axs2[0,0].set_ylabel("Theta Output")

axs2[0,1].hist2d(out_phi, in_phi, bins=100, range=[phirange,phirange], cmap="seismic", norm=LogNorm(), label="Output vs Target")
axs2[0,1].set_xlabel("Phi Target")
axs2[0,1].set_ylabel("Phi Output")

axs2[0,2].hist2d(out_mag, in_mag, bins=100, range=[magrange,magrange], cmap="seismic", norm=LogNorm(), label="Output vs Target")
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

axs2[2,0].hist2d(in_theta, thetadiff, bins=100, range=[thetarange,thetadiffrange], cmap="seismic", norm=LogNorm(), label="Difference vs Target")
axs2[2,0].set_xlabel("Theta Target")
axs2[2,0].set_ylabel("Theta Difference")

axs2[2,1].hist2d(in_phi, phidiff, bins=100, range=[phirange,phidiffrange], cmap="seismic", norm=LogNorm(), label="Difference vs Target")
axs2[2,1].set_xlabel("Phi Target")
axs2[2,1].set_ylabel("Phi Difference")

axs2[2,2].hist2d(in_mag, magdiff, bins=100, range=[magrange,magdiffrange], cmap="seismic", norm=LogNorm(), label="Difference vs Target")
axs2[2,2].set_xlabel("Mag Target")
axs2[2,2].set_ylabel("Mag Difference")

plt.show()
plt.savefig(outGraphFile2)

# Create histograms where the theta value has been cut at less than 3.14
fig3, axs3 = plt.subplots(3, 3, figsize=(12, 12))

out_theta_cut = out_theta[out_theta < 3.14]
in_theta_cut = in_theta[out_theta < 3.14]
thetadiff_cut = thetadiff[out_theta < 3.14]

out_phi_cut = out_phi[out_theta < 3.14]
in_phi_cut = in_phi[out_theta < 3.14]
phidiff_cut = phidiff[out_theta < 3.14]

out_mag_cut = out_mag[out_theta < 3.14]
in_mag_cut = in_mag[out_theta < 3.14]
magdiff_cut = magdiff[out_theta < 3.14]

axs3[0,0].hist2d(out_theta_cut, in_theta_cut, bins=100, range=[thetarange,thetarange], cmap="seismic", norm=LogNorm(), label="Output vs Target")
axs3[0,0].set_xlabel("Theta Target")
axs3[0,0].set_ylabel("Theta Output")

axs3[0,1].hist2d(out_phi_cut, in_phi_cut, bins=100, range=[phirange,phirange], cmap="seismic", norm=LogNorm(), label="Output vs Target")
axs3[0,1].set_xlabel("Phi Target")
axs3[0,1].set_ylabel("Phi Output")

axs3[0,2].hist2d(out_mag_cut, in_mag_cut, bins=100, range=[magrange,magrange], cmap="seismic", norm=LogNorm(), label="Output vs Target")
axs3[0,2].set_xlabel("Mag Target")
axs3[0,2].set_ylabel("Mag Output")

axs3[1,0].hist(thetadiff_cut, bins=100, alpha=0.5, range=thetadiffrange, label="Difference")
axs3[1,0].set_xlabel("Theta Difference")
axs3[1,0].set_ylabel("Counts")

axs3[1,1].hist(phidiff_cut, bins=100, alpha=0.5, range=phidiffrange, label="Difference")
axs3[1,1].set_xlabel("Phi Difference")
axs3[1,1].set_ylabel("Counts")

axs3[1,2].hist(magdiff_cut, bins=100, alpha=0.5, range=magdiffrange, label="Difference")
axs3[1,2].set_xlabel("Mag Difference")
axs3[1,2].set_ylabel("Counts")

axs3[2,0].hist2d(in_theta_cut, thetadiff_cut, bins=100, range=[thetarange,thetadiffrange], cmap="seismic", norm=LogNorm(), label="Difference vs Target")    
axs3[2,0].set_xlabel("Theta Target")
axs3[2,0].set_ylabel("Theta Difference")

axs3[2,1].hist2d(in_phi_cut, phidiff_cut, bins=100, range=[phirange,phidiffrange], cmap="seismic", norm=LogNorm(), label="Difference vs Target")
axs3[2,1].set_xlabel("Phi Target")
axs3[2,1].set_ylabel("Phi Difference")

axs3[2,2].hist2d(in_mag_cut, magdiff_cut, bins=100, range=[magrange,magdiffrange], cmap="seismic", norm=LogNorm(), label="Difference vs Target")
axs3[2,2].set_xlabel("Mag Target")
axs3[2,2].set_ylabel("Mag Difference")

plt.show()
plt.savefig(outGraphFile3)

# Create plots where a Gaussian has been fitted to the data
# Function to fit a Gaussian and plot the results
def plot_gaussian_fit(ax, data, range, xlabel, ylabel):
    def gaussian(x, mu, sigma, amplitude):
        return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    hist, bin_edges = np.histogram(data, bins=100, range=range, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    popt, _ = curve_fit(gaussian, bin_centers, hist, p0=[0, 0.01, 1])
    mu, sigma, amplitude = popt
    print(f"mu={mu}, sigma={sigma}, amplitude={amplitude}")

    x = np.linspace(range[0], range[1], 100)
    ax.plot(x, gaussian(x, *popt), 'k', linewidth=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.hist(data, bins=100, alpha=0.5, range=range, edgecolor='black', density=True)
    ax.legend([f'Fit: $\mu$={mu:.5f}, $\sigma$={sigma:.5f}'])

# Create histograms with Gaussian fits
fig4, axs4 = plt.subplots(3, 1, figsize=(8, 12))

plot_gaussian_fit(axs4[0], thetadiff, thetadiffrange, "Theta Difference", "Density")
plot_gaussian_fit(axs4[1], phidiff_cut, phidiffrange, "Phi Difference",   "Density")
plot_gaussian_fit(axs4[2], magdiff, magdiffrange,     "Mag Difference",   "Density")

plt.show()
plt.savefig(outGraphFile4)


# Look at the correlations between all of the variables in target_data and projected_inputs
# Create a comparison variable array
comparisson_variables = np.concatenate((input_data, target_data), axis=1)

# Project inputs onto the X0 plane
projected_inputs = ProjectToX0Plane().project_numpy(input_data)

# Concatenate the projected inputs and target data
projected_comparisson_variables = np.concatenate((projected_inputs, target_data), axis=1)
diff_cut = diff[(abs(projected_comparisson_variables[:, 3]) < 0.02)]
projected_comparisson_variables_cut = projected_comparisson_variables[(abs(projected_comparisson_variables[:, 3]) < 0.02)]
diff_cut = diff_cut[(abs(projected_comparisson_variables_cut[:, 2]+0.025) < 0.028)]
projected_comparisson_variables_cut = projected_comparisson_variables_cut[(abs(projected_comparisson_variables_cut[:, 2]+0.025) < 0.028)]  # Filter for px < 0.1

# Calculate limits for each variable based on the data
limits = {
    "ox": [np.min(comparisson_variables[:, 0]), np.max(comparisson_variables[:, 0])],
    "oy": [np.min(comparisson_variables[:, 1]), np.max(comparisson_variables[:, 1])],
    "y": [-200, 200],
    "z": [-17000,-9000],
    "px": [np.min(projected_comparisson_variables_cut[:, 2]), np.max(projected_comparisson_variables_cut[:, 2])],
    "py": [np.min(projected_comparisson_variables_cut[:, 3]), np.max(projected_comparisson_variables_cut[:, 3])],
    "opx": [np.min(comparisson_variables[:, 3]), np.max(comparisson_variables[:, 3])],
    "opy": [np.min(comparisson_variables[:, 4]), np.max(comparisson_variables[:, 4])],
    "opz": [np.min(comparisson_variables[:, 5]), np.max(comparisson_variables[:, 5])],
    "Px": [np.min(projected_comparisson_variables_cut[:, 4]), np.max(projected_comparisson_variables_cut[:, 4])],
    "Py": [np.min(projected_comparisson_variables_cut[:, 5]), np.max(projected_comparisson_variables_cut[:, 5])],
    "Pz": [np.min(projected_comparisson_variables_cut[:, 6]), np.max(projected_comparisson_variables_cut[:, 6])],
}


labels = ["ox","oy", "z", "opx", "opy", "opz", "Px", "Py", "Pz"]
fig5, axs5 = plt.subplots(9, 9, figsize=(30, 30))
for i in range(9):
    for j in range(9):
        if i == j:
            axs5[j, i].hist(comparisson_variables[:, i], range=limits[labels[i]], bins=100, alpha=0.5, label=labels[i])
            axs5[j, i].set_xlabel(labels[i])
            axs5[j, i].set_ylabel("Counts")
            #set log scale for y-axis if the data is skewed
            axs5[j, i].set_yscale('log')
        else:
            axs5[j, i].hist2d(comparisson_variables[:, i], comparisson_variables[:, j], range=[limits[labels[i]],limits[labels[j]]], bins=100, cmap="seismic", norm=LogNorm())
            axs5[j, i].set_xlabel(labels[i])
            axs5[j, i].set_ylabel(labels[j])
plt.tight_layout()
plt.savefig(correlationFile)
plt.show()



# Plot the correlations between all of the variables in target_data and projected_inputs
projected_labels = ["y", "z", "px", "py", "Px", "Py", "Pz"]
fig6, axs6 = plt.subplots(7, 7, figsize=(30, 30))
for i in range(7):
    for j in range(7):
        if i == j:
            axs6[j, i].hist(projected_comparisson_variables_cut[:, i], range=limits[projected_labels[i]], bins=100, alpha=0.5, label=projected_labels[i])
            axs6[j, i].set_xlabel(projected_labels[i])
            axs6[j, i].set_ylabel("Counts")
            #set log scale for y-axis if the data is skewed
            axs6[j, i].set_yscale('log')
        else:
            axs6[j, i].hist2d(projected_comparisson_variables_cut[:, i], projected_comparisson_variables_cut[:, j], range=[limits[projected_labels[i]],limits[projected_labels[j]]], bins=100, cmap="seismic", norm=LogNorm())
            axs6[j, i].set_xlabel(projected_labels[i])
            axs6[j, i].set_ylabel(projected_labels[j])


plt.tight_layout()
plt.savefig(projectedCorrelationFile)
plt.show()

# Plot the correlations between the output diferences and projected inputs
output_labels = ["pred_PX", "pred_PY", "pred_PZ"]
fig7, axs7 = plt.subplots(3, 7, figsize=(15, 8))
for i in range(3):
    for j in range(7):
        axs7[i, j].hist2d(projected_comparisson_variables_cut[:, j], diff_cut[:, i], range=[limits[projected_labels[j]],diffrange[i]], bins=100, cmap="seismic", norm=LogNorm())
        axs7[i, j].set_xlabel(projected_labels[j])
        axs7[i, j].set_ylabel(output_labels[i])


plt.tight_layout()
plt.savefig(differenceCorrelationFile)
plt.show()
