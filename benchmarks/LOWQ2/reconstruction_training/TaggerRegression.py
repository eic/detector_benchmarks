# import edm4hep
import torch
import argparse
from ProcessData import create_arrays

from RegressionModel import makeModel, trainModel

# Parse arguments
parser = argparse.ArgumentParser(description='Train a regression model for the Tagger.')
parser.add_argument('dataFiles', type=str, nargs='+', help='Path to the data files')
parser.add_argument('beamEnergy', type=float, help='Electron beam energy')
args = parser.parse_args()
dataFiles = args.dataFiles
beamEnergy = args.beamEnergy

input_data, target_data = create_arrays(dataFiles, beamEnergy)
print(input_data)
print(target_data)

torch_input_data = torch.tensor(input_data)
torch_target_data = torch.tensor(target_data)

# Split data into training and validation sets
validation_fraction = 0.1
split_index = int(len(torch_input_data) * (1 - validation_fraction))
val_input_data = torch_input_data[split_index:]
val_target_data = torch_target_data[split_index:]
torch_input_data = torch_input_data[:split_index]
torch_target_data = torch_target_data[:split_index] 

epochs = 1000
model = trainModel(epochs, torch_input_data, torch_target_data, val_input_data, val_target_data)

# Save the trained model to ONNX format
dummy_input = torch_input_data[0].unsqueeze(0)  # Create a dummy input for the model
onnx_file_path = "regression_model.onnx"
torch.onnx.export(model, dummy_input, onnx_file_path, 
                  input_names=['input'], output_names=['output'],
                  dynamic_axes={'input': {0: 'batch_size'}, 'output': {0: 'batch_size'}})

print(f"Model has been saved to {onnx_file_path}")