import torch
import argparse
from ProcessData import create_arrays
from torch.utils.data import DataLoader, TensorDataset
from RegressionModel  import trainModel

# Parse arguments
parser = argparse.ArgumentParser(description='Train a regression model for the Tagger.')
parser.add_argument('--dataFiles', type=str, nargs='+', help='Path to the data files')
parser.add_argument('--outModelFile', type=str, default="regression_model.onnx", help='Output file for the trained model')
parser.add_argument('--batchSize', type=int, default=1024, help='Batch size for training')
parser.add_argument('--epochs', type=int, default=1000, help='Number of epochs for training')
parser.add_argument('--entries', type=int, default=None, help='Number of entries to process from the data files')

args   = parser.parse_args()

input_data, target_data = create_arrays(args.dataFiles, entries=args.entries)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")
if device.type == 'cuda':
    print("Device:", torch.cuda.get_device_name(0))
    
torch_input_data  = torch.tensor(input_data, dtype=torch.float32)
torch_target_data = torch.tensor(target_data, dtype=torch.float32)

print(f"Input data shape: {torch_input_data.shape}")
print(f"Target data shape: {torch_target_data.shape}")

# Split data into training and validation sets
validation_fraction = 0.25
split_index         = int(len(torch_input_data) * (1 - validation_fraction))

val_input_data    = torch_input_data[split_index:]
val_target_data   = torch_target_data[split_index:]
train_input_data  = torch_input_data[:split_index]
train_target_data = torch_target_data[:split_index] 

# Create TensorDatasets
train_dataset = TensorDataset(train_input_data, train_target_data)
val_dataset   = TensorDataset(val_input_data, val_target_data)

# Create DataLoaders
train_loader = DataLoader(train_dataset, batch_size=args.batchSize, shuffle=True )
val_loader   = DataLoader(val_dataset,   batch_size=args.batchSize, shuffle=False)

# Train the model
print(f"Training data: {len(train_input_data)} samples")
model  = trainModel(args.epochs, train_loader, val_loader, device)

# Save the trained model to ONNX format
dummy_input = torch_input_data[0].unsqueeze(0).to(device)  # Create a dummy input for the model
model.to(device)  # Ensure the model is on the correct device

torch.onnx.export(model, dummy_input, args.outModelFile, 
                  input_names=['input'], output_names=['output'],
                  dynamic_axes={'input': {0: 'batch_size'}, 'output': {0: 'batch_size'}})

print(f"Model has been saved to {args.outModelFile}")