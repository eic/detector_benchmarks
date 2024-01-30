import pandas as pd
import numpy as np
import seaborn as sb
import torch
import torch.nn as nn
import torch.optim as optim
import torch.optim.lr_scheduler as lr_scheduler
import matplotlib.pyplot as plt
import argparse
import sys
import hashlib

torch.set_default_dtype(torch.float32)

if torch.cuda.is_available():
  device = torch.device("cuda")
  print("GPU is available!")
else:
  device = torch.device("cpu")
  print("GPU not found. Using CPU.")

class NeuralNet(nn.Module):
    def __init__(self, size_input, size_output, n_layers, size_first_hidden_layer=128, multiplier=0.5, leak_rate=0.025):
        super().__init__()
        self.layers = nn.ModuleList()

        size_current_hidden_layer = size_first_hidden_layer

        self.layers.append(nn.Linear(size_input, size_current_hidden_layer))
        for i in range(n_layers - 2):
            self.layers.append(nn.LeakyReLU(leak_rate))
            self.layers.append(nn.Linear(size_current_hidden_layer, int(size_current_hidden_layer * multiplier)))
            size_current_hidden_layer = int(size_current_hidden_layer * multiplier)
        self.layers.append(nn.LeakyReLU(leak_rate))
        self.layers.append(nn.Linear(size_current_hidden_layer, size_output))

        print("Create a network with the following layers:")
        for layer in self.layers:
            print(layer)

    def forward(self, x):
        for layer in self.layers:
            x = layer(x)
        return x

def standardize(x):
  mean = torch.mean(x, axis=0)
  std = torch.std(x, axis=0)
  standardized_tensor = (x - mean) / std
  return standardized_tensor, mean, std

def train_model(input_tensor, target_tensor, model, hyperparameters):
  # Send model to device
  model=model.to(device)
  
  # Define the loss function and optimizer
  criterion = torch.nn.HuberLoss(reduction='mean', delta=1.0)
  optimizer = torch.optim.Adam(model.parameters(), lr=hyperparameters.learning_rate)

  # Create a learning rate scheduler
  scheduler = lr_scheduler.ReduceLROnPlateau(optimizer,'min',patience=100,cooldown=100,factor=0.5,threshold=1e-4,verbose=True)

  # Track the losses
  losses = []

  # Train the model
  for epoch in range(hyperparameters.num_epochs):
    # Forward pass
    inputs, targets = input_tensor.to(device), target_tensor.to(device)
    predictions = model(inputs)
    loss = criterion(predictions, targets)

    # Backward and optimize
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    # Track the loss value at each epoch
    losses.append(loss.item())

    # Step the learning rate scheduler
    scheduler.step(loss)

    # Print progress
    if (epoch + 1) % 10 == 0:
      print("Epoch "+str(epoch+1)+"/"+str(hyperparameters.num_epochs)+", Loss: "+"{0:0.10f}".format(loss.item()))

  # Plot the loss values
  plt.figure()
  plt.plot(range(1, hyperparameters.num_epochs+1), losses)
  plt.xlabel('Epoch')
  plt.ylabel('Loss')
  plt.title('Loss as a Function of Epoch')
  plt.yscale('log')
  plt.savefig(hyperparameters.model_dir+"/LossVsEpoch_"+hyperparameters.model_name+".png")

  torch.jit.script(model).save(hyperparameters.model_dir+"/"+hyperparameters.model_name+".pt")
  return

def run_experiment(hyperparameters):
  
  # Load training data in tensors
  training_data = pd.DataFrame()

  for i in hyperparameters.input_files:
    temp_training_data = pd.read_csv(i, delimiter='\t', header=None)
    training_data = pd.concat([training_data, temp_training_data], ignore_index=True)

  training_RP_pos_tensor = torch.tensor(training_data.iloc[:,3:7].values, dtype=torch.float32)
  training_MC_mom_tensor = torch.tensor(training_data.iloc[:,0:3].values, dtype=torch.float32)

  # Standardize training data
  match hyperparameters.model_name:
    case "model_pz":
      source = training_RP_pos_tensor
      scaled_source, mean_source, std_source = standardize(source)
      target = training_MC_mom_tensor[:,2].unsqueeze(1)
    
    case "model_py":
      source = torch.cat((training_RP_pos_tensor[:,2:4], training_MC_mom_tensor[:,2].unsqueeze(1)), 1)
      scaled_source, mean_source, std_source = standardize(source)
      target = training_MC_mom_tensor[:,1].unsqueeze(1)
   
    case "model_px":
      source = torch.cat((training_RP_pos_tensor[:,0:2], training_MC_mom_tensor[:,2].unsqueeze(1)), 1)
      scaled_source, mean_source, std_source = standardize(source)
      target = training_MC_mom_tensor[:,0].unsqueeze(1)

    case _:
      print("No model name provided. Stop further processing")
      return

  # Initialize models
  initial_model = NeuralNet(size_input=int(hyperparameters.size_input),
                               size_output=int(hyperparameters.size_output), 
                               n_layers=int(hyperparameters.n_layers),
                               size_first_hidden_layer=int(hyperparameters.size_first_hidden_layer),
                               multiplier=float(hyperparameters.multiplier),
                               leak_rate=float(hyperparameters.leak_rate)) 
 
  # Train models
  train_model(scaled_source, target, initial_model, hyperparameters)
  
  # Print end statement
  print("Training completed using "+str(len(hyperparameters.input_files))+" files with "+str(training_RP_pos_tensor.shape[0])+" eligible events")

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Train neural network model for roman pots")
  parser.add_argument('--input_files', type=str, nargs='+', required=True, help='Specify a location of input files.')  
  parser.add_argument('--model_name', type=str, required=True, help='Specify model name.')
  parser.add_argument('--model_dir', type=str, required=True, help='Specify location to save model')
  parser.add_argument('--num_epochs', type=int, required=True, help='Specify number of epochs')
  parser.add_argument('--learning_rate', type=float, required=True, help='Specify learning rate')
  parser.add_argument('--size_input', type=int, required=True, help='Specify input size')
  parser.add_argument('--size_output', type=int, required=True, help='Specify output size')
  parser.add_argument('--n_layers', type=int, required=True, help='Specify number of layers')
  parser.add_argument('--size_first_hidden_layer', type=int, required=True, help='Size of first hidden layer')
  parser.add_argument('--multiplier', type=float, required=True, help='Specify mutilplier to calculate size of subsequent hidden layers')
  parser.add_argument('--leak_rate', type=float, required=True, help='Specify leak rate')
  hyperparameters = parser.parse_args()
  run_experiment(hyperparameters)


