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

def train_model(input_tensor, target_tensor, model, num_epochs, learning_rate):
  # Send model to device
  model=model.to(device)
  
  # Define the loss function and optimizer
  criterion = torch.nn.HuberLoss(reduction='mean', delta=1.0)
  optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

  # Create a learning rate scheduler
  scheduler = lr_scheduler.ReduceLROnPlateau(optimizer,'min',patience=100,cooldown=100,factor=0.5,threshold=1e-4,verbose=True)

  # Track the losses
  losses = []

  # Train the model
  for epoch in range(num_epochs):
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
      print("Epoch "+str(epoch+1)+"/"+str(num_epochs)+", Loss: "+"{0:0.10f}".format(loss.item()))

  # Plot the loss values
  plt.plot(range(1, num_epochs+1), losses)
  plt.xlabel('Epoch')
  plt.ylabel('Loss')
  plt.title('Loss as a Function of Epoch')
  plt.savefig('Loss vs Epoch')

  return model

def run_experiment(hyperparameters):
  
  # Load input and target training data in tensors
  training_RP_pos = pd.DataFrame()
  training_MC_mom = pd.DataFrame()

  for i in range(1,int(hyperparameters.num_training_inputs)+1):
    temp_training_RP_pos = pd.read_csv(hyperparameters.input_files+str(i)+'.txt', delimiter='\t', header=None)
    training_RP_pos = pd.concat([training_RP_pos, temp_training_RP_pos], ignore_index=True)
    temp_training_MC_mom = pd.read_csv(hyperparameters.target_files+str(i)+'.txt', delimiter='\t', header=None)
    training_MC_mom = pd.concat([training_MC_mom, temp_training_MC_mom], ignore_index=True)

  training_RP_pos_tensor = torch.tensor(training_RP_pos.values, dtype=torch.float32)
  training_MC_mom_tensor = torch.tensor(training_MC_mom.values, dtype=torch.float32)

  # Standardize training data
  source_pz = training_RP_pos_tensor
  scaled_source_pz, mean_source_pz, std_source_pz = standardize(source_pz)
  target_pz = training_MC_mom_tensor[:,2].unsqueeze(1)

  source_py = torch.cat((training_RP_pos_tensor[:,2:4], training_MC_mom_tensor[:,2].unsqueeze(1)), 1)
  scaled_source_py, mean_source_py, std_source_py = standardize(source_py)
  target_py = training_MC_mom_tensor[:,1].unsqueeze(1)

  source_px = torch.cat((training_RP_pos_tensor[:,0:2], training_MC_mom_tensor[:,2].unsqueeze(1)), 1)
  scaled_source_px, mean_source_px, std_source_px = standardize(source_px)
  target_px = training_MC_mom_tensor[:,0].unsqueeze(1)

  # Initialize models
  initial_model_pz = NeuralNet(size_input=hyperparameters.size_input_pz,
                               size_output=hyperparameters.size_output_pz, 
                               n_layers=hyperparameters.n_layers_pz,
                               size_first_hidden_layer=hyperparameters.size_first_hidden_layer_pz,
                               multiplier=hyperparameters.multiplier_pz,
                               leak_rate=hyperparameters.leak_rate_pz)
  initial_model_py = NeuralNet(size_input=hyperparameters.size_input_py,
                               size_output=hyperparameters.size_output_py, 
                               n_layers=hyperparameters.n_layers_py,
                               size_first_hidden_layer=hyperparameters.size_first_hidden_layer_py,
                               multiplier=hyperparameters.multiplier_py,
                               leak_rate=hyperparameters.leak_rate_py)
  initial_model_px = NeuralNet(size_input=hyperparameters.size_input_px,
                               size_output=hyperparameters.size_output_px, 
                               n_layers=hyperparameters.n_layers_px,
                               size_first_hidden_layer=hyperparameters.size_first_hidden_layer_px,
                               multiplier=hyperparameters.multiplier_px,
                               leak_rate=hyperparameters.leak_rate_px)
  
  # Train models
  model_pz = train_model(scaled_source_pz, target_pz, initial_model_pz, num_epochs=hyperparameters.num_epochs_pz, learning_rate=hyperparameters.learning_rate_pz)
  model_py = train_model(scaled_source_py, target_py, initial_model_py, num_epochs=hyperparameters.num_epochs_py, learning_rate=hyperparameters.learning_rate_py)
  model_px = train_model(scaled_source_px, target_px, initial_model_px, num_epochs=hyperparameters.num_epochs_px, learning_rate=hyperparameters.learning_rate_px)

  # Save models
  torch.jit.script(model_pz).save('model_pz.pt')
  torch.jit.script(model_py).save('model_py.pt')
  torch.jit.script(model_px).save('model_px.pt')


if __name__ == "__main__":
  parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
  hyperparameters_list = ['--input_files', '--target_files', '--num_training_inputs', 
                   '--num_epochs_pz', '--learning_rate_pz', '--size_input_pz', '--size_output_pz', '--n_layers_pz', '--size_first_hidden_layer_pz', '--multiplier_pz', '--leak_rate_pz',
                   '--num_epochs_py', '--learning_rate_py', '--size_input_py', '--size_output_py', '--n_layers_py', '--size_first_hidden_layer_py', '--multiplier_py', '--leak_rate_py',
                   '--num_epochs_px', '--learning_rate_px', '--size_input_px', '--size_output_px', '--n_layers_px', '--size_first_hidden_layer_px', '--multiplier_px', '--leak_rate_px']
  for hyperparameter in hyperparameters_list:
    parser.add_argument(hyperparameter)
  hyperparameters = parser.parse_args(['@'+str(sys.argv[1])])
  run_experiment(hyperparameters)


