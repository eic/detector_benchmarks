import pandas as pd
import numpy as np
import seaborn as sb
import torch
import torch.nn as nn
import torch.optim as optim
import torch.optim.lr_scheduler as lr_scheduler
import matplotlib.pyplot as plt
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


