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
  def __init__(self, size_input, size_output, n_layers, size_firsthiddenlayer=128, multiplier=0.5, leak_rate=0.025):
    super().__init__()
    self.fc = []
    self.relu = []
    self.n_layers = n_layers
    self.fc.append(nn.Linear(size_input,size_firsthiddenlayer))
    self.relu.append(nn.LeakyReLU(leak_rate))
    for i in range(1,n_layers-1):
      size_currenthiddenlayer = int(size_firsthiddenlayer*multiplier**i)
      self.fc.append(nn.Linear(int(size_currenthiddenlayer/multiplier), size_currenthiddenlayer))
      self.relu.append(nn.LeakyReLU(leak_rate))
    self.fc.append(nn.Linear(size_currenthiddenlayer, size_output))
    self.fc=nn.ParameterList(self.fc)
    self.relu=nn.ParameterList(self.relu)
    print("Create a network with the linear layers "+str(self.fc))
    print("and leaky relu activation layers "+str(self.relu))

    def forward(self, x):
      for i in range(0,self.n_layers-1):
        x = self.fc[i](x)
        x = self.relu[i](x)
      x = self.fc[self.n_layers-1](x)
      return x

def standardize(tensor):
  mean = torch.mean(tensor, axis=0)
  std = torch.std(tensor, axis=0)
  standardized_tensor = (tensor - mean) / std
  return standardized_tensor

def train_model(input_tensor, target_tensor, model):

  # Define the loss function and optimizer
  criterion = torch.nn.HuberLoss(reduction='mean', delta=1.0)
  optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

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
      print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.3e}')

  # Plot the loss values
  plt.plot(range(1, num_epochs+1), losses)
  plt.xlabel('Epoch')
  plt.ylabel('Loss')
  plt.title('Loss as a Function of Epoch')
  plt.savefig('Loss vs Epoch')

  return model


