import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

class RegressionModel(nn.Module):
    def __init__(self):
        super(RegressionModel, self).__init__()
        self.fc1  = nn.Linear(4, 512)
        self.fc2  = nn.Linear(512, 64)
        self.fc4  = nn.Linear(64, 3)
        self.input_mean       = torch.tensor([0.0, 0.0, 0.0, 0.0])
        self.input_std        = torch.tensor([1.0, 1.0, 1.0, 1.0])
        self.input_covariance = torch.tensor([[1.0, 0.0, 0.0, 0.0],
                                              [0.0, 1.0, 0.0, 0.0],
                                              [0.0, 0.0, 1.0, 0.0],
                                              [0.0, 0.0, 0.0, 1.0]])
        self.output_mean = torch.tensor([0.0, 0.0, 0.0])
        self.output_std  = torch.tensor([1.0, 1.0, 1.0])
        self.output_correlation = torch.tensor([[1.0, 0.0, 0.0],
                                                [0.0, 1.0, 0.0],
                                                [0.0, 0.0, 1.0]])

    def forward(self, x):
        x = (x-self.input_mean)/self.input_std
        x = torch.tanh(self.fc1(x))
        x = torch.tanh(self.fc2(x))
        x = self.fc4(x)
        x = x*self.output_std + self.output_mean
        return x
    
    def adapt(self, input_data, output_data):
        in_mean = input_data.mean(axis=0)
        in_std  = input_data.std (axis=0)
        self.input_mean  = torch.tensor(in_mean)
        self.input_std   = torch.tensor(in_std)

        # Calculate the correlation matrix of the input data
        input_normalized  = (input_data-in_mean)/in_std   
        input_correlation = np.corrcoef(input_normalized, rowvar=False)         
        # Invert the correlation matrix and convert into float tensor
        self.input_covariance = torch.tensor(np.linalg.inv(input_correlation).astype(np.float32))

        self.output_mean = torch.tensor(output_data.mean(axis=0))
        self.output_std  = torch.tensor(output_data.std (axis=0))

def makeModel():
    # Create the model
    model = RegressionModel()
    # Define the optimizer
    optimizer = optim.Adam(model.parameters(), lr=0.0001)

    # Define the loss function
    criterion = nn.MSELoss()

    return model, optimizer, criterion

def trainModel(epochs, train_loader, val_loader):
    
    # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # print(f"Using device: {device}")
    
    model, optimizer, criterion = makeModel()
    # model.to(device)
    
    # Verify that the model parameters are on the GPU
    # for name, param in model.named_parameters():
    #     print(f"{name} is on {param.device}")

    # Adapt the model using the training data from the training loader
    model.adapt(train_loader.dataset.tensors[0].detach().numpy(), train_loader.dataset.tensors[1].detach().numpy())

    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        for inputs, targets in train_loader:
            # inputs, targets = inputs.to(device), targets.to(device)
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
            running_loss += loss.item() * inputs.size(0)
        
        epoch_loss = running_loss / len(train_loader.dataset)
        # print(f"Epoch [{epoch+1}/{epochs}], Loss: {epoch_loss:.4f}")

        
        # Validation step
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for val_inputs, val_targets in val_loader:
                # val_inputs, val_targets = val_inputs.to(device), val_targets.to(device)
                val_outputs = model(val_inputs)
                val_loss += criterion(val_outputs, val_targets).item() * val_inputs.size(0)
            # val_outputs = model(val_input)
            # val_loss = criterion(val_outputs, val_target)

        val_loss /= len(val_loader.dataset)

        print(f"Epoch [{epoch+1}/{epochs}], Loss: {loss}, Val Loss: {val_loss}")

    return model