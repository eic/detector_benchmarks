import torch
import torch.nn as nn
import torch.optim as optim

class RegressionModel(nn.Module):
    def __init__(self):
        super(RegressionModel, self).__init__()
        self.fc1 = nn.Linear(4, 1024)
        self.fc2 = nn.Linear(1024, 128)
        self.fc3 = nn.Linear(128, 64)
        self.fc4 = nn.Linear(64, 32)
        self.fc5 = nn.Linear(32, 3)
        self.mean = torch.tensor([0.0, 0.0, 0.0, 0.0])
        self.std  = torch.tensor([1.0, 1.0, 1.0, 1.0])

    def forward(self, x):
        x = (x-self.mean)/self.std
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = torch.relu(self.fc3(x))
        x = torch.relu(self.fc4(x))
        x = self.fc5(x)
        return x
    
    def adapt(self, input_data):
        self.mean = input_data.mean(axis=0)
        self.std = input_data.std(axis=0)
    

def makeModel():
    # Create the model
    model = RegressionModel()
    # Define the optimizer
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # Define the loss function
    criterion = nn.MSELoss()

    return model, optimizer, criterion

def trainModel(epochs, input_data, target_data, val_input, val_target):
    model, optimizer, criterion = makeModel()
    model.adapt(input_data)
    for epoch in range(epochs):
        model.train()
        # Zero the parameter gradients
        optimizer.zero_grad()

        # Forward pass
        output = model(input_data)
        loss = criterion(output, target_data)

        # Backward pass
        loss.backward()
        optimizer.step()
        
        # Validation step
        model.eval()
        with torch.no_grad():
            val_outputs = model(val_input)
            val_loss = criterion(val_outputs, val_target)

        print(f"Epoch [{epoch+1}/{epochs}], Loss: {loss.item()}, Val Loss: {val_loss.item()}")

    return model