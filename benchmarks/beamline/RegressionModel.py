import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

class ProjectToX0Plane(nn.Module):
    def forward(self, x):
        # x shape: (batch, 6) -> [x, y, z, px, py, pz]
        x0 = x[:, 0]
        y0 = x[:, 1]
        z0 = x[:, 2]
        px = x[:, 3]
        py = x[:, 4]
        pz = x[:, 5]

        # Avoid division by zero for px
        eps = 1e-8
        px_safe = torch.where(px.abs() < eps, eps * torch.sign(px) + eps, px)
        t = -x0 / px_safe

        y_proj = y0 + py * t
        z_proj = z0 + pz * t

        # Output: [y_proj, z_proj, px, pz]
        return torch.stack([y_proj, z_proj, px, pz], dim=1)

class RegressionModel(nn.Module):
    def __init__(self):
        super(RegressionModel, self).__init__()
        self.project_to_x0 = ProjectToX0Plane()
        self.fc1  = nn.Linear(4, 512)
        self.fc2  = nn.Linear(512, 64)
        self.fc4  = nn.Linear(64, 3)

        # Normalization parameters
        self.input_mean = nn.Parameter(torch.zeros(4), requires_grad=False)
        self.input_std = nn.Parameter(torch.ones(4), requires_grad=False)
        self.output_mean = nn.Parameter(torch.zeros(3), requires_grad=False)
        self.output_std = nn.Parameter(torch.ones(3), requires_grad=False)

    def forward(self, x):
        # Apply projection and normalization
        x = self.project_to_x0(x)
        x = (x - self.input_mean) / self.input_std

        # Pass through the fully connected layers
        x = self._core_forward(x)

        # Denormalize outputs
        x = x * self.output_std + self.output_mean
        return x
    
    def _core_forward(self, x):
        # Core fully connected layers
        x = torch.tanh(self.fc1(x))
        x = torch.tanh(self.fc2(x))
        x = self.fc4(x)
        return x
    
    def adapt(self, input_data, output_data):
        # Compute normalization parameters from training data
        self.input_mean.data = torch.tensor(input_data.mean(axis=0), dtype=torch.float32)
        self.input_std.data = torch.tensor(input_data.std(axis=0), dtype=torch.float32)
        self.output_mean.data = torch.tensor(output_data.mean(axis=0), dtype=torch.float32)
        self.output_std.data = torch.tensor(output_data.std(axis=0), dtype=torch.float32)

def preprocess_data(model, data_loader):
    inputs = data_loader.dataset.tensors[0]
    targets = data_loader.dataset.tensors[1]

    # Apply projection
    projected_inputs = ProjectToX0Plane()(inputs)

    # Compute normalization parameters
    model.adapt(projected_inputs.detach().numpy(), targets.detach().numpy())

    # Normalize inputs and targets
    normalized_inputs = (projected_inputs - model.input_mean) / model.input_std
    normalized_targets = (targets - model.output_mean) / model.output_std

    # Replace the dataset with preprocessed data
    data_loader.dataset.tensors = (normalized_inputs, normalized_targets)

def makeModel():
    # Create the model
    model = RegressionModel()
    # Define the optimizer
    optimizer = optim.Adam(model.parameters(), lr=0.0001)
    # Define the loss function
    criterion = nn.MSELoss()

    return model, optimizer, criterion

def trainModel(epochs, train_loader, val_loader):
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    model, optimizer, criterion = makeModel()
    model.to(device)
    
    # Verify that the model parameters are on the GPU
    # for name, param in model.named_parameters():
    #     print(f"{name} is on {param.device}")
    
    # Preprocess training and validation data
    preprocess_data(model, train_loader)
    preprocess_data(model, val_loader)

    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        for inputs, targets in train_loader:
            # inputs, targets = inputs.to(device), targets.to(device)
            optimizer.zero_grad()
            outputs = model._core_forward(inputs)
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
                val_outputs = model._core_forward(val_inputs)
                val_loss += criterion(val_outputs, val_targets).item() * val_inputs.size(0)
            # val_outputs = model(val_input)
            # val_loss = criterion(val_outputs, val_target)

        val_loss /= len(val_loader.dataset)

        print(f"Epoch [{epoch+1}/{epochs}], Loss: {loss}, Val Loss: {val_loss}")

    return model