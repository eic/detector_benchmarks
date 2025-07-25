import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

class ProjectToX0Plane(nn.Module):
    def forward(self, x):
        # x shape: (batch, 6) -> [x, y, z, px, py, pz]
        x0, y0, z0, px, py, pz = x.unbind(dim=1)

        # Normalize momentum components
        momentum = torch.sqrt(px**2 + py**2 + pz**2)
        px_norm = px / momentum
        py_norm = py / momentum
        pz_norm = pz / momentum

        # Avoid division by zero for px
        # eps = 1e-8
        # px_safe = torch.where(px_norm.abs() < eps, eps * torch.sign(px_norm) + eps, px_norm)
        t = -x0 / px_norm

        y_proj = y0 + py_norm * t
        z_proj = z0 + pz_norm * t

        # Output: [y_proj, z_proj, px_norm, py_norm]
        return torch.stack([y_proj, z_proj, px_norm, py_norm], dim=1)
    
    def project_numpy(self, arr):
        """
        Projects a numpy array of shape (N, 6) using the forward method,
        returns a numpy array of shape (N, 4).
        """
        device = next(self.parameters()).device if any(p.device.type != 'cpu' for p in self.parameters()) else 'cpu'
        x = torch.from_numpy(arr).float().to(device)
        with torch.no_grad():
            projected = self.forward(x)
        return projected.cpu().numpy()

class RegressionModel(nn.Module):
    def __init__(self):
        super(RegressionModel, self).__init__()
        self.project_to_x0 = ProjectToX0Plane()
        self.fc1  = nn.Linear(4, 512)
        self.fc2  = nn.Linear(512, 64)
        self.fc3  = nn.Linear(64, 3)  # Output layer for

        # Normalization parameters
        self.input_mean = nn.Parameter(torch.zeros(4), requires_grad=False)
        self.input_std = nn.Parameter(torch.ones(4), requires_grad=False)
        self.output_mean = nn.Parameter(torch.zeros(3), requires_grad=False)
        self.output_std = nn.Parameter(torch.ones(3), requires_grad=False)
        

    def forward(self, x):
        # Apply projection
        x = self.project_to_x0(x)        
        # Normalize inputs
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
        x = self.fc3(x)
        return x
    
    def adapt(self, input_data, output_data):
        # Normalization
        self.input_mean.data = input_data.mean(dim=0)
        self.input_std.data = input_data.std(dim=0)
        self.output_mean.data = output_data.mean(dim=0)
        self.output_std.data = output_data.std(dim=0)

def preprocess_data(model, data_loader, adapt=True):
    inputs  = data_loader.dataset.tensors[0]
    targets = data_loader.dataset.tensors[1]

    # Apply projection
    projected_inputs = ProjectToX0Plane()(inputs)

    # Compute normalization parameters
    if adapt:
        model.adapt(projected_inputs, targets)

    # Normalize inputs and targets
    normalized_inputs  = (projected_inputs - model.input_mean ) / model.input_std
    normalized_targets = (targets          - model.output_mean) / model.output_std

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

def trainModel(epochs, train_loader, val_loader, device):
    
    model, optimizer, criterion = makeModel()
    
    model.to(device)

    # Preprocess training and validation data
    preprocess_data(model, train_loader, adapt=True)

    # Preprocess validation data without adapting
    preprocess_data(model, val_loader, adapt=False)

    # Move data to the GPU
    train_loader.dataset.tensors = (train_loader.dataset.tensors[0].to(device), train_loader.dataset.tensors[1].to(device))
    val_loader.dataset.tensors = (val_loader.dataset.tensors[0].to(device), val_loader.dataset.tensors[1].to(device))

    # Verify that the model parameters are on the GPU
    for name, param in model.named_parameters():
        print(f"{name} is on {param.device}")

    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        for inputs, targets in train_loader:
            optimizer.zero_grad()
            outputs = model._core_forward(inputs)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
            running_loss += loss.item() * inputs.size(0)
        
        epoch_loss = running_loss / len(train_loader.dataset)

        
        # Validation step
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for val_inputs, val_targets in val_loader:
                val_outputs = model._core_forward(val_inputs)
                val_loss += criterion(val_outputs, val_targets).item() * val_inputs.size(0)

        val_loss /= len(val_loader.dataset)

        print(f"Epoch [{epoch+1}/{epochs}], Loss: {epoch_loss}, Val Loss: {val_loss}")

    return model