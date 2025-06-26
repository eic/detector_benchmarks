import uproot
import awkward as ak

def create_arrays(dataFiles,beamEnergy=18):

    # List of branches to load
    branches = ["features","targets"]

    # Load data from concatenated list of files
    data = uproot.concatenate([f"{file}:events" for file in dataFiles], branches, library="ak")
    
    input_data = data["features"]
    target_data = data["targets"]/beamEnergy

    return input_data, target_data