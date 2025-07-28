import uproot
import awkward as ak

def create_arrays(dataFiles,featureName="_TaggerTrackerFeatureTensor_floatData",targetName="_TaggerTrackerTargetTensor_floatData", entries=None, treeName="events"):

    # List of branches to load
    branches = [featureName,targetName]

    # Load data from concatenated list of files
    data = uproot.concatenate([f"{file}:{treeName}" for file in dataFiles], branches, entry_stop=entries, library="ak")

    input_data  = data[featureName]
    target_data = data[targetName]

    return input_data, target_data