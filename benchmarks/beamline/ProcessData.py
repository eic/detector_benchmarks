import uproot
import awkward as ak

def create_arrays(dataFiles,entries=None):

    # List of branches to load
    branches = ["_TaggerTrackerFeatureTensor_floatData","_TaggerTrackerTargetTensor_floatData"]

    # Load data from concatenated list of files
    data = uproot.concatenate([f"{file}:events" for file in dataFiles], branches, entry_stop=entries, library="ak")

    input_data  = data["_TaggerTrackerFeatureTensor_floatData"]
    target_data = data["_TaggerTrackerTargetTensor_floatData"]

    return input_data, target_data