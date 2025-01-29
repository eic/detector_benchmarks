import uproot
import awkward as ak

def create_arrays(dataFiles):

    # List of branches to load
    branches = ["_TaggerTrackerFeatureTensor_shape","_TaggerTrackerFeatureTensor_floatData","_TaggerTrackerTargetTensor_floatData"]

    # Load data from concatenated list of files
    data = uproot.concatenate([f"{file}:events" for file in dataFiles], branches, library="ak")
    
    # Filter events with at least one track
    num_tracks = data["_TaggerTrackerFeatureTensor_shape"][:,0]
    filtered_data = data[num_tracks == 1]

    input_data = filtered_data["_TaggerTrackerFeatureTensor_floatData"]
    target_data = filtered_data["_TaggerTrackerTargetTensor_floatData"]

    return input_data, target_data

