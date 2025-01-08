import uproot
import awkward as ak
import numpy as np
def create_arrays(dataFiles,beamEnergy):

    # List of branches to load
    branches = ["TaggerTrackerProjectedTracks","MCParticles","MCScatteredElectrons_objIdx"]

    # Load data from concatenated list of files
    # file = uproot.open(dataFiles)
    # tree = file['events']
    # data = tree.arrays(branches, library="ak")
    data = uproot.concatenate([f"{file}:events" for file in dataFiles], branches, library="ak")

    # Filter events with at least one track
    tracks = data["TaggerTrackerProjectedTracks"]
    num_tracks = ak.num(tracks["TaggerTrackerProjectedTracks.pdg"])
    filtered_data = data[num_tracks == 1]

    # Filter tracks for intput data
    filtered_tracks    = filtered_data["TaggerTrackerProjectedTracks"]
    ia = filtered_tracks["TaggerTrackerProjectedTracks.loc.a"]
    ib = filtered_tracks["TaggerTrackerProjectedTracks.loc.b"]
    ipx = np.sin(filtered_tracks["TaggerTrackerProjectedTracks.theta"]) * np.cos(filtered_tracks["TaggerTrackerProjectedTracks.phi"])
    ipy = np.sin(filtered_tracks["TaggerTrackerProjectedTracks.theta"]) * np.sin(filtered_tracks["TaggerTrackerProjectedTracks.phi"])

    # Filter particle array to select scattered electron for target data
    electron_idx = filtered_data["MCScatteredElectrons_objIdx"]["MCScatteredElectrons_objIdx.index"][:,0]
    filtered_particles = filtered_data["MCParticles"][np.arange(len(electron_idx)), electron_idx]

    # Normalize the target variables
    tpx = filtered_particles["MCParticles.momentum.x"]/beamEnergy
    tpy = filtered_particles["MCParticles.momentum.y"]/beamEnergy
    tpz = filtered_particles["MCParticles.momentum.z"]/beamEnergy

    input_data = ak.concatenate([ia,ib,ipx,ipy], axis=1)
    target_data = ak.concatenate([tpx[:,None], tpy[:,None], tpz[:,None]], axis=1)

    return input_data, target_data

