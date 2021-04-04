#######################################
# EMCAL Barrel detector Reconstruction
# J.KIM 04/02/2021
#######################################

from Gaudi.Configuration import *
import os
  
from GaudiKernel.DataObjectHandleBase import DataObjectHandleBase
from Configurables import ApplicationMgr, EICDataSvc, PodioOutput, GeoSvc
from GaudiKernel import SystemOfUnits as units

detector_name = "topside"
if "JUGGLER_DETECTOR" in os.environ :
  detector_name = str(os.environ["JUGGLER_DETECTOR"])+"/"+detector_name

input_sim_file = "jug_input.root"
if "JUGGLER_SIM_FILE" in os.environ :
  input_sim_file  = str(os.environ["JUGGLER_SIM_FILE"])
else :
  print(" ERROR : JUGGLER_SIM_FILE not set" )

output_rec_file = "jug_rec.root"
if "JUGGLER_REC_FILE" in os.environ :
  output_rec_file = str(os.environ["JUGGLER_REC_FILE"])
else :
  print(" ERROR : JUGGLER_REC_FILE not set" )

n_events = 100
if "JUGGLER_N_EVENTS" in os.environ :
  n_events = str(os.environ["JUGGLER_N_EVENTS"])

geo_service  = GeoSvc("GeoSvc", detectors=["{}.xml".format(detector_name)])
podioevent   = EICDataSvc("EventDataSvc", inputs=["sim_output/{}".format(input_sim_file)], OutputLevel=DEBUG)

from Configurables import PodioInput
from Configurables import Jug__Base__InputCopier_dd4pod__Geant4ParticleCollection_dd4pod__Geant4ParticleCollection_ as MCCopier
from Configurables import Jug__Base__InputCopier_dd4pod__CalorimeterHitCollection_dd4pod__CalorimeterHitCollection_ as CalCopier

from Configurables import Jug__Digi__EcalTungstenSamplingDigi as EcalTungstenSamplingDigi

from Configurables import Jug__Reco__EcalTungstenSamplingReco as EcalTungstenSamplingReco
from Configurables import Jug__Reco__SamplingECalHitsMerger as SamplingECalHitsMerger
from Configurables import Jug__Reco__CalorimeterIslandCluster as IslandCluster
from Configurables import Jug__Reco__ClusterRecoCoG as RecoCoG

podioinput = PodioInput("PodioReader", collections=["mcparticles","EcalBarrelHits"], OutputLevel=DEBUG)

# Thrown Information
copier = MCCopier("MCCopier", 
        inputCollection="mcparticles", 
        outputCollection="mcparticles2",
        OutputLevel=DEBUG) 
# Geant4 Information
embarrelcopier = CalCopier("CalBarrelCopier", 
        inputCollection="EcalBarrelHits", 
        outputCollection="EcalBarrelHits2",
        OutputLevel=DEBUG)
# Digitization
embarreldigi = EcalTungstenSamplingDigi("ecal_barrel_digi", 
        inputHitCollection="EcalBarrelHits", 
        outputHitCollection="RawEcalBarrelHits",
        inputEnergyUnit=units.GeV,
        inputTimeUnit=units.ns,
        OutputLevel=DEBUG)
# Reconstruction
embarrelreco = EcalTungstenSamplingReco("ecal_barrel_reco", 
        inputHitCollection="RawEcalBarrelHits", 
        outputHitCollection="RecoEcalBarrelHits",
        OutputLevel=DEBUG)
# 2D+1 Clusterings
# readout id definition for barrel ecal
# <id>system:8,barrel:3,module:4,layer:6,slice:5,x:32:-16,y:-16</id>
# xy_merger sum layers/slices, masking (8+3+4, 8+3+4+5+6-1)
embarrelxymerger = SamplingECalHitsMerger("ecal_barrel_xy_merger",
        cellIDMaskRanges=[(15, 25)],
        inputHitCollection="RecoEcalBarrelHits",
        outputHitCollection="RecoEcalBarrelHitsXY")
# xy_merger sum modules, masking (8+3+4+5+6, 8+3+4+5+6+32-1)
embarrelzmerger = SamplingECalHitsMerger("ecal_barrel_z_merger",
        cellIDMaskRanges=[(26, 57)],
        inputHitCollection="RecoEcalBarrelHits",
        outputHitCollection="RecoEcalBarrelHitsZ")
# Clustering
embarrelcluster = IslandCluster("ecal_barrel_cluster",
        inputHitCollection="RecoEcalBarrelHitsXY",
        outputClusterCollection="EcalBarrelClusters",
        minClusterCenterEdep=5.0*units.MeV,
        splitCluster=False,
        groupRange=5.0)
# Reconstruct the cluster with Center of Gravity method
embarrelclusterreco = RecoCoG("ecal_barrel_clusterreco",
        clusterCollection="EcalBarrelClusters", 
        logWeightBase=6.2) 

out = PodioOutput("out", filename=output_rec_file)

out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, copier, embarrelcopier, embarreldigi, 
        embarrelreco, embarrelxymerger, embarrelzmerger, embarrelcluster, embarrelclusterreco, out],
    EvtSel = 'NONE',
    EvtMax   = n_events,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )
