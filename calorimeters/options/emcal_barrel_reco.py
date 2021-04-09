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

podioinput = PodioInput("PodioReader", collections=["mcparticles","EcalBarrelAstroPixHits"], OutputLevel=DEBUG)

# Thrown Information
copier = MCCopier("MCCopier", 
        inputCollection="mcparticles", 
        outputCollection="mcparticles2",
        OutputLevel=DEBUG) 
# Geant4 Information
embarrelcopier = CalCopier("CalBarrelCopier", 
        inputCollection="EcalBarrelAstroPixHits", 
        outputCollection="EcalBarrelAstroPixHits2",
        OutputLevel=DEBUG)

out = PodioOutput("out", filename=output_rec_file)

out.outputCommands = ["keep *"]

ApplicationMgr(
    TopAlg = [podioinput, copier, embarrelcopier, out],
    EvtSel = 'NONE',
    EvtMax   = n_events,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
 )
