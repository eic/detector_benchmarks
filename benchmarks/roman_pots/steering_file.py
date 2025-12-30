#!/usr/bin/env python
"""
DD4hep simulation with some argument parsing
Based on M. Frank and F. Gaede runSim.py
   @author  A.Sailer
   @version 0.1
Modified with standard EIC EPIC requirements.
"""
from __future__ import absolute_import, unicode_literals
import logging
import sys
import math

#from DDSim import SystemOfUnits
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV, radian

SIM = DD4hepSimulation()
################################################################
#  Steering file to shoot particles into the Roman Pots
#  used to generate the training data for RP ML algorithm.
#
#  The default magnetic field settings are for 18x275 GeV.
#################################################################

#############################################################
#   Particle Gun simulations use these options
#   --> comment out if using MC input
#############################################################

#SIM.enableDetailedShowerMode = False
SIM.enableG4GPS = False
SIM.enableG4Gun = False
SIM.enableGun = True
#SIM.gun.direction = (math.sin(-0.025*radian), 0, math.cos(-0.025*radian))
SIM.crossingAngleBoost = -0.025*radian


## InputFiles for simulation .stdhep, .slcio, .HEPEvt, .hepevt, .hepmc, .pairs files are supported
SIM.inputFiles = []
## Macro file to execute for runType 'run' or 'vis'
SIM.macroFile = ""
## number of events to simulate, used in batch mode
##SIM.numberOfEvents = 5000
## Outputfile from the simulation,only lcio output is supported
##SIM.outputFile = "RP_eDep_and_thresholds_10_6_2023_1k_events.edm4hep.root"
SIM.runType = "run"
## FourVector of translation for the Smearing of the Vertex position: x y z t
SIM.vertexOffset = [0.0, 0.0, 0.0, 0.0]
## FourVector of the Sigma for the Smearing of the Vertex position: x y z t
SIM.vertexSigma = [0.0, 0.0, 0.0, 0.0]



################################################################################
## Configuration for the magnetic field (stepper)
################################################################################
SIM.field.delta_chord = 1e-03
SIM.field.delta_intersection = 1e-03
SIM.field.delta_one_step = .5e-01*mm
SIM.field.eps_max = 1e-03
SIM.field.eps_min = 1e-04
#SIM.field.equation = "Mag_UsualEqRhs"
SIM.field.largest_step = 100.0*mm
SIM.field.min_chord_step = 1.e-2*mm
SIM.field.stepper = "HelixSimpleRunge"

## choose the distribution of the random direction for theta
##
##     Options for random distributions:
##
##     'uniform' is the default distribution, flat in theta
##     'cos(theta)' is flat in cos(theta)
##     'eta', or 'pseudorapidity' is flat in pseudorapity
##     'ffbar' is distributed according to 1+cos^2(theta)
##
##     Setting a distribution will set isotrop = True
##
SIM.gun.distribution = 'uniform'
SIM.gun.energy = 275.0*GeV ## default energy value is MeV
#SIM.gun.momentumMin = 274.0*GeV
#SIM.gun.momentumMax = 275.0*GeV

##  isotropic distribution for the particle gun
##
##     use the options phiMin, phiMax, thetaMin, and thetaMax to limit the range of randomly distributed directions
##     if one of these options is not None the random distribution will be set to True and cannot be turned off!
##
SIM.gun.isotrop = True
#SIM.gun.multiplicity = 1
SIM.gun.particle = "proton"

## Minimal azimuthal angle for random distribution
##SIM.gun.phiMin = 0

##  position of the particle gun, 3 vector. unit mm
##SIM.gun.position = (0, 0, 0)
SIM.gun.thetaMin = 0.000*radian
SIM.gun.thetaMax = 0.005*radian


################################################################################
## Configuration for the output levels of DDG4 components
################################################################################

## Output level for input sources
SIM.output.inputStage = 3

## Output level for Geant4 kernel
SIM.output.kernel = 3

## Output level for ParticleHandler
SIM.output.part = 3

## Output level for Random Number Generator setup
SIM.output.random = 6


################################################################################
## Configuration for the Particle Handler/ MCTruth treatment
################################################################################

## Enable lots of printout on simulated hits and MC-truth information
SIM.part.enableDetailedHitsAndParticleInfo = False

##  Keep all created particles
SIM.part.keepAllParticles = True

## Minimal distance between particle vertex and endpoint of parent after
##     which the vertexIsNotEndpointOfParent flag is set
##
SIM.part.minDistToParentVertex = 2.2e-14

## MinimalKineticEnergy to store particles created in the tracking region
#SIM.part.minimalKineticEnergy = 1*MeV

##  Printout at End of Tracking
SIM.part.printEndTracking = False

##  Printout at Start of Tracking
SIM.part.printStartTracking = False

## List of processes to save, on command line give as whitespace separated string in quotation marks
SIM.part.saveProcesses = ['Decay']


################################################################################
## Configuration for the PhysicsList
################################################################################
SIM.physics.decays = True
SIM.physics.list = "FTFP_BERT" # "FTFP_BERT"

################################################################################
## Properties for the random number generator
################################################################################

## If True, calculate random seed for each event based on eventID and runID
## allows reproducibility even when SkippingEvents
SIM.random.enableEventSeed = False
SIM.random.file = None
SIM.random.luxury = 1
SIM.random.replace_gRandom = True
SIM.random.seed = None
SIM.random.type = None
