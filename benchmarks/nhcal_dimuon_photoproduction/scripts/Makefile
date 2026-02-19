#===============================================================================
#   File: Makefile
#
#   Example Makefile for Pythia8 written for HF in STAR
#   Author: Thomas Ullrich
#   Last modified: September 9, 2008
#
#   This needs 3 packages installed: Pythia8, LHAPDF, and ROOT
#
#   This is setup for Mac OS X 10.5.4 but should be easy to 
#   adapt for other (Unix) platforms. In principle changing the
#   program name (PROGRAM), the location of PYTHIA (PYTHIAPATH),
#   and the LHAPDF libraries (LHAPDFPATH) should do the job.
#   Note that the environment variable ROOTSYS needs to be set.
#   Otherwise define it here in the makefile.
#===============================================================================

PROGRAM  =  pythiaDiMuon
PYTHIAPATH   = /opt/local
LHAPDFPATH = $(shell lhapdf-config --prefix)
ROOTSYS = $(shell root-config --prefix)

FASTJETCXX = `fastjet-config --cxxflags`
FASTJETLIB = `fastjet-config --libs --plugins`

DICTHEADERS = PythiaEvent.h LinkDef.h

ROOTLIB = `root-config --libdir`

CXX      =  g++
CXXFLAGS = -O -W -Wall -m64 -std=c++20 $(FASTJETCXX)
CPPFLAGS = -I$(PYTHIAPATH)/include -I$(ROOTSYS)/include/root -I$(LHAPDFPATH)/include

LDFLAGS  = -L$(PYTHIAPATH)/lib -L$(ROOTLIB) -L$(LHAPDFPATH)/lib -lpythia8 -lHepMC3 -lstdc++ -lCore -lThread -lRIO -lNet -lHist -lMathCore -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lfreetype -lpthread -lm -ldl

$(PROGRAM):	$(PROGRAM).cpp DictOutput.cxx libPythiaEvent.so
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(PROGRAM).cpp DictOutput.cxx $(LDFLAGS) $(FASTJETLIB) libPythiaEvent.so -o $(PROGRAM)
		
DictOutput.cxx : PythiaEvent.h
	rootcint -f $@ -c $(DICTHEADERS)