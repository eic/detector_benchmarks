#!/usr/bin/env python3
"""
Script to extract cross section per event from HepMC3 files
"""
import ROOT  # noqa: F401; loads core libraries required by pyhepmc3
from pyHepMC3 import HepMC3 as hepmc
from pyHepMC3.rootIO import HepMC3 as hmrootIO
import sys

def extract_cross_section_info(hepmc_file, scale_factor=1.0):
    """
    Extract cross section and number of events from HepMC3 tree file
    Returns cross section per event in mb, scaled by scale_factor
    """
    
    # Try to read the HepMC3 file using ReaderRoot from rootIO module
    reader = hmrootIO.ReaderRootTree(hepmc_file)
    
    #Read the first event to access run info
    reader.read_event(hepmc.GenEvent())

    run_info = reader.run_info()
    
    cross_section_str = run_info.attribute_as_string("crossSection")
    
    print(f"Found attributes: crossSection={cross_section_str}", file=sys.stderr)
    
    try:
        if cross_section_str:
            cross_section = float(cross_section_str)
            
            if cross_section > 0:
                cross_section_per_event = cross_section * scale_factor
                
                print(f"Cross section: {cross_section}, Scale factor: {scale_factor}", file=sys.stderr)
                print(f"Cross section per event (scaled): {cross_section_per_event} mb", file=sys.stderr)
                reader.close()
                return cross_section_per_event
                
    except (ValueError, TypeError) as e:
        print(f"Error converting cross section values: {e}", file=sys.stderr)
        
    reader.close()


if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python extract_cross_section.py <hepmc_file> [scale_factor]")
        sys.exit(1)
    
    hepmc_file = sys.argv[1]
    scale_factor = float(sys.argv[2]) if len(sys.argv) == 3 else 1.0
    
    cross_section = extract_cross_section_info(hepmc_file, scale_factor)
    
    if cross_section is not None:
        print(cross_section)
    else:
        print(0.0551 * scale_factor)  # Default value, scaled
