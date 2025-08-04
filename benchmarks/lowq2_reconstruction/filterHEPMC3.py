
from pyHepMC3 import HepMC3 as hm
from pyHepMC3.rootIO import HepMC3 as hmrootIO
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description='Train a regression model for the Tagger.')
parser.add_argument('--outFile', type=str, default="temp.hepmc3.tree.root", help='Path to the output file')
parser.add_argument('--inFile', type=str, nargs='+', help='Path to the input files')

args = parser.parse_args()

input_file = args.inFile[0]   # Change to your input file
output_file = args.outFile

# Initialize reader and writer
reader = hm.deduce_reader(input_file)
if not reader:
    print(f"Error: Could not open input file {input_file}", file=sys.stderr)
    sys.exit(1)
writer = hmrootIO.WriterRootTree(output_file)
if not writer:
    print(f"Error: Could not create output file {output_file}", file=sys.stderr)
    sys.exit(1)

event = hm.GenEvent()
while not reader.failed():
    reader.read_event(event)
    for p in list(event.particles()):
        if p.pid() != 11:
            event.remove_particle(p)
    # Only write events that still have electrons
    if any(p.pid() == 11 for p in event.particles()):
        writer.write_event(event)
    event.clear()

reader.close()
writer.close()
print(f"Filtered file written to {output_file}")