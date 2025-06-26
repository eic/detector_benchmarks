import os
from pyHepMC3 import HepMC3 as hm
import numpy as np
import argparse


PARTICLES = [
#    (111, 0.134.9766),      # pi0
    (211, 0.13957018),      # pi+
#    (-211, 0.13957018),     # pi-
#    (311, 0.497648),        # K0
    (321, 0.493677),        # K+
#    (-321, 0.493677),       # K-
    (2212, 0.938272),       # proton
#    (2112, 0.939565),       # neutron
#    (11, 0.51099895e-3),    # electron
#    (-11, 0.51099895e-3),   # positron
]


# p in GeV, angle in degree, vertex in mm
def gen_event(prange=(8, 100), arange=(0, 20)):
    evt = hm.GenEvent(hm.Units.MomentumUnit.GEV, hm.Units.LengthUnit.MM)
    pid, mass = PARTICLES[np.random.randint(len(PARTICLES))]

    # final state
    state = 1

    # momentum, angles, energy
    p = np.random.uniform(*prange)
    theta = np.random.uniform(*arange)*np.pi/180.
    phi = np.random.uniform(0., 2.*np.pi)
    e0 = np.sqrt(p*p + mass*mass)

    px = np.cos(phi)*np.sin(theta)
    py = np.sin(phi)*np.sin(theta)
    pz = np.cos(theta)

    # beam
    pbeam = hm.GenParticle(hm.FourVector(0, 0, 0, 0.938272), 2212, 4)
    ebeam = hm.GenParticle(hm.FourVector(0, 0, e0, np.sqrt(e0*e0 + 0.511e-3*0.511e-3)), -11, 4)

    hout = hm.GenParticle(hm.FourVector(px*p, py*p, pz*p, e0), pid, state)
    # evt.add_particle(part)
    vert = hm.GenVertex()
    vert.add_particle_in(ebeam)
    vert.add_particle_in(pbeam)
    vert.add_particle_out(hout)
    evt.add_vertex(vert)
    return evt


if __name__ == "__main__":
    parser = argparse.ArgumentParser('RICH dataset generator')

    parser.add_argument('output', help='path to the output file')
    parser.add_argument('-n', type=int, default=1000, dest='nev', help='number of events to generate')
    parser.add_argument('--pmin', type=float, default=8.0, dest='pmin', help='minimum momentum in GeV')
    parser.add_argument('--pmax', type=float, default=100.0, dest='pmax', help='maximum momentum in GeV')
    parser.add_argument('--angmin', type=float, default=0.0, dest='angmin', help='minimum angle in degree')
    parser.add_argument('--angmax', type=float, default=20.0, dest='angmax', help='maximum angle in degree')
    parser.add_argument('-s', type=int, default=-1, dest='seed', help='seed for random generator')

    args = parser.parse_args()

    # random seed
    if args.seed < 0:
        args.seed = os.environ.get('SEED', int.from_bytes(os.urandom(4), byteorder='big', signed=False))
    np.random.seed(args.seed)

    output = hm.WriterAscii(args.output);
    if output.failed():
        print("Cannot open file \"{}\"".format(args.output))
        sys.exit(2)

    count = 0
    while count < args.nev:
        if (count % 1000 == 0):
            print("Generated {} events".format(count), end='\r')
        evt = gen_event((args.pmin, args.pmax), (args.angmin, args.angmax))
        output.write_event(evt)
        evt.clear()
        count += 1

    print("Generated {} events".format(args.nev))
    output.close()

