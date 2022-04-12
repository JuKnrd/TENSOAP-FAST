#!/usr/bin/env python3

import argparse
import sys,os
import numpy as np
from ase.data import atomic_numbers

parser = argparse.ArgumentParser(description="Extra LODE parameters for TENSOAP-FAST model")
parser.add_argument("-nn",   "--nonorm",                            action='store_true', help="Do not normalize power spectrum")
parser.add_argument("-sew",  "--sigewald", type=float, default=1.0,                      help="Gaussian width for Ewald splitting")
parser.add_argument("-rad",  "--radsize",  type=int,   default=50,                       help="Dimension of the Gauss-Legendre grid needed for numerical radial integration of the direct-space potential")
parser.add_argument("-leb",  "--lebsize",  type=int,   default=146,                      help="Dimension of Lebedev grid used for numerical angular integration of the direct-space potential")
parser.add_argument("-o",    "--outfile",  type=str,   required=True,                    help="Output file")
args = parser.parse_args()

# Check whether Lebedev grid size is an allowed value
if not (args.lebsize in [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]):
  print("ERROR: Lebedev grid size not allowed!")
  sys.exit(0)

model = np.zeros(4,float)

model[0] = args.nonorm
model[1] = args.sigewald
model[2] = args.radsize
model[3] = args.lebsize

# Save file
output = args.outfile
if (output.count('/')>0):
  op = output.split('/')
  folder = '/'.join(op[:len(op)-1])
  if (not os.path.isdir(folder)):
    split_folder = op[:len(op)-1]
    make_folder = ''
    for i in range(len(split_folder)):
      make_folder += split_folder[i] + '/'
      if (not os.path.isdir(make_folder)):
        os.mkdir(make_folder)
mfl = output + '.mdl'
model.tofile(mfl)
