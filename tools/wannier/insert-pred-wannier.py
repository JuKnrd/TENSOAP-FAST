#!/usr/bin/env python3

from ase.io import read,write
import numpy as np
import argparse,sys

# INPUT ARGUMENTS.
parser = argparse.ArgumentParser(description="Wannier polarization")
parser.add_argument("-f", "--file", required=True, help="Input frame")
parser.add_argument("-o", "--output", required=True, help="Output file")
parser.add_argument("-w", "--wannier", required=True, help="Wannier centres")
parser.add_argument("-e", "--elements", nargs="+", default=["O"], help="Elements to be assigned Wannier centres")
args = parser.parse_args()

frames = read(args.file,":")

wannier = np.loadtxt(args.wannier)

# Remove Wannier distances if they are present
for i in range(len(frames)):
  if 'wannier_dist' in frames[i].arrays:
    del frames[i].arrays['wannier_dist']

# Put in new Wannier distances
wdist = np.array([np.zeros((len(frames[i]),3)) for i in range(len(frames))])

k = -1
for i in range(len(frames)):
  for j in range(len(frames[i])):
    if frames[i][j].symbol in args.elements:
      k += 1
      wdist[i,j] = wannier[k]
  frames[i].arrays['wannier_dist'] = wdist[i]

write(args.output,frames)
