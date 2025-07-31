#!/usr/bin/env python3

from ase.io import read,write
import numpy as np
import argparse,sys

def MIC(vec,cell,icell):
  # Apply minimum image convention to a vector to account for
  # periodic boundary conditions
  ivc = np.dot(icell,vec)
  rvc = np.round(ivc)
  cvc = np.dot(cell,rvc)
  ovc = vec - cvc
  return ovc

# INPUT ARGUMENTS.
parser = argparse.ArgumentParser(description="Wannier polarization")
parser.add_argument("-f", "--file", required=True, help="Input frame")
parser.add_argument("-o", "--output", required=True, help="Output file")
parser.add_argument("-n", "--nwannier", nargs="+", default=["O","4"], help="Number of Wannier centres")
parser.add_argument("-q", "--charges", nargs="+", default=["H","1","O","6","Al","3","Si","4"], help="List of charges")
parser.add_argument("-e", "--elements", default=["O"], nargs="+", help="List of elements to be assigned centres")
parser.add_argument("-tq", "--totalcharge", type=float, default=0.0, help="Expected total charge")
args = parser.parse_args()

frames = read(args.file,":")

nwannier = {}
for i in range(0,len(args.nwannier),2):
  nwannier[args.nwannier[i]] = int(args.nwannier[i+1])

charges = {}
for i in range(0,len(args.charges),2):
  charges[args.charges[i]] = float(args.charges[i+1])

for key in charges.keys():
  if not (key in nwannier):
    nwannier[key] = 0.0

for i in range(len(frames)):
  qpol = frames[i].get_cell().T * 4.8032047
  pol = np.sum([frames[i][j].position*charges[frames[i][j].symbol] for j in range(len(frames[i]))],axis=0) * 4.8032047
  for j in range(len(frames[i])):
    if (frames[i][j].symbol in args.elements):
      pol -= 4.8032047 * 2 * frames[i].info["n_wannier"] * (frames[i][j].position + frames[i].arrays["wannier_dist"][j])
  pol -= np.dot(qpol,np.round(np.dot(np.linalg.inv(qpol),pol),0))
  frames[i].info["mu"] = pol
  print(pol)
  print(frames[i].info["mu"],i,frames[i].get_cell())

# Write output
write(args.output,frames)

sys.exit(0)


for i in range(len(frames)):
  xyz = frames[i]
  cl  = xyz.get_cell().T
  qp  = 4.8032047 * cl
  icl = np.linalg.inv(cl)
  iqp = np.linalg.inv(qp)
  # Check that total charge is correct
  q = sum([charges[xyz[j].symbol] - 2*nwannier[xyz[j].symbol] for j in range(len(xyz))])
  if (q!=args.totalcharge):
    print("ERROR! q = ",q)
#    sys.exit(0)
  # Get polarizaztion
  pol = np.zeros(3,dtype=float)
  for j in range(len(xyz)):
    pol += (charges[xyz[j].symbol] - 2*nwannier[xyz[j].symbol]) * xyz[j].position
    pol -= 2*xyz.arrays["wannier_dist"][j]
  pol *= 4.8032047
  xyz.info["unwrapped_mu"] = pol
  pol = MIC(pol,qp,iqp)
  xyz.info["mu"] = pol
  frames[i] = xyz
  
