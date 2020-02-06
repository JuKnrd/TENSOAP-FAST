#!/usr/bin/python

import argparse
import numpy as np
from ase.io import read,write
import random
import os

parser = argparse.ArgumentParser(description="Split coordinate file")
parser.add_argument("-f", "--fname", required=True,  type=str,                help="File name")
parser.add_argument("-n", "--num",   required=True,  type=int,                help="Number of files")
args = parser.parse_args()

fname = args.fname
num   = args.num

xyz = read(args.fname,':')

split_xyz = np.array_split(xyz,num)

print "Created %i sets"%num
for i in xrange(num):
    fld = 'run_' + str(i+1)
    if not os.path.exists(fld)):
        os.mkdir(fld)
    write(fld + '/' + fname,split_xyz[i])
