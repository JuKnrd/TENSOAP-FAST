#!/usr/bin/python

import argparse
import sys,os
import numpy as np
from ase.data import atomic_numbers

parser = argparse.ArgumentParser(description="Put together SOAPFAST model for prediction program")
parser.add_argument("-ps",   "--power",     type=str,   required=True,                     help="Training power spectrum file")
parser.add_argument("-sf",  "--sparse",     type=str,   required=True,                     help="Sparsification file")
parser.add_argument("-w",   "--weights",    type=str,   required=True,                     help="Weights file")
parser.add_argument("-p0",  "--power0",     type=str,   required=False,                    help="Scalar training power spectrum file")
parser.add_argument("-sf0", "--sparse0",    type=str,   required=False,                    help="Scalar sparsification file")
parser.add_argument("-hp",  "--hyperparam", type=str,   required=False,        nargs='+',  help="Hyperparameter string")
parser.add_argument("-o",   "--outfile",    type=str,   required=True,                     help="Output file")
parser.add_argument("-n",   "--nmax",       type=int,   default=-1,            nargs='+',  help="nmax")
parser.add_argument("-l",   "--lmax",       type=int,   default=-1,            nargs='+',  help="lmax")
parser.add_argument("-rc",  "--rcut",       type=float, default=-1.0,          nargs='+',  help="rcut")
parser.add_argument("-pr",  "--periodic",   action='store_true',                           help="Periodic system")
parser.add_argument("-lm",  "--lambdaval",  type=int,   default=0,                         help="Spherical order")
parser.add_argument("-z",   "--zeta",       type=int,   default=1,                         help="zeta")
parser.add_argument("-sg",  "--sigma",      type=float, default=-1.0,          nargs='+',  help="sigma")
parser.add_argument("-rs",  "--radial",     type=float, default=[0.0,0.0,0.0], nargs='+',  help="radial scaling")
parser.add_argument("-c",   "--centres",    type=str,                          nargs='+',  help="centres")
parser.add_argument("-s",   "--species",    type=str,                          nargs='+',  help="species")
args = parser.parse_args()

do_scalar = False

if (args.power0!=None or args.sparse0!=None):
    # Check that all of them are specified, otherwise exit
    if (not (args.power0!='' and args.sparse0!='')):
        print "ERROR: all scalar files must be specified!"
        sys.exit(0)
    do_scalar = True

if (len(args.radial)!=3 and len(args.radial)!=6):
    print "ERROR: three or six parameters are required for radial scaling!"
    sys.exit(0)

if (args.centres==None):
    numcen = 0
else:
    numcen = len(args.centres)
if (args.species==None):
    numspec = 0
else:
    numspec = len(args.species)

if (not do_scalar):
    pp = np.load(args.power)
    fp = np.load(args.sparse + '_fps.npy')
    am = np.load(args.sparse + '_Amat.npy')
    weights = np.load(args.weights)
    wt = weights[4]
    if (len(weights)>5):
        mv = np.load(args.weights)[5]
    else:
        mv = 0.0
    # Make model array
    nmol = np.shape(pp)[0]
    ncut = np.shape(pp)[-1]
    if (len(np.shape(pp))==3):
        degen = 1
    else:
        degen = np.shape(pp)[-2]
    ln = 5 + np.size(pp) + 1 + np.size(fp) + 2*np.size(am) + 1 + np.size(wt) + 7 + 1 + numcen + 1 + numspec
    model = np.zeros(ln,float)
    p1 = np.zeros((nmol,1,degen,ncut),float)
    if (degen==1):
        for j in xrange(nmol):
            for k in xrange(ncut):
                p1[j,0,0,k] = pp[j,0,k]
    else:
        for j in xrange(nmol):
            for k in xrange(degen):
                for l in xrange(ncut):
                    p1[j,0,k,l] = pp[j,0,k,l]
    model[0] = 0.0
    model[1] = args.lambdaval
    model[2] = args.zeta
    if (args.periodic):
        model[3] = 1.0
    else:
        model[3] = 0.0
    model[3] = nmol
    model[4] = ncut
    i = 5
    for j in xrange(nmol):
        for k in xrange(degen):
            for l in xrange(ncut):
                i+=1
                model[i] = p1[j,0,k,l]
    i+=1
    model[i] = ncut
    for j in xrange(ncut):
        i+=1
        model[i] = fp[j]
    for j in xrange(ncut):
        for k in xrange(ncut):
            i+=1
            model[i] = np.real(am[j,k])
            i+=1
            model[i] = np.imag(am[j,k])
    i+=1
    model[i] = mv
    for j in xrange(nmol*degen):
        i+=1
        model[i] = wt[j]
    i+=1
    print len(args.nmax)
    model[i] = float(args.nmax)
    i+=1
    model[i] = float(args.lmax)
    i+=1
    model[i] = args.rcut[0]
    i+=1
    model[i] = args.sigma[0]
    i+=1
    model[i] = args.radial[0]
    i+=1
    model[i] = args.radial[1]
    i+=1
    model[i] = args.radial[2]
    i+=1
    model[i] = numcen
    if (numcen>0):
        for j in xrange(numcen):
            i+=1
            model[i] = atomic_numbers[args.centres[j]]
    i+=1
    model[i] = numspec
    if (numspec>0):
        for j in xrange(numspec):
            i+=1
            model[i] = atomic_numbers[args.species[j]]

    # Save files
    output = args.outfile
    if (output.count('/')>0):
        op = output.split('/')
        folder = '/'.join(op[:len(op)-1])
        if (not os.path.isdir(folder)):
            split_folder = op[:len(op)-1]
            make_folder = ''
            for i in xrange(len(split_folder)):
                make_folder += split_folder[i] + '/'
                if (not os.path.isdir(make_folder)):
                    os.mkdir(make_folder)
    mfl = output + '.mdl'
    hfl = output + '.hyp'
    model.tofile(mfl)
    fl = open(hfl,'w')
    if (args.hyperparam):
        print >> fl, args.hyperparam[0]
    else:
        print >> fl, ''
    fl.close()

else:
    pp = np.load(args.power)
    fp = np.load(args.sparse + '_fps.npy')
    am = np.load(args.sparse + '_Amat.npy')
    weights = np.load(args.weights)
    wt = weights[4]
    if (len(weights)>5):
        mv = np.load(args.weights)[5]
    else:
        mv = 0.0
    # Repeat for scalar
    p0 = np.load(args.power0)
    f0 = np.load(args.sparse0 + '_fps.npy')
    a0 = np.load(args.sparse0 + '_Amat.npy')
    # Make model array
    nmol = np.shape(pp)[0]
    ncut = np.shape(pp)[-1]
    ncut0 = np.shape(p0)[-1]
    if (len(np.shape(pp))==3):
        degen = 1
    else:
        degen = np.shape(pp)[-2]
    ln = 7 + np.size(pp) + np.size(p0) + 1 + np.size(fp) + 2*np.size(am) + 1 + np.size(f0) + 2*np.size(a0) + 1 + np.size(wt) + 14 + 1 + numcen + 1 + numspec
    model = np.zeros(ln,float)
    p1 = np.zeros((nmol,1,degen,ncut),float)
    if (degen==1):
        for j in xrange(nmol):
            for k in xrange(ncut):
                p1[j,0,0,k] = pp[j,0,k]
    else:
        for j in xrange(nmol):
            for k in xrange(degen):
                for l in xrange(ncut):
                    p1[j,0,k,l] = pp[j,0,k,l]
    model[0] = 1.0
    model[1] = args.lambdaval
    model[2] = args.zeta
    if (args.periodic):
        model[3] = 1.0
    else:
        model[3] = 0.0
    model[4] = nmol
    model[5] = ncut
    model[6] = ncut0
    i = 6
    for j in xrange(nmol):
        for k in xrange(degen):
            for l in xrange(ncut):
                i+=1
                model[i] = p1[j,0,k,l]
    for j in xrange(nmol):
        for k in xrange(ncut0):
            i+=1
            model[i] = p0[j,0,k]
    i+=1
    model[i] = ncut
    for j in xrange(ncut):
        i+=1
        model[i] = fp[j]
    for j in xrange(ncut):
        for k in xrange(ncut):
            i+=1
            model[i] = np.real(am[j,k])
            i+=1
            model[i] = np.imag(am[j,k])
    i+=1
    model[i] = ncut0
    for j in xrange(ncut0):
        i+=1
        model[i] = f0[j]
    for j in xrange(ncut0):
        for k in xrange(ncut0):
            i+=1
            model[i] = np.real(a0[j,k])
            i+=1
            model[i] = np.imag(a0[j,k])
    i+=1
    model[i] = mv
    for j in xrange(nmol*degen):
        i+=1
        model[i] = wt[j]
    i+=1
    model[i] = float(args.nmax[0])
    i+=1
    model[i] = float(args.lmax[0])
    i+=1
    model[i] = args.rcut[0]
    i+=1
    model[i] = args.sigma[0]
    i+=1
    model[i] = args.radial[0]
    i+=1
    model[i] = args.radial[1]
    i+=1
    model[i] = args.radial[2]
    i+=1
    if (len(args.nmax)>1):
        model[i] = float(args.nmax[1])
    else:
        model[i] = float(args.nmax[0])
    i+=1
    if (len(args.lmax)>1):
        model[i] = float(args.lmax[1])
    else:
        model[i] = float(args.lmax[0])
    i+=1
    if (len(args.rcut)>1):
        model[i] = args.rcut[1]
    else:
        model[i] = args.rcut[0]
    i+=1
    if (len(args.sigma)>1):
        model[i] = args.sigma[1]
    else:
        model[i] = args.sigma[0]
    i+=1
    if (len(args.radial)>3):
        model[i] = args.radial[3]
    else:
        model[i] = args.radial[0]
    i+=1
    if (len(args.radial)>3):
        model[i] = args.radial[4]
    else:
        model[i] = args.radial[1]
    i+=1
    if (len(args.radial)>3):
        model[i] = args.radial[5]
    else:
        model[i] = args.radial[2]
    i+=1
    model[i] = numcen
    if (numcen>0):
        for j in xrange(numcen):
            i+=1
            model[i] = atomic_numbers[args.centres[j]]
    i+=1
    model[i] = numspec
    if (numspec>0):
        for j in xrange(numspec):
            i+=1
            model[i] = atomic_numbers[args.species[j]]

    # Save files
    output = args.outfile
    if (output.count('/')>0):
        op = output.split('/')
        folder = '/'.join(op[:len(op)-1])
        if (not os.path.isdir(folder)):
            split_folder = op[:len(op)-1]
            make_folder = ''
            for i in xrange(len(split_folder)):
                make_folder += split_folder[i] + '/'
                if (not os.path.isdir(make_folder)):
                    os.mkdir(make_folder)
    mfl = output + '.mdl'
    hfl = output + '.hyp'
    model.tofile(mfl)
    fl = open(hfl,'w')
    if (args.hyperparam):
        print >> fl, args.hyperparam[0]
    else:
        print >> fl, ''
    fl.close()
