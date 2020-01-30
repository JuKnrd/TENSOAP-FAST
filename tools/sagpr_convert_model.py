#!/usr/bin/python

import argparse
import sys,os
import numpy as np

parser = argparse.ArgumentParser(description="Put together SOAPFAST model for prediction program")
parser.add_argument("-p",   "--power",      type=str, required=True,             help="Training power spectrum file")
parser.add_argument("-sf",  "--sparse",     type=str, required=True,             help="Sparsification file")
parser.add_argument("-w",   "--weights",    type=str, required=True,             help="Weights file")
parser.add_argument("-p0",  "--power0",     type=str, required=False,            help="Scalar training power spectrum file")
parser.add_argument("-sf0", "--sparse0",    type=str, required=False,            help="Scalar sparsification file")
parser.add_argument("-w0",  "--weights0",   type=str, required=False,            help="Weights file")
parser.add_argument("-hp",  "--hyperparam", type=str, required=False, nargs='+', help="Hyperparameter string")
parser.add_argument("-o",   "--outfile",    type=str, required=True,             help="Output file")
args = parser.parse_args()

do_scalar = False

if (args.power0!=None or args.sparse0!=None or args.weights0!=None):
    # Check that all of them are specified, otherwise exit
    if (not (args.power0!='' and args.sparse0!='' and args.weights0!='')):
        print "ERROR: all scalar files must be specified!"
        sys.exit(0)
    do_scalar = True

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
    ln = 3 + np.size(pp) + 1 + np.size(fp) + 2*np.size(am) + 1 + np.size(wt)
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
    model[1] = nmol
    model[2] = ncut
    i = 2
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
    print "SCALAR POWER SPECTRA FOR SPHERICAL KERNELS NOT IMPLEMENTED!"
    sys.exit(0)
