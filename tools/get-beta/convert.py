#!/usr/bin/env python3
import numpy as np
import argparse,os

parser = argparse.ArgumentParser(description="Convert hyperpolarizabilities")
parser.add_argument("-s", "--spherical", required=True, nargs="+", help="Name of spherical predictions")
parser.add_argument("-c", "--cartesian", required=True, help="Output file name")
args = parser.parse_args()

py_dir = os.path.dirname(os.path.abspath(__file__))

outvec = [np.loadtxt(args.spherical[0]).reshape(-1),np.loadtxt(args.spherical[1]).reshape(-1)]

CR = [np.load(py_dir + "/CR_0.npy"),np.load(py_dir + "/CR_1.npy")]
CS = np.load(py_dir + "/CS.npy")
degen = [3,7]
size = int(len(outvec[0])/degen[0])
keep_cols = [[], [[1]], [[1, 0], [1, 1], [1, 2]], [[1, 0, 1], [1, 1, 0], [1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 2, 2], [1, 2, 3]]]
keep_list = [[1, 0, 1], [1, 2, 3]]
lin_dep_list = [[0, 4, 0.8944271909999162]]
sym_list = [True, True]

full_degen = [2*keep_cols[-1][i][-1] + 1 for i in range(len(keep_cols[-1]))]

full_outvec = [np.zeros((size,full_degen[i]),dtype=complex) for i in range(len(full_degen))]

# Fill the full_outvec array where we have values; anywhere else we leave it at zero. Do the transformation back to spherical harmonics here
for i in range(len(keep_list)):
    scalfac = 1.0
    if (not sym_list[i]):
        scalfac = 1.0j
    fill_col = keep_cols[-1].index(keep_list[i])
    if (degen[i] == 1):
        full_outvec[fill_col] = outvec[i] * scalfac
    else:
        ov = np.split(outvec[i],len(outvec[i]) / degen[i])
        for j in range(len(ov)):
            full_outvec[fill_col][j] = np.dot(np.conj(CR[i]).T,ov[j]) * scalfac

# Check the list of linear dependency, and include these dependencies if present
for i in range(len(lin_dep_list)):
    full_outvec[lin_dep_list[i][1]] = full_outvec[lin_dep_list[i][0]] * lin_dep_list[i][2]

# Now concatenate these outputs to give one array
concat_vec = np.zeros((int(len(outvec[0])/degen[0]),sum(full_degen)),dtype=complex)
for i in range(int(len(outvec[0])/degen[0])):
    full_outvecs = []
    for j in range(len(full_degen)):
        if (full_degen[j]==1):
            full_outvecs.append(full_outvec[j][i])
        else:
            for k in range(len(full_outvec[j][i])):
                full_outvecs.append(full_outvec[j][i][k])
    concat_vec[i] = full_outvecs

# Transform array back to Cartesian tensor
cartesian = np.real(np.dot(concat_vec,np.conj(CS).T))

np.savetxt(args.cartesian,cartesian)
