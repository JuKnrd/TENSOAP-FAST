This example provides a tutorial in which a SA-GPR moedl with LODE descriptors is produced using the TENSOAP code, converted into TENSOAP-FAST format, and used to make predictions.

Before starting, the user should download and compile TENSOAP (https://github.com/dilkins/TENSOAP), including the `make LODE` option.

Step 1: Create TENSOAP model
============================

Firstly, a model for the energy of random NaCl snapshots is created. SFT is a variable containing the path to the base folder of TENSOAP.

::

  $ source ${SFT}/env.sh

We split the data set into a testing and training set, and calculate the SOAP power spectrum for the training set.

::

  $ python3 -c $'from ase.io import read,write;xyz=read("coords_with_energies.xyz",":");import random;random.shuffle(xyz);write("train.xyz",xyz[:1600]);write("test.xyz",xyz[1600:])'
  $ sagpr_get_PS -f train.xyz -p -ele -sg 0.3 -rc 2.0 -l 0 -n 8 -nc 400 -o PS0_train -c Na Cl -s Na Cl

The command `-nc 400` specifies that the power spectrum will be sparsified over spherical harmonic features, keeping the 400 furthest points. This creates four files, `PS0_train.npy`, the training power spectrum, `PS0_natoms.npy`, an array containing the number of atoms in each frame, and two files, `PS0_train_Amat.npy` and `PS0_train_fps.npy` that contain the sparsification details. These are then used to build the testing power spectrum:

::

  $ sagpr_get_PS -f test.xyz -p -ele -sg 0.3 -rc 2.0 -l 0 -n 8 -sf PS0_train -o PS0_test -c Na Cl -s Na Cl

The next step is to choose an active set of points to use when training out models. For this, we begin by getting an atomic power spectrum for the training set -- this is simply an array where each row corresponds to a single environment in the set, rather than each row corresponding to a frame in the set (which is the case in `PS1_train.npy`).

::

  $ get_atomic_power_spectrum.py -f train.xyz -p PS0_train.npy -o PS0_train_atomic

We now choose that active set using furthest point sampling, and get the active power spectrum:

::

  $ do_fps.py -p PS0_train_atomic.npy -n 1000
  $ apply_fps.py -p PS0_train_atomic.npy -o PS0_train_sparse

Now, three kernels are created: K_NM (the kernel between the training set and the active set), K_MM (kernel between the active set and itself) and K_TM (kernel between testing set and active set):

::

  $ sagpr_get_kernel -z 1 -ps PS0_train.npy PS0_train_sparse.npy -s PS0_train_natoms.npy NONE -o K_NM
  $ sagpr_get_kernel -z 1 -ps PS0_train_sparse.npy -s NONE -o K_MM
  $ sagpr_get_kernel -z 1 -ps PS0_test.npy PS0_train_sparse.npy -s PS0_test_natoms.npy NONE -o K_TM

Using the first two matrices, an SA-GPR model is trained (using the Moore-Penrose pseudoinverse instead of matrix inversion):

::

  $ sagpr_train -f train.xyz -p energy -r 0 -sf K_NM.npy K_MM.npy -reg 1e-6 -sel 0 1600 -m pinv -perat

The files `weights_0.npy` contains the weight vector. We can make predictions for the testing set with:

:: 

  $ sagpr_prediction -r 0 -k K_TM.npy -w weights -o prediction

The file `prediction_L0.txt` contains the predicted energies per atom for the training set. The following one-line bash command finds the root mean squared error (RMSE) for this model:

::

  $ paste <(cat test.xyz | sed "s/=/ /g" | awk '!/Properties/{if (NF==1){natom=$1}}/Properties/{for (i=1;i<=NF;i++){if ($i=="energy"){print $(i+1)/natom}}}') <(cat prediction_L0.txt) | awk 'BEGIN{n=m1=m2=0}{n++;m1+=$1^2;m2+=($1-$2)^2}END{print "RMSE=",(m2/n)^0.5;print "Intrinsic deviation=",(m1/n)^0.5;print "Relative RMSE=",100*(m2/m1)^0.5,"%"}'

This gives a relative RMSE of 11% in testing (this can be improved by choice of hyperparameters, and by not normalizing the power spectrum).

Step 2: Convert to TENSOAP-FAST model
=====================================

The next step is to create a file, MODEL.mdl, which contains the weights, training power spectrum and SOAP hyperparameters, and a file LODE.mdl containing the additional LODE hyperparameters, telling TENSOAP-FAST to carry out a LODE calcu.ation. TNS is a variable containing the root folder of the TENSOAP-FAST code.

::

  $ ${TNS}/bin/sagpr_convert -ps PS0_train_sparse.npy -sf PS0_train -w weights_0.npy -o MODEL -c Na Cl -s Na Cl -sg 0.3 -rc 2.0 -l 0 -n 8 -pr
  $ ${TNS}/bin/lode_convert -o LODE

Note that this latter file is needed even if all of the LODE hyperparameters are at their default value.

Step 3: Apply TENSOAP-FAST model
================================

This model is now applied to calculate the predictions for the test set.

::

  $ ${TNS}/bin/sagpr_apply -f test.xyz -m MODEL.mdl -l LODE.mdl -o PREDICTION_L0.txt

This will create a file, `PREDICTION_L0.txt`, which should match as closely as possible with `prediction_L0.txt`, the output of TENSOAP. Note that here we do not check whether `PREDICTION_L0.txt` matches the testing set, but whether the two codes give the same answers.

::

  $ paste PREDICTION_L0.txt prediction_L0.txt | awk 'BEGIN{n=m=0}{n++;m+=($1-$2)^2}END{print "RMSE=",(m/n)^0.5}'

Upon running this, we found an RMSE of ~1e-13, indicating that the two codes match very well -- as they should!
