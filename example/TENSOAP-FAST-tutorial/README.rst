This example provides a tutorial in which a SA-GPR model is produced using the TENSOAP code, converted into TENSOAP-FAST format, and used to make predictions.

Before starting, the user should download and compile TENSOAP (https://github.com/dilkins/TENSOAP).

Step 1: Create TENSOAP model
============================

Firstly, a model for the dipole moments of water monomers is created. SFT is a variable containing the path to the base folder of TENSOAP.

::

  $ source ${SFT}/env.sh

We split the data set into a testing and a training set, and calculate the lambda-SOAP (lambda=1) power spectrum for the training set.

::

  $ python3 -c $'from ase.io import read,write;xyz=read("coords.xyz",":");import random;random.shuffle(xyz);write("train.xyz",xyz[:800]);write("test.xyz",xyz[800:])'
  $ sagpr_get_PS -f train.xyz -lm 1 -nc 600 -o PS1_train -c H O -s H O

The command `-nc 600` specifies that the power spectrum will be sparsified over spherical harmonic features, keeping the 600 furthest points. This creates four files, `PS1_train.npy`, the training power spectrum, `PS1_natoms.npy`, an array containing the number of atoms in each frame, and two files, `PS1_train_Amat.npy` and `PS1_train_fps.npy` that contain the sparsification details. These details are then used to build the testing power spectrum:

::

  $ sagpr_get_PS -f test.xyz -lm 1 -sf PS1_train -o PS1_test -c H O -s H O

The next step is to choose an active set of points to use when training out models. For this, we begin by getting an atomic power spectrum for the training set -- this is simply an array where each row corresponds to a single environment in the set, rather than each row corresponding to a molecule in the set (which is the case in `PS1_train.npy`).

::

  $ get_atomic_power_spectrum.py -f train.xyz -p PS1_train.npy -o PS1_train_atomic

We now choose the active set using furthest point sampling, and get the active power spectrum:

::

  $ do_fps.py -p PS1_train_atomic.npy -n 1000
  $ apply_fps.py -p PS1_train_atomic.npy -o PS1_train_sparse

Now, three kernels are created: K_NM (the kernel between the training set and the active set), K_MM (kernel between the active set and itself) and K_TM (kernel between testing set and active set):

::

  $ sagpr_get_kernel -z 1 -ps PS1_train.npy PS1_train_sparse.npy -s PS1_train_natoms.npy NONE -o K_NM
  $ sagpr_get_kernel -z 1 -ps PS1_train_sparse.npy -s NONE -o K_MM
  $ sagpr_get_kernel -z 1 -ps PS1_test.npy PS1_train_sparse.npy -s PS1_test_natoms.npy NONE -o K_TM

Using the first two matrices, an SA-GPR model is trained (using the Moore-Penrose pseudoinverse instead of matrix inversion):

::

  $ sagpr_train -f train.xyz -p mu -r 1 -sf K_NM.npy K_MM.npy -reg 1e-10 -sel 0 800 -m pinv

The file `weights_1.npy` contains the weight vector. We can make predictions for the testing set with:

::

  $ sagpr_prediction -r 1 -k K_TM.npy -w weights -o prediction

The file `prediction_cartesian.txt` contains the predicted dipole moments for the training set, and `prediction_L1.txt` contains the lambda=1 spherical tensor version of the predictions (related to the Cartesian vector by a simple transposition). The following one-line bash command finds the root mean squared error (RMSE) for this model:

::

  $ paste <(cat test.xyz | sed "s/\"/ /g" | awk '/Properties/{for (i=1;i<=NF;i++){if ($i=="mu="){print $(i+1),$(i+2),$(i+3)}}}') prediction_cartesian.txt | awk 'BEGIN{n=m1=m2=0}{n++;m1+=($1^2 + $2^2 + $3^2);m2+=($1-$4)^2 + ($2-$5)^2 + ($3-$6)^2}END{print "RMSE=",(m2/n)^0.5;print "Intrinsic deviation=",(m1/n)^0.5;print "Relative RMSE=",100*(m2/m1)^0.5,"%"}'

This gave a relative RMSE of 0.6% in testing.

Step 2: Convert to TENSOAP-FAST model
=====================================

The next step is to create a file, MODEL.mdl, which contains all of the information TENSOAP-FAST needs to apply the model. TNS is a variable containing the root folder of the TENSOAP-FAST code.

::

  $ ${TNS}/bin/sagpr_convert -ps PS1_train_sparse.npy -sf PS1_train -w weights_1.npy -o MODEL -c H O -s H O -lm 1

This specifies that the active set power spectrum is contained in `PS1_train_sparse.npy`, the details for sparsification over spherical harmonic components in `PS1_train_{Amat,fps}.npy`, and the weights in `weights_1.npy`. The various extra arguments to this script are used either to specify the hyperparameters of the model.

Step 3: Apply TENSOAP-FAST model
================================

This model is now applied to recalculate the predictions for the test set.

::

  $ ${TNS}/bin/sagpr_apply -f test.xyz -m MODEL.mdl -o PREDICTION_L1.txt

This will create a file, `PREDICTION_L1.txt`, which should match as closely as possible with `prediction_L1.txt`, the output of TENSOAP. Note that here we do not check whether `PREDICTION_L1.txt` matches the testing set, but whether the two codes give the same answers.

::

  $ paste PREDICTION_L1.txt prediction_L1.txt | awk 'BEGIN{n=m=0}{n++;m+=($1-$4)^2 + ($2-$5)^2 + ($3-$6)^2}END{print "RMSE=",(m/n)^0.5}'
