# TENSOAP with FORTRAN

This code allows the use of models trained using `TENSOAP`/`SOAPFAST`; the key parts of this code needed to apply a symmetry-adapted Gaussian process regression (SA-GPR) model are re-implemented in FORTRAN, making them faster and more parallelizable.

# Requirements and Installation

I have tested this code with `gfortran`, version `7.4.0`, and `ifort`, version `18.0.5`. In order to compile this code, it should suffice to run the `configure` script in the topmost directory. This will search for the required libraries, and compile them from scratch if needed; it will then produce a `Makefile`, so that you can compile the code by running `make` in the `src` directory.

The command `make install` installs the executables and libraries. The default root folder is `/usr/local`; to use a different directory either run the configuration script with the command `SAGPR_PREFIX=/path/to/folder configure` or modify the `SAGPR_PREFIX` variable in the Makefile.

It should be noted that there is no guarantee this code is as optimized as it can be; further, if you already have the required libraries (`BLAS`, `lapack`) on your system these are recommended, rather than compiling your own version.

# Use

In order to use this code, the first thing needed is a model produced using TENSOAP. There are a couple of requirements:

1. Power spectra should be sparsified on spherical harmonic components. This means that there will be two files, with names like `PS_fps.npy` and `PS_Amat.npy` providing sparsification details.
2. Power spectra should also be sparsified over environments. There will be a single sparsified power spectra, with a name like `PS.npy`
3. Finally, there should be a weights file, with a name like `weights_0.npy`.

To create model files, run `/path/to/soapfast_fortran/bin/sagpr_convert -ps PS.npy -sf PS -w weights_0.npy -o fname`, where the `PS` in `-sf PS` is the prefix to `_fps.npy` and `_Amat.npy`. Any additional (lambda-)SOAP hyperparameters are specified via command line arguments (`/path/to/soapfast_fortran/bin/sagpr_convert -h` gives a list of these). For example, `-p -n 5 -l 3 -z 2` is for a periodic system, with `nmax=5`, `lmax=3` and `zeta=2`.

This will create the file `fname.mdl`, a binary file containing the training power spectrum, sparsification details, weights and hyperparameters.

To use this model, `sagpr_apply` is called as `/path/to/TENSOAP-FAST/bin/sagpr_apply -f file_name.xyz -m fname.mdl -o prediction.out`, which uses the model `fname.mdl` to make predictions for `file_name.xyz`, storing them in `prediction.out`. *Note: these models can only be for a single spherical component; a separate prediction must be made for each component.*

# Sockets

It is also possible to use this program in combination with i-PI (https://github.com/i-PI), which can send configurations to TENSOAP-FAST via a socket interface. An example of this is given by the contents of the `example` folder.

# Committee Models

Committee models can also be used to make predictions. To use a committee model, the `-w` flag of `sagpr_convert` should be given a *list* of weights files (e.g. `-w WEIGHTS.*.npy`), and if a calibration factor is desired, this is specified by the `-a` flag.

# Long-Range Descriptors

Later versions of TENSOAP-FAST can use the LODE descriptor for long-range information. In this case, as well as the standard SA-GPR model, an additional model file is needed containing the extra LODE hyperparameters. These are specified using the `-l` flag when running `sagpr_apply`. If the user knows the the unit cell will stay the same with every input frame given to the model, the `-fc` flag can also be used to save time.

# Examples

TO BE WRITTEN

# Maintenance

Please contact d.wilkins@qub.ac.uk with any problems.
