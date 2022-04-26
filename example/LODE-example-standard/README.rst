This example predicts the energies for random NaCl structures (from TENSOAP repository) using LODE models.

To apply a LODE model, two XXX.mdl files are needed: the first of these is used as for standard SA-GPR models, giving the SOAP hyperparameters, weights and active set power spectrum, and the second is used to give extra LODE hyperparameters. The `-l` option to `sagpr_apply` is used to specify this extra model file. If no LODE hyperparameter file is specified, or `-l NONE` is used, the result will be a standard SA-GPR model.

TNS is the base folder of TENSOAP-FAST:

${TNS}/bin/sagpr_apply -m MODEL.mdl -l LODE.mdl -f NaCl.xyz -o prediction_L0.out
