This example applies a committee model for the dipole moment of water dimers. To convert a committee model to TENSOAP-FAST format, the `sagpr_convert` script is used (as in the tutorial), but with multiple weight files, one for each member of the committee. All members have the same active set. If a scaling factor alpha is needed to rescale the uncertainty (as in JCTC 15, 906 [2019]), this is specified with `-a ALPHA`.

To run this code, with TNS the base folder of TENSOAP-FAST:

${TNS}/bin/sagpr_apply -m MODEL_C.mdl -f frames.xyz -o prediction_committee.txt

The output file, `prediction_committee.txt` will have (D*NC) columns, where D is the degeneracy of spherical tensor (in this case, D=3 for a rank-1 tensor) and NC the number of members of the committee, in this case 5. The output is arranged:

PRED_1_1 PRED_1_2 PRED_1_3 PRED_2_1 PRED_2_2 PRED_2_2 ... PRED_5_1 PRED_5_2 PRED_5_3

where PRED_X_Y is the Yth component of the prediction of the Xth member of the committee.

TENSOAP-FAST provides a script that can find the average of the committee predictions and their uncertainties (note that, so long as the parameter alpha was specified when converting the model, this is included in the uncertainty prediction):

::

  $ ~/source/TENSOAP-FAST/tools/find_uncertainty.sh prediction_committee.txt 3 > prediction_average.txt

The file `prediction_average.txt` has six columns, the first three of which are the average of the five committees for each component of the tensor, and the second three are the uncertainties.
