This example predicts the response properties for bulk water.

TNS is the base folder of TENSOAP-FAST:

${TNS}/bin/sagpr_apply -m ${TNS}/example/models/bulk-water/{mu-H2O.mdl,alpha-H2O-0.mdl,alpha-H2O-2.mdl} -f frames.xyz -o prediction-L{1,0,2}.out

To apply each model separately:

${TNS}/bin/sagpr_apply -m ${TNS}/example/models/bulk-water/mu-H2O.mdl      -f frames.xyz -o prediction-L1.out
${TNS}/bin/sagpr_apply -m ${TNS}/example/models/bulk-water/alpha-H2O-0.mdl -f frames.xyz -o prediction-L0.out
${TNS}/bin/sagpr_apply -m ${TNS}/example/models/bulk-water/alpha-H2O-2.mdl -f frames.xyz -o prediction-L2.out
