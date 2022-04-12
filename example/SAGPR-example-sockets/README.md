This example is based on the pes/qtip4pf/h2o-piglet example from the i-pi code (https://https://github.com/i-pi/i-pi), for bulk water with a q-TIP4P/F forcefield. The mu-H2O model is used to predict the polarization of each MD frame.

To run this code, i-pi must be used. IPI is the base folder of the i-pi code and TNS the base folder of TENSOAP-FAST:

${IPI}/bin/i-pi input.xml &> ipi.out & sleep 30
${IPI}/bin/i-pi-driver -u -h driver -m qtip4pf &> driver.out &
${TNS}/bin/sagpr_apply -m ${TNS}/example/models/bulk-water/mu-H2O.mdl -u -s sagpr 31401 -f init.xyz -o prediction.out &> sagpr.out
