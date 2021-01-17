An example, based on the pes/qtip4pf/h2o-piglet example from the i-pi code (https://https://github.com/i-pi/i-pi) of bulk water with a q-TIP4P/F forcefield, and the mu-H2O model used to predict the polarization.

To run
------

Where IPI is the base folder of your i-pi code and TNS the base folder of TENSOAP-FAST

${IPI}/bin/i-pi input.xml &> ipi.out & sleep 30
${IPI}/bin/i-pi-driver -u -h driver -m qtip4pf &> driver.out &
${TNS}/bin/sagpr_apply -m ${TNS}/models/bulk-water/mu-H2O.mdl -u -s sagpr 31401 -f init.xyz &> sagpr.out
