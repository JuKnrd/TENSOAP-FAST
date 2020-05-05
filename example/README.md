An example for the Partridge-Schwenke water monomer potential (adapted from the corresponding i-pi example)

To run
------

${IPI}/bin/i-pi input.xml &> ipi.out & sleep 30
${IPI}/bin/i-pi-driver -u -h driver -m pswater &> driver.out &
../bin/sagpr_apply -m MODEL/MODEL0.mdl -u -s sagpr 31401 -f init.xyz &> sagpr.out
