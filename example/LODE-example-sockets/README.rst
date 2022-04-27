This example shows the use of a LODE model with sockets. We use the i-pi code (https://https://github.com/i-pi/i-pi) in replay mode, with a trajectory output by running the example from SAGPR-example-sockets. Here, the polarization is recalculated using LODE.

To run this code, i-pi must be used. IPI is the base folder of the i-pi code and TNS the base folder of TENSOAP-FAST:

${IPI}/bin/i-pi input.xml &> ipi.out & sleep 30
${TNS}/bin/sagpr_apply -m MODEL.mdl -u -s sagpr 31401 -f init.xyz -o prediction.out &> sagpr.out

Note that the results are not expected to be close to those of the SAGPR example; the LODE model has not (yet) been optimized in any way.
