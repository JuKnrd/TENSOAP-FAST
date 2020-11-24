This folder contains the alpha-H2O and mu-H2O models from the paper:

V. Kapil, D. M. Wilkins, J. Lan, M. Ceriotti, "Inexpensive modeling of quantum dynamics using path integral generalized Langevin equation thermostats", J. Chem. Phys. 152, 124104 (2020)

To use these models, see instructions below:

1. mu-H2O:
	# Get lambda=1 prediction
	${TENSOAP}/bin/sagpr_apply_process -f ${FILE}.xyz -m ${TENSOAP}/models/bulk-water/mu-H2O.mdl -o dipole_mu1.out
	# Convert to prediction of polarization
	${TENSOAP}/tools/spherical_to_cartesian.sh 1 dipole_mu1.out dipole_prediction.out

2. alpha-H2O:
