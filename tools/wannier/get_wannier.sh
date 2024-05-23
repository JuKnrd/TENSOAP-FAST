#!/bin/bash

RUNDIR=$(dirname $(readlink -f "$0"))

MODEL=""
INFILE=""
VERBOSE=""
OUTFILE=""
for arg in $(seq 1 $#);do
 if [ "${!arg}" == "-m" ];then arg1=$((arg+1));MODEL=${!arg1};fi
 if [ "${!arg}" == "-f" ];then arg1=$((arg+1));INFILE=${!arg1};fi
 if [ "${!arg}" == "-v" ];then VERBOSE="-v";fi
 if [ "${!arg}" == "-o" ];then arg1=$((arg+1));OUTFILE=${!arg1};fi
done

# Get raw outputs
${RUNDIR}/../../bin/sagpr_apply -m ${MODEL} -f ${INFILE} -o TFAST_PRED.txt -a ${VERBOSE}

# Get only outputs for oxygen atoms
cat ${INFILE} | awk '{if (NF==1){print "NATM",$1}}' > NATOM.txt
cat NATOM.txt TFAST_PRED.txt | awk '/NATM/{n++;natom[n]=$2}/#/{m++}/O/{print $4*natom[m],$2*natom[m],$3*natom[m]}' > TFAST_CART_PRED.txt

#TODO: More general version needed for other atom types!

${RUNDIR}/insert-pred-wannier.py -f ${INFILE} -w TFAST_CART_PRED.txt -o PRED_WANNIER.xyz

#TODO: Need to tell the script what the total charge is that we expect!

${RUNDIR}/get-wannier-polarization.py -f PRED_WANNIER.xyz -o ${OUTFILE}

rm TFAST_PRED.txt NATOM.txt TFAST_CART_PRED.txt PRED_WANNIER.xyz
