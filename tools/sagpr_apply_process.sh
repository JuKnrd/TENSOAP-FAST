#!/bin/bash

RUNDIR=$(dirname $(readlink -f "$0"))/../bin

${RUNDIR}/sagpr_apply "${@}"

fname=""
model=""
oname=""

for arg in $(seq 1 $#);do
	if [ "${!arg}" == "-f" ];then arg1=$((arg+1));fname=${!arg1};fi
	if [ "${!arg}" == "-m" ];then arg1=$((arg+1));model=${!arg1};fi
	if [ "${!arg}" == "-o" ];then arg1=$((arg+1));oname=${!arg1};fi
done

# Multiply by number of atoms

cat ${fname} | awk '{if (NF==1){print}}' > natoms.txt

if [ $(wc -l natoms.txt | awk '{print $1}') -ne $(wc -l ${oname} | awk '{print $1}') ];then echo "ERROR: natoms.txt and ${oname} have different numbers of lines";fi

paste ${oname} natoms.txt | awk '{for (i=1;i<NF;i++){printf "%f ",$i*$NF};printf "\n"}' > ${oname}_natoms.txt;mv ${oname}_natoms.txt ${oname}

rm natoms.txt
