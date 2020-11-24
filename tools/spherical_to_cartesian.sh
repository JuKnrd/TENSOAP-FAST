#!/bin/bash

lval=${1}
if [ ${lval} == 1 ];then
	export fname1=${2}
	export degen=3
	export ncommittee=$(head -n 1 ${fname1} | awk '{print NF/ENVIRON["degen"]}')
	for i in $(seq 1 ${ncommittee});do
		export ii=${i}
		cat ${fname1} | awk 'BEGIN{dg=ENVIRON["degen"];ii=ENVIRON["ii"]}{i1=(ii-1)*dg + 1;i2=ii*dg;for (i=i1;i<=i2;i++){printf "%f ",$i};printf "\n"}' | awk '{printf "%f %f %f\n",$3,$1,$2}' > MODEL_${fname1}_${i}
	done
	oname=${3}
	paste MODEL_${fname1}_* > ${oname}
elif [ ${lval} == 2 ];then
	fname0=${2}
	fname2=${3}
	export degen=5
	export ncommittee=$(head -n 1 ${fname0} | awk '{print NF/ENVIRON["degen"]}')
	export ncommitte2=$(head -n 1 ${fname2} | awk '{print NF/ENVIRON["degen"]}')
	if [ ${ncommittee} -ne ${ncommitte2} ];then
		echo "ERROR: different numbers of committees for two models"
		echo ${ncommittee} ${ncommitte2}
		exit
	fi
elif [ ${lval} == 0 ];then
	echo "No conversion needed"
	exit
else
	echo "This l value is not supporte"
	exit
fi
