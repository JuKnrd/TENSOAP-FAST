#!/bin/bash

lval=${1}
if [ ${lval} == 1 ];then
	export fname1=${2}
	export degen=3
	export oname=${3}
	cat ${fname1} | awk '/#/{print}!/#/{if (NF==1){print}else{print $1,$4,$2,$3}}' > ${oname}
elif [ ${lval} == 2 ];then
	fname0=${2}
	fname2=${3}
	oname=${4}
	paste ${fname0} ${fname2} | awk 'BEGIN{f2=(1./2.)**0.5;f3=(1./3.)**0.5;f6=(1./6.)**0.5;f23=(2./3.)**0.5}/#/{print "#"}!/#/{if (NF==2){print $1}else{a0=$2;a2m2=$4;a2m1=$5;a20=$6;a2p1=$7;a2p2=$8;axy=ayx=f2*a2m2;ayz=azy=f2*a2m1;axz=azx=f2*a20;axx=(-f3*a0 - f6*a2p1 + f2*a2p2);ayy=(-f3*a0 - f6*a2p1 - f2*a2p2);azz=(-f3*a0 + f23*a2p1);printf "%f %f %f %f %f %f %f %f %f\n",axx,axy,axz,ayx,ayy,ayz,azx,azy,azz}}' > ${oname}
elif [ ${lval} == 0 ];then
	echo "No conversion needed"
	exit
else
	echo "This l value is not supported"
	exit
fi
#rm MODEL_*
