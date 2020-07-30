#!/bin/bash

# Get error bars in pure committee prediction

export skip=${2}
if [ "${skip}" == "" ];then export skip=1;fi

cat ${1} | awk '{m=n=0;for (i=1;i<=NF;i++){n++;m+=$i};m/=n;l=0;for (i=1;i<=NF;i++){l+=($i-m)**2};l/=(n-1);printf "%.16f %.16f\n",m,l**0.5}'
