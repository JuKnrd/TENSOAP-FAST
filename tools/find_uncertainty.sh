#!/bin/bash

# Get error bars in pure committee prediction

export skip=${2}
if [ "${skip}" == "" ];then export skip=1;fi

cat ${1} | awk 'BEGIN{sk=ENVIRON["skip"]}{for (j=1;j<=sk;j++){m=n=0;for (i=j;i<=NF;i+=sk){n++;m+=$i};m/=n;l=0;for (i=j;i<=NF;i+=sk){l+=($i-m)**2};l/=(n-1);printf "%.16f %.16f ",m,l**0.5};printf "\n"}' | awk '{for (i=1;i<=NF;i+=2){printf "%.16f ",$i};for (i=2;i<=NF;i+=2){printf "%.16f ",$i};printf "\n"}'
