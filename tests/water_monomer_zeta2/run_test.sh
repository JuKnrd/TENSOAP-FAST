#!/bin/bash

../../bin/sagpr_apply -f frames.xyz -m MODEL/MODEL

paste prediction.out test_prediction.out | awk 'BEGIN{n=m=0}{n++;m+=($1-$4)**2 + ($2-$5)**2 + ($3-$6)**2}END{printf "Difference = %.16f\nThis should be as close to zero as possible\n",(m/n)**0.5}'

rm prediction.out
