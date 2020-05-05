#!/bin/bash

../../bin/sagpr_apply -f test.xyz -m MODEL/MODEL0.mdl -o prediction_L0.txt
../../bin/sagpr_apply -f test.xyz -m MODEL/MODEL2.mdl -o prediction_L2.txt

echo "RESULTS SHOULD BE AS CLOSE TO ZERO AS POSSIBLE:"

paste test_L0.txt prediction_L0.txt | awk 'BEGIN{n=m=0}{n++;m+=($1-$2)**2}END{print (m/n)**0.5}'
paste test_L2.txt prediction_L2.txt | awk 'BEGIN{n=m=0}{n++;m+=($1-$6)**2 + ($2-$7)**2 + ($3-$8)**2 + ($4-$9)**2 + ($5-$10)**2}END{print (m/n)**0.5}'
