#!/bin/bash

set -ueo pipefail
SAMPLES="AAef1
AAef3
AAem1
AAem3
AAf1
AAf2
AAf3
AAf4
AAf5
AAf6
AAf7
AAm1
AAm2
AAm3
AAm4
AAm5
AAm6
AAm7
"

count=1
for i in $SAMPLES
do
echo $i > input_$count
count=$((count+1))
done