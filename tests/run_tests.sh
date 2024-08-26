#!/bin/bash

output1=$(python3 casp.py -a "target-7YR6.pdb" "test" -r)
echo "$output1" | tail -n 1

output2=$(python3 casp.py -a "target-7YR6.pdb" "mR1190TS035_1od.pdb" -r)
echo "$output2" | tail -n 1

output3=$(python3 casp.py -a "target-7YR6.pdb" "mR1190TS035_1od.pdb" -r -c --model_renum "A:1-54>A:1-54,A:56-56>A:55-55,D:1-54>D:1-54,D:56-56>D:55-55")
echo "$output3" | tail -n 1

