#!/bin/bash
find . -name "*.csv" -exec rm {} \;
python3 compare.py -r --target_path "examples/2lxc/2lxc#1.pdb" --model_path "examples/2lxc/conformations"
