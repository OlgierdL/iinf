#!/bin/bash
find . -name "*.csv" -exec rm {} \;
python3 compare.py --target_path "examples/7qde_as/7qde#1.pdb" --model_path "examples/7qde_as/conformations" -a
