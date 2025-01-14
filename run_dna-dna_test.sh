#!/bin/bash
find . -name "*.csv" -exec rm {} \;
python3 compare.py --target_path "examples/1lwa/1lwa#1.pdb" --model_path "examples/1lwa/conformations"
