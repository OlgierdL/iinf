#!/bin/bash
find . -name "*.csv" -exec rm {} \;
python3 compare.py --target_path "examples/2jyj/2jyj#1.pdb" --model_path "examples/2jyj/conformations"
