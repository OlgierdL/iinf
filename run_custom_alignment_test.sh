#!/bin/bash
find . -name "*.csv" -exec rm {} \;
python3 compare.py --target_path "examples/7qde_ca/7qde#1.pdb" --model_path "examples/7qde_ca/conformations" -r -c "7qde#2|C:103-113>A:103-113;7qde#10|A:1-11>A:103-113"
