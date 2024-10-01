#!/bin/bash
find . -name "*.csv" -exec rm {} \;
python3 compare.py --target_path "examples/negative/target.pdb" --model_path "examples/negative/model.pdb" -r -c "target|G:-3.32>G:1.36;model|G:-3.32>G:1.36"  
