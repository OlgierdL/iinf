#!/bin/bash

python3 casp.py -a --target_path "7qde/7qde#1.pdb" --model_path "7qde/conformations" -r --out_name "ranking1"

python3 casp.py -a --target_path "7qde_ca/7qde#1.pdb" --model_path "7qde_ca/conformations" -r -c "7qde#2|A:103-113>C:103-113;7qde#10|A:103-113>A:1-11" --out_name "ranking2"
