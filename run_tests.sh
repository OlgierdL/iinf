#!/bin/bash

python3 compare.py --target_path "examples/7qde/7qde#1.pdb" --model_path "examples/7qde/conformations"

python3 casp.py --target_path "examples/7qde_ca/7qde#1.pdb" --model_path "examples/7qde_ca/conformations" -r -c "7qde#2|A:103-113>C:103-113;7qde#10|A:103-113>A:1-11"
