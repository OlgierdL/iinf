#!/bin/bash
find . -name "*.csv" -exec rm {} \;
python3 compare.py --target_path "examples/7qde_ca/7qde#1.pdb" --model_path "examples/7qde_ca/conformations" -d "7qde#1|A:112-113,B:281-282;7qde#2|C:112-113,B:281-282;7qde#10|A:10-11,B:281-282" -r -c "7qde#2|C:103-111>A:103-111;7qde#10|A:1-9>A:103-111"
