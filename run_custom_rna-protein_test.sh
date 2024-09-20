#!/bin/bash
find . -name "*.csv" -exec rm {} \;
python3 compare.py --target_path "examples/complex/target.pdb" --model_path "examples/complex/model.pdb" -d "target|C:110-110;model|N:68-68" -r -c "model|P:1-94>V:1-94,M:91-109>C:91-109,N:69-86>B:69-86"
