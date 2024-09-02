# A tool to assess a quality of RNAP predictions

## Features:

A tool that allows users to assess the interface accuracy of in-silico 3D models of RNAPs using INF-based score (Interaction Network Fidelity]. It uses the following packages: rna-tools and HBPlus.
To do that, it executes the following steps:
 * Preprocess the model(s) to ensure consistency with the reference structure (target),
 * Analyze both target and model(s), rename and renumber chains if necessary (using the rna-tools),
 * Run HBPlus on all PDB files considered to identify hydrogen bonds,
 * Compare hydrogen bonds found in the model within the context of the reference structure to compute model-target interface similarity.

## Output:

The tool allows users to compute the INF (Interaction Network Fidelity) score between two or more PDB files.

The INF score (a range between 0.0 and 1.0 - the higher value the better) is the measure of similarity between two molecule interfaces considered.
This is given by the formula:

$$
\text{inf} = \sqrt{\frac{tp}{tp + fp} \cdot \frac{tp}{tp + fn}}
$$

Where:

$$\text{tp}$$
is the number of pairs (bonds) that are found both in the target and the prediction files

$$\text{fn}$$
is the number of pairs (bonds) that are found both in the target but not in the prediction file

$$\text{fp}$$
is the number of pairs (bonds) that are found both in the prediction but not in the target file

## Installation:
The script ```compare.py``` is viable to be run in an IDE or directly executed.

Beforehand, it requires installing rna-tools and HBPlus.
 rna-tools can be installed using the following installation guide:
 https://rna-tools.readthedocs.io/en/latest/install.html
 Next, its path should be added to the PATH environment variable in Linux (export PATH="$PATH:<path_to_rna_tools_scripts>").

 HBPlus can be installed from the website:
 https://www.ebi.ac.uk/thornton-srv/software/HBPLUS
 where installation instructions are provided. 
 Next, its path should be added to the PATH environment variable in Linux (export PATH="$PATH:<path_to_hbplus>")

Using Conda:

```conda create -n iinf-env python=3.12```

```conda activate iinf-env```

```pip install -r requirements.txt```

To set environment variables, the ```set-exports.sh``` script can be used. Make sure to change the ```TODO``` to your HBPlus and rna-tools installation paths. Run using:

```source ./set-exports.sh```

## Manual:

To run:

```python3 compare.py --target_path "target.pdb" --model_path "prediction.pdb"```

The prediction.pdb can be replaced by a directory containing multiple PDB files.

Parameters:

```--target_path``` - path to the reference structure (e.g., target.pdb).

```--model_path``` - path to the 3D model (s) (e.g., model.pdb).

```-a, --adjust_inf``` - adjust the INF value, if the amount of residues differs between target and prediction, the resultant INF score is multiplied by a fraction residues_no_in_predicition/residues_no_in_target.

```-r, --renumber_structures``` - renumber all chains in both target and model(s), so that residue numbering is consistant for comparison.

```-c, --custom_alignement``` - use user own (custom) alignement format for renumbering. Example below.

```--target_renum``` - custom target alignement.

```--model_renum``` - custom model(s) alignement.

Alignement format example (protein chain A):

```
A:1-54>A:1-54,A:56-56>A:55-55
```
## Example Usage:
Example usage is provided in the ```run_tests.sh``` file.

## Dependancies:
 * rna-tools ([rna-tools.readthedocs.io](https://rna-tools.readthedocs.io/en/latest/))
 * hbplus (https://www.ebi.ac.uk/thornton-srv/software/HBPLUS)

Both tools should be added to PATH or otherwise accesible.

## Tests:
To test, run ```run_tests.sh```. Expected result:
```
model,score
7qde#10,0.712
7qde#2,0.632
7qde#3,0.516
7qde#4,0.624
7qde#5,0.485
7qde#6,0.505
7qde#7,0.562
7qde#8,0.650
7qde#9,0.500

```
