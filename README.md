# A tool to assess a quality of RNAP predictions

## Features:

The IINF tool is designed to evaluate *in silico* predictions of RNA-protein
3D structures by quantifying hydrogen bonds at RNA-protein binding sites.
It compares the predicted model of the RNA-protein assembly with a
reference structure and computes the Intermolecular Interaction Network
Fidelity (I-INF), a normalized similarity measure with values ranging from
0 to 1. IINF relies on two external packages, rna-tools and HBPLUS, and
operates through the following steps:
 * Preprocess the model(s) to ensure consistency with the reference structure (target),
 * Analyze both target and model(s), rename, renumber and delete parts of chains if necessary (using the rna-tools),
 * Run HBLUS on all PDB files considered to identify hydrogen bonds,
 * Compare hydrogen bonds found in the model within the context of the reference structure to compute model-target interface similarity.

## Output:

The tool allows users to compute the I-INF (Intermolecular Interaction Network Fidelity) score between two or more PDB files.

The I-INF score (a range between 0.0 and 1.0 - the higher value the better) is the measure of similarity between two molecule interfaces considered.
This is given by the formula:

$$
\text{inf} = \sqrt{\frac{tp}{tp + fp} \cdot \frac{tp}{tp + fn}}
$$

Where:

$$\text{tp}$$
is the number of pairs (binding residues) that are found both in the target and the prediction files

$$\text{fn}$$
is the number of pairs (binding residues) that are found both in the target but not in the prediction file

$$\text{fp}$$
is the number of pairs (binding residues) that are found both in the prediction but not in the target file

## Installation:

The script ```compare.py``` is viable to be run in an IDE or directly executed.

Beforehand, it requires installing rna-tools and HBPLUS.
 The rna-tools can be installed using the following installation guide:
 https://rna-tools.readthedocs.io/en/latest/install.html
 Next, its path should be added to the PATH environment variable in Linux (export PATH="$PATH:<path_to_rna_tools_scripts>").

 HBPLUS can be installed from the website:
 https://www.ebi.ac.uk/thornton-srv/software/HBPLUS
 where installation instructions are provided. 
 Next, its path should be added to the PATH environment variable in Linux (export PATH="$PATH:<path_to_hbplus>")

Using Conda:

```conda create -n iinf-env python=3.12```

```conda activate iinf-env```

```pip install -r requirements.txt```

To set environment variables, the ```set-exports.sh``` script can be used. Make sure to change the ```TODO``` to your HBPLUS and rna-tools installation paths. Run using:

```source ./set-exports.sh```

## Manual:

To run:

```python3 compare.py --target_path "target.pdb" --model_path "prediction.pdb"```

The prediction.pdb can be replaced by a directory containing multiple PDB files.

Parameters:

```--target_path``` - path to the reference structure (e.g., target.pdb).

```--model_path``` - path to the 3D model (s) (e.g., model.pdb).

```-a, --adjust_inf``` - adjust the I-INF value, if the amount of residues differs between target and prediction, the resultant I-INF score is multiplied by a fraction residues_no_in_predicition/residues_no_in_target.

```-r, --renumber_structures``` - renumber all chains in both target and model(s), so that residue numbering is consistant for comparison.

```-o, --own_mapping``` - chain identifier mapping will not be generated based on order of chains in files.

```-c, --custom_alignement``` - use user own (custom) alignement format for renumbering. Can be used without -r. Example below.

```-d, --custom_removal``` - use user own removal template for deleting parts of chains. Example below.

Alignement format example (filename|alignement;[next]):

```
7qde#2|C:103-113>A:103-113;7qde#10|A:1-11>A:103-113
```

Note, that if there are negative residue numbers to renumber, use "." insted of "-".

No custom alignment introduced eliminates residue renumbering.

Removal format example (filename|removal;[next]):

```
7qde#1|A:112-113,B:281-282;7qde#2|C:112-113,B:281-282;7qde#10|A:10-11,B:281-282
```

## Usage scenarios:

Usage scenarios are provided in the ```run_*_test.sh``` files.

## Dependancies:

 * rna-tools ([rna-tools.readthedocs.io](https://rna-tools.readthedocs.io/en/latest/))
 * hbplus (https://www.ebi.ac.uk/thornton-srv/software/HBPLUS)

Both tools should be added to PATH or otherwise accesible.

## Tests:

A) To execute very simple case, run ```run_basic_test.sh```. Expected result for '7qde' subdirectory is presented below:

```
model,score
7qde#10,0.712
7qde#8,0.650
7qde#2,0.632
7qde#4,0.624
7qde#7,0.562
7qde#3,0.516
7qde#6,0.505
7qde#9,0.500
7qde#5,0.485

```

B) To experiment with custom alignment, run ```run_custom_alignment_test.sh```. Expected result for '7qde_ca' subdirectory is presented below:

```
model,score
7qde#10,0.712
7qde#2,0.632
```

C) To experiment with both custom alignment and custom residue deletion at once, run ```run_complex_custom_alignment_test.sh```. Expected result for '7qde_ca' subdirectory is presented below:

```
model,score
7qde#10,0.816
7qde#2,0.730
```

D) To run the code on a larger dataset, use ```run_large_dataset_test.sh```. Expected result, from the combined.csv file, is presented below:

```
model,I_INF
1C0A,0.456
1DFU,0.434
1E8O,0.258
1F7U,0.417
1F7Y,0.706
1FFY,0.582
1G1X,0.704
1GAX,0.655
1H3E,0.000
1H4S,0.690
1HQ1,0.966
1IL2,0.353
1J1U,0.882
1JBS,0.000
1JID,0.340
1K8W,0.242
1KOG,0.514
1LNG,0.533
1MMS,0.308
1N78,0.531
1OOA,0.000
1Q2R,0.478
1QTQ,0.765
1R3E,0.764
1R9F,0.274
1RC7,0.094
1S03,0.492
1SER,0.000
1SJ3,0.719
1T0K,0.000
1U0B,0.000
1UN6,0.000
1YVP,0.806
2AKE,0.641
2ANR,1.000
2AZ0,0.500
2BH2,0.520
2BTE,0.167
2CSX,0.516
2CZJ,0.105
2DU3,0.236
2FK6,0.000
2FMT,0.169
2GJW,0.612
2HW8,0.261
2IPY,0.070
2NUG,0.406
2QUX,0.788
2R8S,0.106
2RFK,0.441
2UWM,0.245
2V3C,0.000
2VPL,0.250
2XDB,0.886
2ZKO,0.338
2ZM5,0.785
2ZNI,0.317
2ZUE,0.936
2ZZM,0.557
3ADD,0.263
3CIY,0.183
3DD2,0.566
3EPH,0.842
3FOZ,0.850
3FTF,0.316
3HHZ,0.309
3HL2,0.000
3LRR,0.408
3LWR,0.793
3MOJ,0.667
3OL9,0.930
3OVB,0.878
```

E) To experiment with correcting a file to be compared, run ```run_custom_rna-protein_test.sh```. Expected result from the 'complex' directory is presented below:

```
model,score
model,0.434
```
F) To experiment with negative numbering correction, run ```run_negative_renumbering_test.sh```. Expected result from the negative directory is presented below:

```
model,score
model,0.886
```

G) To experiment with scaling the I-INF, run ```run_adjusted_score_test.sh```. Expected result from the 7qde_as is presented below:

```
model,score
7qde#2,0.625
7qde#5,0.500
```
