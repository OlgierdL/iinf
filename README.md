# A tool to assess a quality of RNAP predictions

## Features:

The I-INF tool is designed to evaluate *in silico* predictions of RNA-protein
3D structures by quantifying hydrogen bonds at RNA-protein binding sites.
It compares the predicted model of the RNA-protein assembly with a
reference structure and computes the Intermolecular Interaction Network
Fidelity (I-INF), a normalized similarity measure with values ranging from
0 to 1. I-INF relies on two external packages, rna-tools and HBPLUS, and
operates through the following steps:
 * Preprocess the model(s) to ensure consistency with the reference structure (target),
 * Analyze both target and model(s), rename, renumber and delete parts of chains if necessary (using the rna-tools),
 * Run HBLUS on all PDB files considered to identify hydrogen bonds,
 * Compare hydrogen bonds found in the model within the context of the reference structure to compute model-target interface similarity.

Additionally, it computes a similar F1 value based on the data.


## Output:

The tool allows users to compute the I-INF (Intermolecular Interaction Network Fidelity) score between two or more PDB files.

The I-INF score (a range between 0.0 and 1.0 - the higher value the better) is the measure of similarity between two molecule interfaces considered.
This is given by the formula:

$$
\text{i-inf} = \sqrt{\frac{tp}{tp + fp} \cdot \frac{tp}{tp + fn}}
$$

The F1 score (a range between 0.0 and 1.0 - the higher value the better) is a second measure implemented:

$$
\text{f1} = \frac{2 * tp}{2 * tp + fp + fn}
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

A) To execute very simple RNAP case, run ```run_basic_test.sh```. Expected result for '7qde' subdirectory is presented below:

```
model,inf_score,f1_score
7qde#10,0.712,0.692
7qde#8,0.650,0.643
7qde#2,0.632,0.615
7qde#4,0.624,0.621
7qde#7,0.562,0.562
7qde#3,0.516,0.516
7qde#6,0.505,0.500
7qde#9,0.500,0.500
7qde#5,0.485,0.483
```

B) To execute another simple case for DNA:RNA hybrid, run ```run_basic_hybrid_test.sh```. Expected result for '1hg9' subdirectory is presented below:

```
model,inf_score,f1_score
1hg9#10,1.000,1.000
1hg9#15,0.961,0.960
```

C) To execute another simple case for protein-protein complex, run ```run_protein-protein_test.sh```. Expected results for '2lxc' subdirectory are presented below:

```
#2lxc#1_ranking.csv:
model,inf_score,f1_score
2lxc#9,0.941,0.941
2lxc#5,0.667,0.667

#2lxc#1_chains_mapping.txt:
Model' best chains mapping:
2lxc#5
target,model
A,A
B,B
C,C
2lxc#9
target,model
A,A
B,B
C,C 
```

To find the best chain mapping between the particular model and the target we enumerate all perfect and maximum matches in the bipartite graph. In the graph, nodes are chains included in the particular model as well as the target. Each edge represents two strands from the model and the target, respectively having identical sequences. The enumeration of all possible chain mappings for the particular case is presented below:

```
Iterating via all possible combinations of strands...

1) Chains mapping:
target,model
A,A
B,B
C,C
Results:
model,I-INF
2lxc#5,0.667
2lxc#9,0.941

2) Chains mapping:
target,model
A,A
B,C
C,B
Results:
model,I-INF
2lxc#5,0.243
2lxc#9,0.235

Done.
```

D) To execute another simple case for RNA-RNA complex, run ```run_rna-rna_test.sh```. Expected results for '2jyj' subdirectory are presented below:

```
#2jyj#1_ranking.csv:
model,inf_score,f1_score
2jyj#8,1.000,1.000
2jyj#2,0.968,0.968

#2jyj#1_chains_mapping.txt:
Model' best chains mapping:
2jyj#8
target,model
A,A
B,B
2jyj#2
target,model
A,A
B,B 
```

The enumeration of all possible chain mappings for the particular case is presented below:

```
Iterating via all possible combinations of strands...

1) Chains mapping:
target,model
A,A
B,B
Results:
model,I-INF
2jyj#8,1.000
2jyj#2,0.968

2) Chains mapping:
target,model
A,B
B,A
Results:
model,I-INF
2jyj#8,0.000
2jyj#2,0.000

Done.
```

E) To execute another simple case for DNA-DNA complex, run ```run_dna-dna_test.sh```. Expected result for '1lwa' subdirectory is presented below:

```
model,inf_score,f1_score
1lwa#2,1.000,1.000
1lwa#18,0.979,0.979
```

F) To experiment with custom alignment, run ```run_custom_alignment_test.sh```. Expected result for '7qde_ca' subdirectory is presented below:

```
model,inf_score,f1_score
7qde#10,0.712,0.692
7qde#2,0.632,0.615
```

G) To experiment with both custom alignment and custom residue deletion at once, run ```run_complex_custom_alignment_test.sh```. Expected result for '7qde_ca' subdirectory is presented below:

```
model,inf_score,f1_score
7qde#10,0.816,0.800
7qde#2,0.730,0.727
```

H) To run the code on a larger dataset, use ```run_large_dataset_test.sh```. Expected result, from the combined.csv file, is presented below:

```
model,I_INF,F1
1C0A,0.456,0.439
1DFU,0.434,0.364
1E8O,0.258,0.211
1F7U,0.417,0.415
1F7Y,0.706,0.706
1FFY,0.582,0.582
1G1X,0.704,0.700
1GAX,0.655,0.652
1H3E,0.000,0.000
1H4S,0.690,0.667
1HQ1,0.966,0.966
1IL2,0.353,0.346
1J1U,0.882,0.875
1JBS,0.000,0.000
1JID,0.340,0.316
1K8W,0.242,0.240
1KOG,0.514,0.480
1LNG,0.533,0.531
1MMS,0.308,0.294
1N78,0.531,0.511
1OOA,0.000,0.000
1Q2R,0.478,0.471
1QTQ,0.765,0.765
1R3E,0.764,0.764
1R9F,0.274,0.273
1RC7,0.094,0.080
1S03,0.492,0.471
1SER,0.000,0.000
1SJ3,0.719,0.710
1T0K,0.000,0.000
1U0B,0.000,0.000
1UN6,0.000,0.000
1YVP,0.806,0.800
2AKE,0.641,0.640
2ANR,1.000,1.000
2AZ0,0.500,0.500
2BH2,0.520,0.511
2BTE,0.167,0.167
2CSX,0.516,0.516
2CZJ,0.105,0.105
2DU3,0.236,0.222
2FK6,0.000,0.000
2FMT,0.169,0.148
2GJW,0.612,0.600
2HW8,0.261,0.256
2IPY,0.070,0.068
2NUG,0.406,0.370
2QUX,0.788,0.788
2R8S,0.106,0.105
2RFK,0.441,0.438
2UWM,0.245,0.240
2V3C,0.000,0.000
2VPL,0.250,0.231
2XDB,0.886,0.880
2ZKO,0.338,0.333
2ZM5,0.785,0.781
2ZNI,0.317,0.314
2ZUE,0.936,0.936
2ZZM,0.557,0.552
3ADD,0.263,0.261
3CIY,0.183,0.182
3DD2,0.566,0.526
3EPH,0.842,0.842
3FOZ,0.850,0.844
3FTF,0.316,0.308
3HHZ,0.309,0.308
3HL2,0.000,0.000
3LRR,0.408,0.400
3LWR,0.793,0.778
3MOJ,0.667,0.615
3OL9,0.930,0.930
3OVB,0.878,0.878
```

I) To experiment with correcting a file to be compared, run ```run_custom_rna-protein_test.sh```. Expected result from the 'complex' directory is presented below:

```
model,inf_score,f1_score
model,0.434,0.364
```
J) To experiment with negative numbering correction, run ```run_negative_renumbering_test.sh```. Expected result from the negative directory is presented below:

```
model,inf_score,f1_score
model,0.886,0.880
```

K) To experiment with scaling the I-INF, run ```run_adjusted_score_test.sh```. Expected result from the 7qde_as is presented below:

```
model,inf_score,f1_score
7qde#2,0.625,0.609
7qde#5,0.500,0.494
```
