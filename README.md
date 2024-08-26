# A tool to asses precision of molecular predictions

## Features:

This is a tool for assesing the correctness of in-silico models of molecular structures containing an RNA molecule and a number of protein chains by INF value. It is based on RNA Tools and HBPlus.
To do that, it executes the following steps:
 * Prepare the model file for usage,
 * Analyze both target and model files, rename and renumber chains if necessary (using rna-tools),
 * Run HBPlus on both files to find hydrogen bonds,
 * Compare the found molecule pairs in target and model to asses similarity by calculating the INF value.

## Output:

The tool allows users to compute the INF value between two .pdb files.

The INF value is the measure of similarity between two molecules.
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
The script ```casp_compare.py``` is viable to be run in an IDE or directly executed.

Beforehand, it requires installing rna-tools and HBPlus.
 RNA Tools:
 RNA tools can be installed using the following installation guide:
 https://rna-tools.readthedocs.io/en/latest/install.html
 It may be needed to add path to rna-tools to the PATH in Linux (export PATH="$PATH:<path_to_rna_tools_scripts>")

 HBPlus:
 HBplus can be installed from the website:
 https://www.ebi.ac.uk/thornton-srv/software/HBPLUS
 where installation instructions are provided. 
 Then, add HBPlus to PATH (export PATH="$PATH:<path_to_hbplus>")

Using Conda:

```conda create -n iinf-env python=3.12```

```conda activate iinf-env```

```pip install -r requirements.txt```

To set environment variables, the ```set-exports.sh``` script can be used. Make sure to change the ```TODO``` to your HBPlus and rna-tools installation paths. Run using:

```source ./set-exports.sh```

## Usage:

To run:

```python3 casp.py```

Then input:

```target.pdb```

```prediction.pdb``` (or ```directory of .pdb files```)

And ```Y/N``` to decide whether the INF is to adjusted for mismatched residue numbers. The default answer is no.

Giving a directory instead of a single prediciton will result in the tool comparing all .pdb files of the directory to the given targed file.

NOTE: Depending on your environment, the job may be stopped upon running hbplus. It will finish correctly when resumed.

In case of mismatched line numbering, the user will be asked whether to auto-renumber. If not, then the user will be aksed to input custom alignement.

Alignement format example (protein chain A):

```
A:1-54>A:1-54,A:56-56>A:55-55
```
## Example Usage:
Example usage is provided in the /examples folder.

## Dependancies:
 * rna-tools ([rna-tools.readthedocs.io](https://rna-tools.readthedocs.io/en/latest/))
 * hbplus (https://www.ebi.ac.uk/thornton-srv/software/HBPLUS)

Both tools should be added to PATH or otherwise accesible.

## Tests:
To test, run ```run_tests.sh```. Expected result:
```
[['R1189TS035_3om_t_t.hb2', '0.000'], ['R1190TS494_3om_t_t.hb2', '0.086'], ['R1189TS035_4om_t_t.hb2', '0.000'], ['R1189TS035_5om_t_t.hb2', '0.000']]
[['mR1190TS035_1odm_t_ren_t.hb2', '0.091']]
[['mR1190TS035_1odm_t_cus_t.hb2', '0.091']]
```
