# A tool to asses precision of molecular predictions

This is a tool for assesing the correctness of in-silico models of molecular structures containing an RNA molecule and a number of protein chains by INF value. It is based on RNA Tools and HBPlus.

The script is viable to be run in an IDE or directly executed.
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

## Dependancies:
 * rna-tools ([rna-tools.readthedocs.io](https://rna-tools.readthedocs.io/en/latest/))
 * hbplus (https://www.ebi.ac.uk/thornton-srv/software/HBPLUS)

Both tools should be added to PATH or otherwise accesible.
