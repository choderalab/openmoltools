[![Build Status](https://travis-ci.org/choderalab/moltools.svg)](https://travis-ci.org/choderalab/moltools)
[![Code Health](https://landscape.io/github/choderalab/moltools/master/landscape.svg)](https://landscape.io/github/choderalab/moltools/master)

## moltools: Tools for Small Molecules, Antechamber, OpenMM, and More.

This set of tools allows users to automate various tasks related to
simulating small molecules using various molecules dynamics techniques.
It also contains several python tools for working with small molecules,
packing boxes (python wrappers for packmol), and parameterizing small
molecules.  It also contains tools for creating OpenMM XML forcefield files
for small molecules.

Our goal with this project is to have modular components that are both 
documented and well-tested.  

This tool is in BETA testing: use at your own risk!

MolTools is a meant to be used with a number of companion tools:
ParmEd, MDTraj, OpenMM, PackMol, OpenEye.


Installation:

```
python setup.py install

or

conda config --add channels http://conda.binstar.org/omnia
conda install moltools

```

To test your installation, use the following command:

```
nosetests moltools -v
```
