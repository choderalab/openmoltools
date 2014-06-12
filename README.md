[![Build Status](https://travis-ci.org/choderalab/gaff2xml.svg)](https://travis-ci.org/choderalab/gaff2xml)

## gaff2xml: Converting Antechamber GAFF files to OpenMM XML format

This set of tools allows users to automate the conversion of ligand
force field parameters from Antechamber formats to OpenMM XML format.

This tool is in BETA testing: use at your own risk!


Installation:

```
python setup.py install
```

Usage:

```
mkdir sustiva
cd sustiva
cp ../chemicals/sustiva/sustiva.sdf ./
generate_example_data.py sustiva
```

This should create a file sustiva.xml, which can be used by openmm to simulate your system.  
Under the hood, `generate_example_data.py` first runs Antechamber then converts
to XML via `processAmberForceField.py`.

A boilerplate script for simulating your protein-ligand system is available in 
examples/test_protein_ligand.py

Finally, we recommend that you check agreement (energies and force parameters)
of your XML-based system as compared to creating your system via tleap.  
generate_example_data should have already created the necessary inpcrd and prmtop files.
You can evaluate agreement between the XML and prmtop files via the 
boilerplate script in examples/test_example.py
