#!/usr/bin/env python
import sys
from openmoltools.utils import run_antechamber

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("""Usage: generate_example_data.py ligand_name
Note: this should be run in the openmoltools/chemicals/ligand_name directory.
""")
    else:
        molecule_name, mol2_filename = sys.argv[1:]
        run_antechamber(molecule_name, mol2_filename, charge_method="bcc")
