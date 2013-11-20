#!/usr/bin/env python
import sys
from gaff2xml.utils import run_antechamber

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("""Usage: generate_example_data.py ligand_name
Note: this should be run in the gaff2xml/chemicals/ligand_name directory.
""")
    else:
        run_antechamber(sys.argv[1])
