#!/usr/bin/env python
import os
import tempfile
import sys

def run_antechamber(ligand_name):
    cmd = "obabel -i sdf %s.sdf -o pdb > %s.pdb" % (ligand_name, ligand_name)
    os.system(cmd)
    cmd = "antechamber -i %s.pdb -fi pdb -o %s.mol2 -fo mol2 -c bcc -s 2" % (ligand_name, ligand_name)
    os.system(cmd)
    cmd = "parmchk -i %s.mol2 -f mol2 -o %s.frcmod" % (ligand_name, ligand_name)
    os.system(cmd)
    cmd = "python ../../scripts/processAmberForceField.py ~/src/amber12/dat/leap/parm/gaff.dat ./%s.mol2 ./%s.frcmod > %s.xml" % (ligand_name, ligand_name, ligand_name)
    os.system(cmd)

    tleap_input = """
source leaprc.ff99SB
source leaprc.gaff
LIG = loadmol2 %s.mol2 
check LIG
loadamberparams %s.frcmod
saveoff LIG %s.lib
saveamberparm LIG %s.prmtop %s.inpcrd
quit

""" % (ligand_name, ligand_name, ligand_name, ligand_name, ligand_name)

    file_handle = tempfile.NamedTemporaryFile()
    file_handle.writelines(tleap_input)
    file_handle.flush()

    cmd = "tleap -f %s " % file_handle.name
    os.system(cmd)

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("""Usage: generate_example_data.py ligand_name
Note: this should be run in the gaff2xml/chemicals/ligand_name directory.
""")
    else:
        run_antechamber(sys.argv[1])
