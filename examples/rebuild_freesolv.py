"""
Use the gaff2xml wrappers of OpenEye and Antechamber to rebuild input files
for the FreeSolv database.
"""

import pickle
import os
import glob
import gaff2xml

FREESOLV_PATH = "/home/kyleb/src/choderalab/FreeSolv/"
FREESOLV_FILENAME = os.path.join(FREESOLV_PATH, "primary-data/primary-data.pickle")

database = pickle.load(open(FREESOLV_FILENAME))

for (key, entry) in database.items():
    print "Processing molecule %s ..." % (key)
    
    tripos_filename = os.path.join('tripos_mol2', key + '.mol2')
    gaff_mol2_filename = os.path.join("./gaff_mol2/", "%s.prmtop" % key)
    frcmod_filename = os.path.join("./gaff_mol2/", "%s.frcmod" % key)
    gaff_mol2_filename = os.path.join("./gaff_mol2/", "%s.mol2" % key)
    frcmod_filename = os.path.join("./gaff_mol2/", "%s.frcmod" % key)
    prmtop_filename = os.path.join("./gaff_mol2/", "%s.prmtop" % key)
    inpcrd_filename = os.path.join("./gaff_mol2/", "%s.inpcrd" % key)

    molecule = gaff2xml.openeye.smiles_to_oemol(entry['smiles'])
    charged = gaff2xml.openeye.get_charges(molecule)
    gaff2xml.openeye.molecule_to_mol2(charged, tripos_filename)

    _, _ = gaff2xml.utils.run_antechamber(key, tripos_filename, charge_method=None, gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename)
    gaff2xml.utils.run_tleap(key, gaff_mol2_filename, frcmod_filename, prmtop_filename=prmtop_filename, inpcrd_filename=inpcrd_filename)
