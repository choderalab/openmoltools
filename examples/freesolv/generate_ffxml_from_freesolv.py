"""
Parameterize the FreeSolv database using gaff2xml to generate an OpenMM ForceField ffxml file from the SMARTS representation of all molecules.

"""

from unittest import skipIf
import tempfile
import os
import numpy as np
import mdtraj as md
from distutils.spawn import find_executable
import tarfile
import pickle
import os
import numpy as np

import gaff2xml


# Unpack the FreeSolv database into a temporary directory.
tar_filename = gaff2xml.utils.get_data_filename("chemicals/freesolv/freesolve_v0.3.tar.bz2")
tar = tarfile.open(tar_filename, mode="r:bz2")
tar.extract('v0.3/smiles_to_cid.pickle')
tar.close()

# Load the database file.
smiles_to_cid = pickle.load(open("./v0.3/smiles_to_cid.pickle"))
smiles_strings = smiles_to_cid.keys()

# Parameterize database using SMARTS representation.
[trajectories, ffxml] = gaff2xml.utils.smiles_to_mdtraj_ffxml(smiles_strings, base_molecule_name="FreeSolv_")

# Write out ffxml file.
outfile = open('FreeSolv.ffxml')
outfile.write(ffxml)
outfile.close()


