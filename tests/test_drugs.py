#!/usr/bin/env python

import openeye.oechem

from gaff2xml import utils


def test_drugs():
    # Test downloaded drug database.
    database_filename = './chemicals/drugs/Zdd.mol2.gz' # mol2 database source
    ifs = openeye.oechem.oemolistream(database_filename)
    for molecule in ifs.GetOEGraphMols():
        yield lambda : utils.test_molecule(molecule)  # Cute trick to iteratively run this test over entire database.

