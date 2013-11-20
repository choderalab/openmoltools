#!/usr/bin/env python

from unittest import skipIf

import openeye.oechem

from gaff2xml import utils

def test_parse_ligand_filename():
    pass
    
def test_run_antechamber():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    gaff_mol2_filename, frcmod_filename = utils.run_antechamber(molecule_name, input_filename, charge_method=None)
