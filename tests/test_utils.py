from unittest import skipIf
import logging

logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")

from mdtraj.testing import eq
from gaff2xml import utils

def test_enter_temp_directory():
    with utils.enter_temp_directory():
        pass

def test_parse_ligand_filename():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    name, ext = utils.parse_ligand_filename(input_filename)
    
    eq(name, "sustiva")
    eq(ext, ".mol2")

def test_run_antechamber():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        gaff_mol2_filename, frcmod_filename = utils.run_antechamber(molecule_name, input_filename, charge_method=None)

def test_run_tleap():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        gaff_mol2_filename, frcmod_filename = utils.run_antechamber(molecule_name, input_filename, charge_method=None)
        prmtop, inpcrd = utils.run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename)

def test_run_test_molecule():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        utils.test_molecule(molecule_name, input_filename)

def test_run_antechamber_charges():
    molecule_name = "acetate"
    input_filename = utils.get_data_filename("chemicals/acetate/acetate.mol2")
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        gaff_mol2_filename, frcmod_filename = utils.run_antechamber(molecule_name, input_filename, charge_method=None, net_charge=-1)