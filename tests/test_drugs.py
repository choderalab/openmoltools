from unittest import skipIf
from gaff2xml import utils
import os

@skipIf(os.environ.get("TRAVIS", None) == 'true', "Skip testing of entire drug database on Travis.")
def test_drugs():
    import openeye.oechem
    database_filename = utils.get_data_filename("chemicals/drugs/Zdd.mol2.gz")
    ifs = openeye.oechem.oemolistream(database_filename)
    for molecule in ifs.GetOEGraphMols():
        with utils.enter_temp_directory():
            molecule_name, tripos_mol2_filename = utils.molecule_to_mol2(molecule)
            yield lambda : utils.test_molecule(molecule_name, tripos_mol2_filename)  # Cute trick to iteratively run this test over entire database.

@skipIf(os.environ.get("TRAVIS", None) == 'true', "Skip testing of entire drug database on Travis.")
def test_drug():
    import openeye.oechem
    database_filename = utils.get_data_filename("chemicals/drugs/Zdd.mol2.gz")
    ifs = openeye.oechem.oemolistream(database_filename)
    for molecule in ifs.GetOEGraphMols():
        with utils.enter_temp_directory():
            molecule_name, tripos_mol2_filename = utils.molecule_to_mol2(molecule)
            yield lambda : utils.test_molecule(molecule_name, tripos_mol2_filename)  # Cute trick to iteratively run this test over entire database.
        break
