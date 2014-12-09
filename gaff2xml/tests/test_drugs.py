from unittest import skipIf
from gaff2xml import utils
import os

try:
    oechem = utils.import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oequacpac = utils.import_("openeye.oequacpac")
    if not oequacpac.OEQuacPacIsLicensed(): raise(ImportError("Need License for oequacpac!"))
    oeiupac = utils.import_("openeye.oeiupac")
    if not oeiupac.OEIUPACIsLicensed(): raise(ImportError("Need License for OEOmega!"))        
    oeomega = utils.import_("openeye.oeomega")
    if not oeomega.OEOmegaIsLicensed(): raise(ImportError("Need License for OEOmega!"))    
    HAVE_OE = True
except:
    HAVE_OE = False


@skipIf(not HAVE_OE, "Cannot run test_drugs() module without OpenEye tools.")
def test_drugs():
    import openeye.oechem
    database_filename = utils.get_data_filename("chemicals/drugs/Zdd.mol2.gz")
    ifs = openeye.oechem.oemolistream(database_filename)
    for molecule in ifs.GetOEGraphMols():
        with utils.enter_temp_directory():
            molecule_name, tripos_mol2_filename = utils.molecule_to_mol2(molecule)
            yield utils.tag_description(lambda : utils.test_molecule(molecule_name, tripos_mol2_filename), "Testing drugs %s" % molecule_name)

@skipIf(not HAVE_OE, "Cannot test test_drug() without OpenEye tools.")
def test_drug():
    import openeye.oechem
    database_filename = utils.get_data_filename("chemicals/drugs/Zdd.mol2.gz")
    ifs = openeye.oechem.oemolistream(database_filename)
    for molecule in ifs.GetOEGraphMols():
        with utils.enter_temp_directory():
            molecule_name, tripos_mol2_filename = utils.molecule_to_mol2(molecule)
            yield utils.tag_description(lambda : utils.test_molecule(molecule_name, tripos_mol2_filename), "Testing drugs %s" % molecule_name)
        break
