"""
Test some common molecules using OpenEye tools.

"""

import openmoltools.openeye
from nose.plugins.attrib import attr
from unittest import skipIf
from openmoltools import utils
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

molecules = [
    'benzene',
    'toluene',
    'phenol',
    'catechol',
    'aspirin',
    ]

@skipIf(not HAVE_OE, "Cannot run test_common_molecules() module without OpenEye tools.")
def test_common_molecules():
    import openeye.oechem
    for molecule_name in molecules:
        molecule = openmoltools.openeye.iupac_to_oemol(molecule_name)
        molecule = openmoltools.openeye.get_charges(molecule, keep_confs=1)
        with utils.enter_temp_directory():
            tripos_mol2_filename = 'molecule.mol2'
            molecule_name, tripos_mol2_filename = utils.molecule_to_mol2(molecule, tripos_mol2_filename=tripos_mol2_filename)
            yield utils.tag_description(lambda : utils.test_molecule('molecule', tripos_mol2_filename), "Testing molecule %s" % molecule_name)
