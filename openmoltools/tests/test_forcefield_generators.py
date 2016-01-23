from nose.plugins.attrib import attr
from unittest import skipIf
import tempfile
import os
from openmoltools import utils, forcefield_generators
from simtk.openmm.app import ForceField, NoCutoff

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
    openeye_exception_message = str()
except Exception as e:
    HAVE_OE = False
    openeye_exception_message = str(e)

def createOEMolFromIUPAC(iupac_name='ibuprofen'):
    from openeye import oechem, oeiupac
    mol = oechem.OEGraphMol()
    oeiupac.OEParseIUPACName(mol, iupac_name)
    oechem.OEAssignAromaticFlags(mol, oechem.OEAroModelOpenEye)
    oechem.OEAddExplicitHydrogens(mol)
    return mol

@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_PerceiveBondOrdersExplicitHydrogens(write_pdf=False):
    molecule = createOEMolFromIUPAC('ibuprofen')
    from openmoltools.forcefield_generators import PerceiveBondOrdersExplicitHydrogens
    mols = PerceiveBondOrdersExplicitHydrogens(mol)

    if write_pdf:
        # DEBUG
        from openeye import oedepict
        options = oedepict.OEReportOptions()
        report = oedepict.OEReport(options)
        for mol in mols:
            oedepict.OEPrepareDepiction(mol)
        for mol in mols:
            cell = report.NewCell()
            opts = oedepict.OE2DMolDisplayOptions()
            disp = oedepict.OE2DMolDisplay(mol, opts)
            oedepict.OERenderMolecule(cell, disp)
        oedepict.OEWriteReport('output.pdf', report)

import unittest
class TestForceFieldGenerators(unittest.TestCase):
    @skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
    def test_generateTopologyFromOEMol():
        from openeye import oechem, oeiupac
        molecule = createOEMolFromIUPAC('ibuprofen')
        from openmoltools.forcefield_generators import generateTopologyFromOEMol
        topology = generateTopologyFromOEMol()
        assertEqual(len(topology.atoms()), len(molecule.GetAtoms()))
        assertEqual(topology.residues()[0].name, molecule.GetTitle())
        for (top_atom, mol_atom) in zip(topology.atoms(), molecule.GetAtoms()):
            assertEqual(top_atom.name, mol_atom.GetName())

@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_generateResidueTemplate():
    """
    Test GAFF residue template generation from OEMol molecules.
    """
    from openeye import oechem, oeiupac
    # TODO: Test more molecules.
    iupac_name = 'ibuprofen'
    mol = oechem.OEGraphMol()
    oeiupac.OEParseIUPACName(mol, iupac_name)
    oechem.OEAssignAromaticFlags(mol, oechem.OEAroModelOpenEye)
    oechem.OEAddExplicitHydrogens(mol)
    # Generate an ffxml residue template.
    from openmoltools.forcefield_generators import generateResidueTemplate
    [template, ffxml] = generateResidueTemplate(mol)
    # Create a ForceField object.
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml', 'gaff.xml')
    # Add the additional parameters and template to the forcefield.
    forcefield.registerResidueTemplate(template)
    forcefield.load_file(StringIO(ffxml))
    # Create a Topology from the molecule.
    from openmoltools.forcefield_generators import generateOpenMMTopology
    topology = generateOpenMMTopology(molecule)
    # Parameterize system.
    system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)

def test_gaffResidueTemplateGenerator():
    """
    Test the GAFF residue template generator.
    """

    #
    # Test where we generate parameters for only a ligand.
    #

    # Load the PDB file.
    pdb_filename = utils.get_data_filename("chemicals/proteins/T4-lysozyme-L99A-p-xylene-implicit.pdb")
    pdb = PDBFile(pdb_filename)
    # Create a ForceField object.
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml', 'gaff.xml')
    # Add the residue template generator.
    from openmoltools.forcefield_generators import gaffTemplateGenerator
    forcefield.registerTemplateGenerator(gaffTemplateGenerator)
    # Parameterize system.
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    # TODO: Test energies are finite?

if __name__ == '__main__':
    test_PerceiveBondOrdersExplicitHydrogens(write_pdf=True)
