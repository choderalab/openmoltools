from nose.plugins.attrib import attr
from unittest import skipIf
import tempfile
import os
from openmoltools import utils, forcefield_generators
from simtk.openmm.app import ForceField, NoCutoff

def test_OEPerceiveBondOrdersExplicitHydrogens(write_pdf=False):
    from openeye import oechem, oeiupac
    iupac_name = 'ibuprofen'
    mol = oechem.OEGraphMol()
    oeiupac.OEParseIUPACName(mol, iupac_name)
    oechem.OEAssignAromaticFlags(mol, oechem.OEAroModelOpenEye)
    oechem.OEAddExplicitHydrogens(mol)
    from openmoltools.forcefield_generators import OEPerceiveBondOrdersExplicitHydrogens
    mols = OEPerceiveBondOrdersExplicitHydrogens(mol)

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
    test_OEPerceiveBondOrdersExplicitHydrogens(write_pdf=True)
