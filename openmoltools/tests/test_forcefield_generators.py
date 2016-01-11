from nose.plugins.attrib import attr
from unittest import skipIf
import tempfile
import os
from openmoltools import utils, forcefield_generators
from simtk.openmm.app import ForceField, NoCutoff

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


