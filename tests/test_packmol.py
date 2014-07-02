import tempfile
import numpy as np
import mdtraj as md
from unittest import skipIf
import logging
from mdtraj.testing import eq
from gaff2xml import utils, packmol
import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
from distutils.spawn import find_executable

HAVE_RDKIT = True
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    HAVE_RDKIT = False

OBABEL_PATH = find_executable("obabel")
SKIP_SMILES = not ((HAVE_RDKIT) & (OBABEL_PATH is not None))

logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")


@skipIf(SKIP_SMILES, "Skipping testing of smiles conversion because openbabel or rdkit not found.")
def test_packmol_simulation_ternary():
    input_filename0 = utils.get_data_filename("chemicals/benzene/benzene.mol2")
    pdb_filename0 = tempfile.mktemp(suffix=".pdb")
    gaff_mol2_filename0, frcmod_filename0 = utils.run_antechamber("benzene", input_filename0)
    md.load(gaff_mol2_filename0).save(pdb_filename0)
    ffxml0 = utils.create_ffxml_file(gaff_mol2_filename0, frcmod_filename0, override_mol2_residue_name="benzene")

    input_filename1 = utils.get_data_filename("chemicals/cyclopropane/cyclopropane.mol2")
    pdb_filename1 = tempfile.mktemp(suffix=".pdb")
    gaff_mol2_filename1, frcmod_filename1 = utils.run_antechamber("cyclopropane", input_filename1)
    md.load(gaff_mol2_filename1).save(pdb_filename1)
    ffxml1 = utils.create_ffxml_file(gaff_mol2_filename1, frcmod_filename1, override_mol2_residue_name="cyclopropane")

    input_filename2 = utils.get_data_filename("chemicals/propene/propene.mol2")
    pdb_filename2 = tempfile.mktemp(suffix=".pdb")
    gaff_mol2_filename2, frcmod_filename2 = utils.run_antechamber("propene", input_filename2)
    md.load(gaff_mol2_filename2).save(pdb_filename2)
    ffxml2 = utils.create_ffxml_file(gaff_mol2_filename2, frcmod_filename2, override_mol2_residue_name="propene")

    trj = packmol.pack_box([pdb_filename0, pdb_filename1, pdb_filename2], [6, 11, 23])

    xyz = trj.openmm_positions(0)
    top = trj.top.to_openmm()
    top.setUnitCellDimensions(mm.Vec3(*trj.unitcell_lengths[0])*u.nanometer)

    forcefield = app.ForceField(ffxml0, ffxml1, ffxml2)

    temperature = 300 * u.kelvin
    friction = 0.1 / u.picosecond
    timestep = 0.1 * u.femtosecond

    system = forcefield.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=None)

    integrator = mm.LangevinIntegrator(temperature, friction, timestep)

    simulation = app.Simulation(top, system, integrator)
    simulation.context.setPositions(xyz)

    simulation.step(25)
