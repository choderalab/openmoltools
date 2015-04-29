import tempfile
import numpy as np
import mdtraj as md
from unittest import skipIf
import logging
from mdtraj.testing import eq
from openmoltools import utils, packmol
import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm

logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")

HAVE_RDKIT = True
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    HAVE_RDKIT = False

@skipIf(not HAVE_RDKIT, "Skipping testing of packmol conversion because rdkit not found.")
@skipIf(packmol.PACKMOL_PATH is None, "Skipping testing of packmol conversion because packmol not found.")
def test_packmol_simulation_ternary():
    smiles_list = ["Cc1ccccc1", "c1ccccc1", "CC"]
    pdb_filenames = []
    ligand_trajectories, ffxml = utils.smiles_to_mdtraj_ffxml(smiles_list)

    for k, ligand_traj in enumerate(ligand_trajectories):
        pdb_filename = tempfile.mktemp(suffix=".pdb")
        ligand_traj.save(pdb_filename)
        pdb_filenames.append(pdb_filename)

    pdb_filenames = pdb_filenames[0:2] + [ligand_traj]  # Test passing BOTH pdb filenames and trajectories as input.
    trj = packmol.pack_box(pdb_filenames, [6, 11, 23])

    xyz = trj.openmm_positions(0)
    top = trj.top.to_openmm()
    top.setUnitCellDimensions(mm.Vec3(*trj.unitcell_lengths[0])*u.nanometer)

    forcefield = app.ForceField(ffxml)

    temperature = 300 * u.kelvin
    friction = 0.1 / u.picosecond
    timestep = 0.1 * u.femtosecond

    system = forcefield.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=None)

    integrator = mm.LangevinIntegrator(temperature, friction, timestep)

    simulation = app.Simulation(top, system, integrator)
    simulation.context.setPositions(xyz)

    simulation.step(25)
