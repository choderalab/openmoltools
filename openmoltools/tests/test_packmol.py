import os
import tempfile
import mdtraj as md
from unittest import skipIf
import logging
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


def test_standardize_water():
    """Test utility function standardize_water.

    The water bonds must be recognized even when residue names do not
    match the standard definition in mdtraj.formats.pdb.data.residues.xml.

    """
    water_filepath = utils.get_data_filename("chemicals/water/water.mol2")
    water_traj = md.load(water_filepath)

    # Store in pdb format and lose CONECT records.
    water_pdb_filepath = tempfile.mktemp(suffix='.pdb')
    water_traj.save_pdb(water_pdb_filepath)
    with open(water_pdb_filepath, 'r') as f:
        pdb_lines = f.readlines()
    with open(water_pdb_filepath, 'w') as f:
        for line in pdb_lines:
            if line[:6] != 'CONECT':
                f.write(line)

    # Test pre-condition: MDTraj cannot detect water bonds automatically.
    water_traj = md.load(water_pdb_filepath)
    assert water_traj.topology.n_bonds == 0

    # The function modifies the Trajectory and bonds are now recognized.
    assert packmol.standardize_water(water_traj) is True
    assert water_traj.topology.n_bonds == 2

    # The second time, the Trajectory object is not modified.
    assert packmol.standardize_water(water_traj) is False

    # Remove temporary file.
    os.remove(water_pdb_filepath)


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
    friction = 91.0 / u.picosecond
    timestep = 0.1 * u.femtosecond

    system = forcefield.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=None)

    integrator = mm.LangevinIntegrator(temperature, friction, timestep)

    simulation = app.Simulation(top, system, integrator)
    simulation.context.setPositions(xyz)
    simulation.minimizeEnergy()
    simulation.step(25)

@skipIf(not HAVE_RDKIT, "Skipping testing of packmol conversion because rdkit not found.")
@skipIf(packmol.PACKMOL_PATH is None, "Skipping testing of packmol conversion because packmol not found.")
def test_packmol_simulation_ternary_bydensity():
    smiles_list = ["Cc1ccccc1", "c1ccccc1", "CC"]
    pdb_filenames = []
    ligand_trajectories, ffxml = utils.smiles_to_mdtraj_ffxml(smiles_list)

    for k, ligand_traj in enumerate(ligand_trajectories):
        pdb_filename = tempfile.mktemp(suffix=".pdb")
        ligand_traj.save(pdb_filename)
        pdb_filenames.append(pdb_filename)

    pdb_filenames = pdb_filenames[0:2] + [ligand_traj]  # Test passing BOTH pdb filenames and trajectories as input.
    size = packmol.approximate_volume_by_density( smiles_list, [12, 22, 46] )
    trj = packmol.pack_box(pdb_filenames, [12, 22, 46], box_size = size)

    xyz = trj.openmm_positions(0)
    top = trj.top.to_openmm()
    top.setUnitCellDimensions(mm.Vec3(*trj.unitcell_lengths[0])*u.nanometer)

    forcefield = app.ForceField(ffxml)

    temperature = 300 * u.kelvin
    friction = 91.0 / u.picosecond
    timestep = 0.1 * u.femtosecond

    system = forcefield.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=None)

    integrator = mm.LangevinIntegrator(temperature, friction, timestep)

    simulation = app.Simulation(top, system, integrator)
    simulation.context.setPositions(xyz)
    simulation.minimizeEnergy()

    simulation.step(25)
