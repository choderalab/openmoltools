import numpy as np
import mdtraj as md
from unittest import skipIf
import logging
from mdtraj.testing import eq
from openmoltools import utils
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

CHECKMOL_PATH = find_executable("checkmol")
SKIP_CHECKMOL = (CHECKMOL_PATH is None)

logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")

def test_enter_temp_directory():
    with utils.enter_temp_directory():
        pass

def test_parse_ligand_filename():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    name, ext = utils.parse_ligand_filename(input_filename)
    
    eq(name, "sustiva")
    eq(ext, ".mol2")

def test_run_test_molecule():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        utils.test_molecule(molecule_name, input_filename)

def test_acpype_conversion():
    molecule_name = 'sustiva'
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    with utils.enter_temp_directory(): # Prevents creating tons of GAFF files everywhere.
        gaff_mol2_filename, frcmod_filename = utils.run_antechamber(molecule_name, input_filename, charge_method=None)
        prmtop, inpcrd = utils.run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename)
        out_top, out_gro = utils.convert_via_acpype( molecule_name, prmtop, inpcrd ) 


@skipIf(SKIP_CHECKMOL, "Skipping testing of checkmol descriptors since checkmol is not found (under that name)." )
def test_checkmol_descriptors():
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    utils.get_checkmol_descriptors( input_filename )


@skipIf(SKIP_SMILES, "Skipping testing of smiles conversion because openbabel or rdkit not found.")
def test_smiles_conversion():
    pdb_filename = utils.get_data_filename("chemicals/proteins/1vii.pdb")
    smiles = 'Cc1ccccc1'  # Also known as toluene.

    temperature = 300 * u.kelvin
    friction = 0.3 / u.picosecond
    timestep = 0.01 * u.femtosecond

    protein_traj = md.load(pdb_filename)
    protein_traj.center_coordinates()

    protein_top = protein_traj.top.to_openmm()
    protein_xyz = protein_traj.openmm_positions(0)

    ligand_trajectories, ffxml = utils.smiles_to_mdtraj_ffxml([smiles])
    ligand_traj = ligand_trajectories[0]
    ligand_traj.center_coordinates()
    
    eq(ligand_traj.n_atoms, 15)
    eq(ligand_traj.n_frames, 1)

    #Move the pre-centered ligand sufficiently far away from the protein to avoid a clash.  
    min_atom_pair_distance = ((ligand_traj.xyz[0] ** 2.).sum(1) ** 0.5).max() + ((protein_traj.xyz[0] ** 2.).sum(1) ** 0.5).max() + 0.3
    ligand_traj.xyz += np.array([1.0, 0.0, 0.0]) * min_atom_pair_distance

    ligand_xyz = ligand_traj.openmm_positions(0)
    ligand_top = ligand_traj.top.to_openmm()

    forcefield = app.ForceField("amber10.xml", ffxml, "tip3p.xml")

    model = app.modeller.Modeller(protein_top, protein_xyz)
    model.add(ligand_top, ligand_xyz)
    model.addSolvent(forcefield, padding=0.4 * u.nanometer)

    system = forcefield.createSystem(model.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=app.HAngles)

    integrator = mm.LangevinIntegrator(temperature, friction, timestep)

    simulation = app.Simulation(model.topology, system, integrator)
    simulation.context.setPositions(model.positions)
    print("running")
    simulation.step(1)
