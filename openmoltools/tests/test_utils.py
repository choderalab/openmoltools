import os
import numpy as np
import mdtraj as md
from unittest import skipIf
import logging
from mdtraj.testing import eq
from openmoltools import utils
from openmoltools import amber
import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
import simtk.openmm.openmm as mmmm
from distutils.spawn import find_executable
import parmed


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


def test_temp_dir_context():
    """Test the context temporary_directory()."""
    with utils.temporary_directory() as tmp_dir:
        assert os.path.isdir(tmp_dir)
    assert not os.path.exists(tmp_dir)


def test_temp_cd_context():
    """Test the context temporary_cd()."""
    with utils.temporary_directory() as tmp_dir:
        with utils.temporary_cd(tmp_dir):
            assert os.getcwd() == os.path.realpath(tmp_dir)
        assert os.getcwd() != os.path.realpath(tmp_dir)


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
        gaff_mol2_filename, frcmod_filename = amber.run_antechamber(molecule_name, input_filename, charge_method=None)
        prmtop, inpcrd = amber.run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename)
        out_top, out_gro = utils.convert_via_acpype( molecule_name, prmtop, inpcrd ) 

def test_parmed_conversion():
    molecule_name = 'sustiva'
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    with utils.enter_temp_directory(): # Prevents creating tons of GAFF files everywhere.
        #Make sure conversion runs
        gaff_mol2_filename, frcmod_filename = amber.run_antechamber(molecule_name, input_filename, charge_method=None)
        prmtop, inpcrd = amber.run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename)
        out_top, out_gro = utils.amber_to_gromacs( molecule_name, prmtop, inpcrd, precision = 8 ) 

        #Test energies before and after conversion
        #Set up amber system
        a = parmed.amber.AmberParm( prmtop, inpcrd )
        ambersys = a.createSystem()
        ambercon = mmmm.Context( ambersys, mm.VerletIntegrator(0.001))
        ambercon.setPositions( a.positions )
        #Set up GROMACS system
        g = parmed.load_file( out_top )
        gro = parmed.gromacs.GromacsGroFile.parse( out_gro ) 
        g.box = gro.box
        g.positions = gro.positions
        gromacssys = g.createSystem()
        gromacscon = mmmm.Context( gromacssys, mm.VerletIntegrator(0.001))
        gromacscon.setPositions( g.positions ) 

        #Check energies
        a_energies = parmed.openmm.utils.energy_decomposition( a, ambercon )    
        g_energies = parmed.openmm.utils.energy_decomposition( g, gromacscon )
        #Check components
        tolerance = 1e-5
        ok = True
        for key in a_energies.keys():
            diff = np.abs(a_energies[key] - g_energies[key] )
            if diff/np.abs(a_energies[key]) > tolerance:
                ok = False
                print("In testing AMBER to GROMACS conversion, %s energy differs by %.5g, which is more than a fraction %.2g of the total, so conversion appears not to be working properly." % ( key, diff, tolerance) )
        if not ok:
            raise(ValueError("AMBER to GROMACS conversion yields energies which are too different.")) 
    

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
