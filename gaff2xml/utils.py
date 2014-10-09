import os
import tempfile
import logging
from pkg_resources import resource_filename
import contextlib
import shutil
import mdtraj as md
from mdtraj.utils import enter_temp_directory
from mdtraj.utils.delay_import import import_

try:
    from subprocess import getoutput  # If python 3
except ImportError:
    from commands import getoutput  # If python 2

import simtk.openmm
from simtk.openmm import app
import simtk.unit as units
from distutils.spawn import find_executable

from gaff2xml import amber_parser, system_checker

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")


def find_gaff_dat():
    AMBERHOME = None
    
    try:
        AMBERHOME = os.environ['AMBERHOME']
    except KeyError:
        pass
    
    if AMBERHOME is None:
        full_path = find_executable("parmchk2")
        try:
            AMBERHOME = os.path.split(full_path)[0]
            AMBERHOME = os.path.join(AMBERHOME, "../")
        except:
            raise(ValueError("Cannot find AMBER GAFF"))

    if AMBERHOME is None:
        raise(ValueError("Cannot find AMBER GAFF"))

    return os.path.join(AMBERHOME, 'dat', 'leap', 'parm', 'gaff.dat')

GAFF_DAT_FILENAME = find_gaff_dat()


def parse_ligand_filename(filename):
    """Split ligand filename into name and extension.  "./ligand.mol2" -> ("ligand", ".mol2")"""
    name, ext = os.path.splitext(os.path.split(filename)[1])
    return name, ext


def run_antechamber(molecule_name, input_filename, charge_method="bcc", net_charge=None):
    """Run AmberTools antechamber and parmchk2 to create GAFF mol2 and frcmod files.

    Parameters
    ----------
    molecule_name : str
        Name of the molecule to be parameterized, will be used in output filenames.
    ligand_filename : str
        The molecule to be parameterized.  Must be tripos mol2 format.
    charge_method : str, optional
        If not None, the charge method string will be passed to Antechamber.
    net_charge : int, optional
        If not None, net charge of the molecule to be parameterized.
        If None, Antechamber sums up partial charges from the input file.

    Returns
    -------
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    """

    ext = parse_ligand_filename(input_filename)[1]

    filetype = ext[1:]
    if filetype != "mol2":
        raise(ValueError("Must input mol2 filename"))

    gaff_mol2_filename = molecule_name + '.gaff.mol2'
    frcmod_filename = molecule_name + '.frcmod'

    cmd = "antechamber -i %s -fi mol2 -o %s -fo mol2 -s 2" % (input_filename, gaff_mol2_filename)
    if charge_method is not None:
        cmd += ' -c %s' % charge_method

    if net_charge is not None:
        cmd += ' -nc %d' % net_charge

    logger.debug(cmd)

    output = getoutput(cmd)
    logger.debug(output)

    cmd = "parmchk2 -i %s -f mol2 -o %s" % (gaff_mol2_filename, frcmod_filename)
    logger.debug(cmd)

    output = getoutput(cmd)
    logger.debug(output)

    return gaff_mol2_filename, frcmod_filename


def convert_molecule(in_filename, out_filename):
    """Use openbabel to convert filenames.  May not work for all file formats!"""

    molecule_name, ext_in = parse_ligand_filename(in_filename)
    molecule_name, ext_out = parse_ligand_filename(out_filename)

    cmd = "obabel -i %s %s -o %s -O %s" % (ext_in[1:], in_filename, ext_out[1:], out_filename)
    print(cmd)
    output = getoutput(cmd)
    logger.debug(output)


def run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename):
    """Run AmberTools tleap to create simulation files for AMBER

    Parameters
    ----------
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk

    Returns
    -------
    prmtop_filename : str
        Amber prmtop file produced by tleap
    inpcrd_filename : str
        Amber inpcrd file produced by tleap
    """

    prmtop_filename = "%s.prmtop" % molecule_name
    inpcrd_filename = "%s.inpcrd" % molecule_name

    tleap_input = """
source leaprc.ff99SB
source leaprc.gaff
LIG = loadmol2 %s
check LIG
loadamberparams %s
saveamberparm LIG %s %s
quit

""" % (gaff_mol2_filename, frcmod_filename, prmtop_filename, inpcrd_filename)

    file_handle = tempfile.NamedTemporaryFile('w')  # FYI Py3K defaults to 'wb' mode, which won't work here.
    file_handle.writelines(tleap_input)
    file_handle.flush()

    cmd = "tleap -f %s " % file_handle.name
    logger.debug(cmd)

    output = getoutput(cmd)
    logger.debug(output)

    file_handle.close()

    return prmtop_filename, inpcrd_filename


def molecule_to_mol2(molecule, tripos_mol2_filename=None):
    """Convert OE molecule to tripos mol2 file.

    Parameters
    ----------
    molecule : openeye.oechem.OEGraphMol
        The molecule to be converted.

    Returns
    -------
    tripos_mol2_filename : str
        Filename of output tripos mol2 file
    
    """
    
    try:
        import openeye.oechem
    except ImportError:
        raise(ImportError("Must install OpenEye tools to process OpenEye MOL2 files."))
    
    # Get molecule name.
    molecule_name = molecule.GetTitle()
    logger.debug(molecule_name)

    # Write molecule as Tripos mol2.
    if tripos_mol2_filename is None:
        tripos_mol2_filename = molecule_name + '.tripos.mol2'

    ofs = openeye.oechem.oemolostream(tripos_mol2_filename)
    ofs.SetFormat(openeye.oechem.OEFormat_MOL2H)
    openeye.oechem.OEWriteMolecule(ofs, molecule)
    ofs.close()

    # Replace <0> substructure names with valid text.
    infile = open(tripos_mol2_filename, 'r')
    lines = infile.readlines()
    infile.close()
    newlines = [line.replace('<0>', 'MOL') for line in lines]
    outfile = open(tripos_mol2_filename, 'w')
    outfile.writelines(newlines)
    outfile.close()

    return molecule_name, tripos_mol2_filename


def create_ffxml_file(gaff_mol2_filenames, frcmod_filenames, ffxml_filename=None, override_mol2_residue_name=None):
    """Process multiple gaff mol2 files and frcmod files using the XML conversion and write to an XML file.

    Parameters
    ----------
    gaff_mol2_filenames : list of str
        The names of the gaff mol2 files
    frcmod_filenames : str
        The names of the gaff frcmod files
    ffxml_filename : str, optional, default=None
        Optional name of output ffxml file to generate.  If None, no file 
        will be generated.
    override_mol2_residue_name : str, default=None
            If given, use this name to override mol2 residue names.        
    
    Returns
    -------
    ffxml_stringio : str
        StringIO representation of ffxml file containing residue entries for each molecule.

    """

    # Generate ffxml file.
    parser = amber_parser.AmberParser(override_mol2_residue_name=override_mol2_residue_name)

    filenames = [GAFF_DAT_FILENAME]
    filenames.extend([filename for filename in gaff_mol2_filenames])
    filenames.extend([filename for filename in frcmod_filenames])

    parser.parse_filenames(filenames)
    
    ffxml_stream = parser.generate_xml()

    if ffxml_filename is not None:
        outfile = open(ffxml_filename, 'w')
        outfile.write(ffxml_stream.read())
        outfile.close()
        ffxml_stream.seek(0)

    return ffxml_stream

def create_ffxml_simulation(molecule_name, gaff_mol2_filename, frcmod_filename):
    """Process a gaff mol2 file and frcmod file using the XML conversion, returning an OpenMM simulation.

    Parameters
    ----------
    molecule_name : str
        The name of the molecule
    gaff_mol2_filename : str
        The name of the gaff mol2 file
    frcmod_filename : str
        The name of the gaff frcmod file

    Returns
    -------
    simulation : openmm.app.Simulation
        A functional simulation object for simulating your molecule
    """

    # Generate ffxml file.
    parser = amber_parser.AmberParser()
    parser.parse_filenames([GAFF_DAT_FILENAME, gaff_mol2_filename, frcmod_filename])

    ffxml_filename = molecule_name + '.ffxml'
    create_ffxml_file([gaff_mol2_filename], [frcmod_filename], ffxml_filename)

    traj = md.load(gaff_mol2_filename)  # Read mol2 file.
    positions = traj.openmm_positions(0)  # Extract OpenMM-united positions of first (and only) trajectory frame
    topology = traj.top.to_openmm()

    # Create System object.
    forcefield = app.ForceField(ffxml_filename)
    system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=None)

    # Create integrator.
    timestep = 1.0 * units.femtoseconds
    integrator = simtk.openmm.VerletIntegrator(timestep)

    # Create simulation.
    platform = simtk.openmm.Platform.getPlatformByName("Reference")
    simulation = app.Simulation(topology, system, integrator, platform=platform)
    simulation.context.setPositions(positions)

    return simulation


def create_leap_simulation(molecule_name, gaff_mol2_filename, frcmod_filename):
    """Create an OpenMM simulation using a Gaff mol2 file and frcmod file.


    Parameters
    ----------
    molecule_name : str
        Name of the molecule
    gaff_mol2_filename : str
        Filename of input (GAFF!) mol2 file
    frcmod_filename : str
        Use this frcmod filename

    """

    # Parameterize system with LEaP.
    (prmtop_filename, inpcrd_filename) = run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename)

    # Create System object.
    prmtop = app.AmberPrmtopFile(prmtop_filename)
    topology = prmtop.topology
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=None)

    # Read positions.
    inpcrd = app.AmberInpcrdFile(inpcrd_filename)
    positions = inpcrd.getPositions()

    # Create integrator.
    timestep = 1.0 * units.femtoseconds
    integrator = simtk.openmm.VerletIntegrator(timestep)

    platform = simtk.openmm.Platform.getPlatformByName("Reference")
    simulation = app.Simulation(topology, system, integrator, platform=platform)
    simulation.context.setPositions(positions)

    return simulation


def test_molecule(molecule_name, tripos_mol2_filename, charge_method="bcc"):
    """Create a GAFF molecule via LEAP and ffXML and compare force terms.


    Parameters
    ----------
    molecule_name : str
        Name of the molecule
    tripos_mol2_filename : str
        Filename of input mol2 file
    charge_method : str, default="bcc"
        If None, use charges in existing MOL2.  Otherwise, use a charge
        model when running antechamber.
    """

    # Generate GAFF parameters.
    (gaff_mol2_filename, frcmod_filename) = run_antechamber(molecule_name, tripos_mol2_filename, charge_method=charge_method)

    # Create simulations.
    simulation_ffxml = create_ffxml_simulation(molecule_name, gaff_mol2_filename, frcmod_filename)
    simulation_leap  = create_leap_simulation(molecule_name, gaff_mol2_filename, frcmod_filename)

    # Compare simulations.
    syscheck = system_checker.SystemChecker(simulation_ffxml, simulation_leap)
    syscheck.check_force_parameters()
    
    groups0, groups1 = syscheck.check_energy_groups()
    energy0, energy1 = syscheck.check_energies()


def get_data_filename(relative_path):
    """Get the full path to one of the reference files shipped for testing

    In the source distribution, these files are in ``gaff2xml/chemicals/*/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the gaff2xml folder).

    """

    fn = resource_filename('gaff2xml', relative_path)

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn



def smiles_to_mdtraj_ffxml(smiles_strings, base_molecule_name="lig"):
    """Generate an MDTraj object from a smiles string.
    
    Parameters
    ----------
    smiles_strings : list(str)
        Smiles strings to create molecules for
    base_molecule_name : str, optional, default='lig'
        Base name of molecule to use inside parameter files.
    
    Returns
    -------
    traj : mdtraj.Trajectory
        MDTraj object for molecule
    ffxml : StringIO
        StringIO representation of ffxml file.
    
    Notes
    -----
    ffxml can be directly input to OpenMM e.g. 
    `forcefield = app.ForceField(ffxml)`
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise(ImportError("Must install rdkit to use smiles conversion."))

    gaff_mol2_filenames = []
    frcmod_filenames = []
    trajectories = []
    for k, smiles_string in enumerate(smiles_strings):
        molecule_name = "%s-%d" % (base_molecule_name, k)
        m = Chem.MolFromSmiles(smiles_string)
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m)
        AllChem.UFFOptimizeMolecule(m)

        pdb_filename = tempfile.mktemp(suffix=".pdb")
        Chem.MolToPDBFile(m, pdb_filename)
        
        mol2_filename = tempfile.mktemp(suffix=".mol2")
        
        convert_molecule(pdb_filename, mol2_filename)  # This is necessary because PDB double bonds are not handled by antechamber...
        print(mol2_filename)

        gaff_mol2_filename, frcmod_filename = run_antechamber(molecule_name, mol2_filename)
        traj = md.load(gaff_mol2_filename)
        print(gaff_mol2_filename)
        print(traj)

        for atom in traj.top.atoms:
            atom.residue.name = molecule_name

        gaff_mol2_filenames.append(gaff_mol2_filename)
        frcmod_filenames.append(frcmod_filename)
        trajectories.append(traj)

    ffxml = create_ffxml_file(gaff_mol2_filenames, frcmod_filenames, override_mol2_residue_name=molecule_name)

    return trajectories, ffxml


def tag_description(lambda_function, description):
    """Add a description flag to a lambda function for nose testing."""
    lambda_function.description = description
    return lambda_function
