import random
import itertools
import string
import os
import tempfile
import logging
import subprocess
from pkg_resources import resource_filename
import contextlib
import shutil
import mdtraj as md
from mdtraj.utils import enter_temp_directory
from mdtraj.utils.delay_import import import_

import simtk.openmm
from simtk.openmm import app
import simtk.unit as units
from distutils.spawn import find_executable

from openmoltools import amber_parser, system_checker

logger = logging.getLogger(__name__)


def getoutput(cmd):
    """Compatibility function to substitute deprecated commands.getoutput in Python2.7"""
    try:
        out = subprocess.getoutput(cmd)
    except AttributeError:
        out = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT,
                               stdout=subprocess.PIPE).communicate()[0]
    try:
        return str(out.decode())
    except:
        return str(out)


def find_gaff_dat():
    print("Warning: find_gaff_dat has been moved to openmoltools.amber.")
    amber = import_("openmoltools.amber")

    return amber.find_gaff_dat()


def parse_ligand_filename(filename):
    """Split ligand filename into name and extension.  "./ligand.mol2" -> ("ligand", ".mol2")"""
    name, ext = os.path.splitext(os.path.split(filename)[1])
    return name, ext


def run_antechamber(*args, **kwargs):
    print("Warning: run_antechamber has been moved to openmoltools.amber.")
    amber = import_("openmoltools.amber")
    return amber.run_antechamber(*args, **kwargs)

def run_tleap(*args, **kwargs):
    print("Warning: run_tleap has been moved to openmoltools.amber.")
    amber = import_("openmoltools.amber")
    return amber.run_tleap(*args, **kwargs)

def convert_via_acpype( molecule_name, in_prmtop, in_crd, out_top = None, out_gro = None, debug = False, is_sorted = False ):
    """Use acpype.py (Sousa Da Silva et al., BMC Research Notes 5:367 (2012)) to convert AMBER prmtop and crd files to GROMACS format using amb2gmx mode. Writes to GROMACS 4.5 (and later) format, rather than the format for earlier GROMACS versions.


    Parameters
    ----------
    molecule_name : str
        String specifying name of molecule
    in_prmtop : str
        String specifying path to AMBER-format parameter/topology (parmtop) file
    in_crd : str
        String specifying path to AMBER-format coordinate file
    out_top : str, optional, default = None
        String specifying path to GROMACS-format topology file which will be written out. If none is provided, created based on molecule_name.
    out_gro : str, optional, default = None
        String specifying path to GROMACS-format coordinate (.gro) file which will be written out. If none is provided, created based on molecule_name.
    debug : bool, optional, default = False
        Print debug info? If not specified, do not.
    is_sorted : bool, optional, default = False
        Sort resulting topology file

    Returns
    -------
    out_top : str
        GROMACS topology file produced by acpype
    out_gro : str
        GROMACS coordinate file produced by acpype

    Notes
    -----
    Deprecated. Please use ParmEd (especially amber_to_gromacs) instead.
    """

    print("WARNING: Use of acpype for conversion is deprecated. ParmEd is preferred; please use amber_to_gromacs instead.")
    acpype = import_('openmoltools.acpype')

    #Create output file names if needed
    if out_top is None:
        out_top = "%s.top" % molecule_name
    if out_gro is None:
        out_gro = "%s.gro" % molecule_name

    #Create temporary output dir for acpype output
    outdir = tempfile.mkdtemp()
    #Define basename for output
    basename = os.path.join( outdir, 'output')


    #Set up acpype
    system = acpype.MolTopol( acFileXyz = in_crd, acFileTop = in_prmtop, basename = basename, is_sorted = is_sorted, gmx45 = True, disam = True )

    #Print debug info if desired
    if debug:
        print(system.printDebug('prmtop and inpcrd files parsed'))

    #Write results
    system.writeGromacsTopolFiles( amb2gmx = True )

    #Acpype names various things in the topology and coordinate file after the base name of the file used as input. Replace these names with an at-least-legible string while writing to desired output
    top_in = open(basename+"_GMX.top", 'r')
    top_out = open( out_top, 'w')
    for line in top_in.readlines():
        top_out.write( line.replace( basename, molecule_name) )
    top_in.close()
    top_out.close()
    gro_in = open(basename+"_GMX.gro", 'r')
    gro_out = open( out_gro, 'w')
    for line in gro_in.readlines():
        gro_out.write( line.replace( basename, molecule_name) )
    gro_in.close()
    gro_out.close()

    #Check if files exist and are not empty; return True if so
    if os.stat( out_top).st_size == 0 or os.stat( out_gro ) == 0:
        raise(ValueError("ACPYPE conversion failed."))

    return out_top, out_gro


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

    amber = import_("openmoltools.amber")
    GAFF_DAT_FILENAME = amber.find_gaff_dat()
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
    amber = import_("openmoltools.amber")
    GAFF_DAT_FILENAME = amber.find_gaff_dat()
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
    amber = import_("openmoltools.amber")
    (gaff_mol2_filename, frcmod_filename) = amber.run_antechamber(molecule_name, tripos_mol2_filename, charge_method=charge_method)

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

    In the source distribution, these files are in ``openmoltools/chemicals/*/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the openmoltools folder).

    """

    fn = resource_filename('openmoltools', relative_path)

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn


@contextlib.contextmanager
def temporary_cd(dir_path):
    """Context to temporary change the working directory."""
    prev_dir = os.getcwd()
    os.chdir(os.path.abspath(dir_path))
    try:
        yield
    finally:
        os.chdir(prev_dir)


@contextlib.contextmanager
def temporary_directory():
    """Context for safe creation of temporary directories."""
    tmp_dir = tempfile.mkdtemp()
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir)


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

        mdl_filename = tempfile.mktemp(suffix=".mdl")
        Chem.MolToMolFile(m, mdl_filename)
        amber = import_("openmoltools.amber")
        gaff_mol2_filename, frcmod_filename = amber.run_antechamber(molecule_name, mdl_filename, input_format='mdl')
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


def molecule_to_mol2(*args, **kwargs):
    print("Warning: molecule_to_mol2 has been moved to openmoltools.openeye.")
    import openmoltools.openeye
    return openmoltools.openeye.molecule_to_mol2(*args, **kwargs)

def get_unique_names(n_molecules):
    """Generate unique random residue names for use in mixture mol2 / pdb files.

    Parameters
    ----------
    n_molecules : int
        Number of unique names to generate

    Notes
    -----
    Names will start with Z to avoid conflicts with common macromolecule
    residue names.  This may be improved in the future.

    THIS FUNCTION will enter an INFINITE LOOP if you request many
    (hundreds) of unique residue names, as it becomes harder or impossible
    to generate many unique names.
    """
    for i in itertools.count():
        names = ["Z" + ''.join(random.choice(string.ascii_uppercase) for _ in range(2)) for i in range(n_molecules)]
        if len(set(names)) == n_molecules:
            return names


def randomize_mol2_residue_names(mol2_filenames):
    """Find unique residue names for a list of MOL2 files.  Then
    re-write the MOL2 files using ParmEd with the unique identifiers.
    """
    import parmed
    names = get_unique_names(len(mol2_filenames))

    for k, filename in enumerate(mol2_filenames):
        struct = parmed.load_file(filename)
        struct.name = names[k]
        mol2file = parmed.formats.Mol2File
        mol2file.write(struct, filename)

def get_checkmol_descriptors( molecule_filename, executable_name = 'checkmol' ):
    """For a specified molecule file, return a list of functional groups as assigned by checkmol for the molecule(s) present. The first entry in the list will correspond to the groups in the first molecule, the second gives groups in the second (if present) and so on. Raises an exception if checkmol is not found.

    Parameters
    ----------
    molecule_filename : str
        Specifies name of file to read
    executable_name : str, default = 'checkmol'
        Specify name (or full path) of execuable for checkmol

    Returns
    -------
    descriptors : list (of lists of strings)
        Checkmol functional group assignments for each molecule(s) in the input file, where descriptors[0] gives the descriptors for the first molecule, etc.

    Notes
    -----
    This should properly handle single-molecule and multiple-molecule files; however, multiple-conformer files may result in each conformer appearing (rather than each molecule) appearing in the list of descriptors, which may or may not be the expected behavior.
    """

    oechem = import_("openeye.oechem")

    status = find_executable( executable_name )
    if status==None:
        raise(ValueError("Cannot find checkmol; cannot assign checkmol descriptors without it."))


    #Open input file
    ifs = oechem.oemolistream( molecule_filename )
    #Input molecule
    mol = oechem.OEGraphMol( )

    #Set up temporary file for molecule output
    fname = tempfile.mktemp( suffix = '.sdf' )

    #Storage for descriptors
    descriptors = []

    #Read/write/run checkmol
    while oechem.OEReadMolecule( ifs, mol ):
        #Dump molecule out
        ofs = oechem.oemolostream( fname )
        oechem.OEWriteMolecule( ofs, mol )
        ofs.close()
        #Run checkmol
        groups = getoutput('%s %s' % (executable_name, fname) )
        #Split to separate groups
        groups = groups.split('\n')
        #Store results
        descriptors.append( groups )

    #Raise an exception if the whole list is empty
    fnd = False
    for elem in descriptors:
        if len(elem)>0:
            fnd = True
    if not fnd:
        raise(ValueError("checkmol only produced empty descriptors for your molecule. Something is wrong; please check your input file and checkmol installation."))

    #Delete temporary file
    os.remove( fname )

    return descriptors

def amber_to_gromacs( molecule_name, in_prmtop, in_crd, out_top = None, out_gro = None, precision = None):
    """Use ParmEd to convert AMBER prmtop and crd files to GROMACS format.

    Requires
    --------
    Currently requires ParmEd v2.0 beta1 or later.


    Parameters
    ----------
    molecule_name : str
        String specifying name of molecule
    in_prmtop : str
        String specifying path to AMBER-format parameter/topology (parmtop) file
    in_crd : str
        String specifying path to AMBER-format coordinate file
    out_top : str, optional, default = None
        String specifying path to GROMACS-format topology file which will be written out. If none is provided, created based on molecule_name.
    out_gro : str, optional, default = None
        String specifying path to GROMACS-format coordinate (.gro) file which will be written out. If none is provided, created based on molecule_name.
    precision : int, optional, default = None
        If not none, set the precision of the coordinates in the written .gro file to the specified number of decimal places.

    Returns
    -------
    out_top : str
        GROMACS topology file produced by ParmEd
    out_gro : str
        GROMACS coordinate file produced by ParmEd

    Notes
    -----
        molecule_name is not currently used except to generate output file names if gro/top file names are not provided. It is an argument partly for API consistency.
    """
    #Create output file names if needed
    if out_top is None:
        out_top = "%s.top" % molecule_name
    if out_gro is None:
        out_gro = "%s.gro" % molecule_name

    #Check precision
    if precision is not None:
        assert isinstance(precision, int), "Precision %s is not an integer." % precision

    #Import ParmEd
    import parmed

    #Read AMBER to ParmEd object
    structure = parmed.amber.AmberParm( in_prmtop, in_crd )
    #Make GROMACS topology
    gromacs_topology = parmed.gromacs.GromacsTopologyFile.from_structure( structure )
    #Write
    parmed.gromacs.GromacsTopologyFile.write( gromacs_topology, out_top )
    if precision == None:
        parmed.gromacs.GromacsGroFile.write( gromacs_topology, out_gro )
    else:
        parmed.gromacs.GromacsGroFile.write( gromacs_topology, out_gro, precision = precision )

    return out_top, out_gro
