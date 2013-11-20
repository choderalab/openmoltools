import os
import os.path
import tempfile
import sys
import logging
from pkg_resources import resource_filename

import openeye.oechem

import simtk.openmm
from simtk.openmm import app
import simtk.unit as units

from gaff2xml import amber_parser, gafftools, system_checker

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")

AMBERHOME = os.environ['AMBERHOME']
GAFF_DAT_FILENAME = os.path.join(AMBERHOME, 'dat', 'leap', 'parm', 'gaff.dat')

def parse_ligand_filename(filename):
    """Split ligand filename into name and extension.  "./ligand.mol2" -> ("ligand", ".mol2")"""
    name, ext = os.path.splitext(os.path.split(filename)[1])
    return name, ext

def run_antechamber(molecule_name, input_filename, charge_method=None):
    """Run AmberTools antechamber and parmchk to create GAFF mol2 and frcmod files.

    Parameters
    ----------
    molecule_name : str
        Name of the molecule to be parameterized, will be used in output filenames.
    ligand_filename : str
        The molecule to be parameterized.  Must be tripos mol2 format.
    charge_method : str, optional
        If not None, the charge method string will be passed to Antechamber.

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

    logger.debug(cmd)
    output = os.system(cmd)
    logger.debug(output)
    cmd = "parmchk -i %s -f mol2 -o %s" % (gaff_mol2_filename, frcmod_filename)
    logger.debug(cmd)
    output = os.system(cmd)
    logger.debug(output)

    return gaff_mol2_filename, frcmod_filename


def convert_molecule(in_filename, out_filename):
    """Use openbabel to convert filenames.  May not work for all file formats!"""
    
    molecule_name, ext_in = parse_ligand_filename(in_filename)
    molecule_name, ext_out = parse_ligand_filename(out_filename)
    
    cmd = "obabel -i %s %s -o %s > %s" % (ext_in, in_filename, ext_out, out_filename)
    os.system(cmd)


def run_tleap(ligand_name, gaff_mol2_filename, frcmod_filename):
    """Run AmberTools tleap to create simulation files for AMBER

    Parameters
    ----------
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk

    Returns
    -------
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    """
    
    prmtop_filename = "%s.prmtop" % ligand_name
    inpcrd_filename = "%s.inpcrd" % ligand_name

    tleap_input = """
source leaprc.ff99SB
source leaprc.gaff
LIG = loadmol2 %s 
check LIG
loadamberparams %s
saveamberparm LIG %s %s
quit

""" % (gaff_mol2_filename, frcmod_filename, prmtop_filename, inpcrd_filename)

    file_handle = tempfile.NamedTemporaryFile()
    file_handle.writelines(tleap_input)
    file_handle.flush()

    cmd = "tleap -f %s " % file_handle.name
    logger.debug(cmd)
    output = os.system(cmd)
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
    newlines = [ line.replace('<0>', 'MOL') for line in lines ]
    outfile = open(tripos_mol2_filename, 'w')
    outfile.writelines(newlines)
    outfile.close()
    
    return tripos_mol2_filename


def create_ffxml_simulation(molecule_name, gaff_mol2_filename, frcmod_filename):
   
    # Generate ffxml file.
    parser = amber_parser.AmberParser()
    parser.parse_filenames([GAFF_DAT_FILENAME, gaff_mol2_filename, frcmod_filename])
    
    ffxml_stream = parser.generate_xml()
    ffxml_filename = molecule_name + '.ffxml'
    outfile = open(ffxml_filename, 'w')
    outfile.write(ffxml_stream.read())
    outfile.close()

    mol2 = gafftools.Mol2Parser(gaff_mol2_filename)  # Read mol2 file.
    (topology, positions) = mol2.to_openmm()

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

def test_molecule(molecule, charge_method=None):
    
    molecule_name = molecule.GetTitle()
    logger.debug(molecule_name)
    
    # Create temporary directory.
    tmp_path = 'tmp' # DEBUG: Create an actual temporary directory.
    if not os.path.exists(tmp_path): os.makedirs(tmp_path)
    logger.debug('temporary directory created in %s' % tmp_path)
    cwd = os.getcwd()
    os.chdir(tmp_path)
    
    tripos_mol2_filename = molecule_to_mol2(molecule)
    
    
    # Generate GAFF parameters.
    (gaff_mol2_filename, frcmod_filename) = run_antechamber(molecule_name, tripos_mol2_filename, charge_method=charge_method)

    # Create simulations.
    simulation_ffxml = create_ffxml_simulation(molecule_name, gaff_mol2_filename, frcmod_filename)
    simulation_leap  = create_leap_simulation(molecule_name, gaff_mol2_filename, frcmod_filename)
    
    # Compare simulations.
    syscheck = system_checker.SystemChecker(simulation_ffxml, simulation_leap)
    syscheck.check_force_parameters()

    # TODO: Clean up temporary directory.

    # Restore current working directory.
    os.chdir(cwd)


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
