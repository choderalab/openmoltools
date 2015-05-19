import mdtraj as md
import tempfile
import logging
import os
from distutils.spawn import find_executable
import openmoltools.utils as utils

try:
    from subprocess import getoutput  # If python 3
except ImportError:
    from commands import getoutput  # If python 2

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")


# http://ambermd.org/tutorials/advanced/tutorial15/Tutorial2.xhtml
# Run tLEaP with input file:
# $ tleap -f commands.in

TLEAP_TEMPLATE = """
source leaprc.gaff
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
setbox box centers
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

#loadmol2_section will look something like this:
#BMI = loadmol2 bmi.mol2
#BF4 = loadmol2 bf4.mol2
#ACN = loadmol2 acn.mol2

#loadamberparams_section looks like this:
#loadamberparams frcmod.bf4
#loadamberparams frcmod.bmi
#loadamberparams frcmod.acn


def build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename):
    """Create a prmtop and inpcrd from a collection of mol2 and frcmod files
    as well as a single box PDB.  We have used this for setting up
    simulations of neat liquids or binary mixtures.  

    Parameters
    ----------
    mol2_filenames : list(str)
        Filenames of GAFF flavored mol2 files.  Each must contain exactly
        ONE ligand. Filenames cannot contain spaces (tleap limitation)
    frcmod_filenames : str
        Filename of input GAFF frcmod filenames.
    box_filename : str
        Filename of PDB containing an arbitrary box of the mol2 molecules.
    prmtop_filename : str
        output prmtop filename.  Should have suffix .prmtop
    inpcrd_filename : str
        output inpcrd filename.  Should have suffix .inpcrd

    Returns
    -------
    tleap_commands : str
        The string of commands piped to tleap for building the prmtop 
        and inpcrd files.  This will *already* have been run, but the
        output can be useful for debugging or archival purposes.
        
    Notes
    -----
    This can be easily broken if there are missing, duplicated, or
    inconsistent ligand residue names in your box, mol2, and frcmod files.
    You can use mdtraj to edit the residue names with something like
    this: trj.top.residue(0).name = "L1"
    """
    
    # Check for one residue name per mol2 file and uniqueness between all mol2 files
    all_names = set()
    for filename in mol2_filenames:
        t = md.load(filename)
        names = set([r.name for r in t.top.residues])
        
        if len(names) != 1:
            raise(ValueError("Must have a SINGLE residue name in each mol2 file."))
        
        all_names = all_names.union(list(names))

    if len(all_names) != len(mol2_filenames):
        raise(ValueError("Must have UNIQUE residue names in each mol2 file."))
    
    #Check for spaces in filenames; AMBER can't handle these.
    for name in mol2_filenames:
        assert ' ' not in name, "Error: tleap cannot process mol2 filenames containing spaces."

    all_names = [md.load(filename).top.residue(0).name for filename in mol2_filenames]
    
    mol2_section = "\n".join("%s = loadmol2 %s" % (all_names[k], filename) for k, filename in enumerate(mol2_filenames))
    amberparams_section = "\n".join("loadamberparams %s" % (filename) for k, filename in enumerate(frcmod_filenames))

    tleap_commands = TLEAP_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section, box_filename=box_filename, prmtop_filename=prmtop_filename, inpcrd_filename=inpcrd_filename)
    print(tleap_commands)
    
    file_handle = tempfile.NamedTemporaryFile('w')  # FYI Py3K defaults to 'wb' mode, which won't work here.
    file_handle.writelines(tleap_commands)
    file_handle.flush()

    cmd = "tleap -f %s " % file_handle.name
    logger.debug(cmd)

    output = getoutput(cmd)
    logger.debug(output)
    check_for_errors( output, other_errors = ['Improper number of arguments'], ignore_errors = ['unperturbed charge of the unit', 'ignoring the error'] )

    file_handle.close()

    return tleap_commands

def check_for_errors( outputtext, other_errors = None, ignore_errors = None ):
    """Check AMBER package output for the string 'ERROR' (upper or lowercase) and (optionally) specified other strings and raise an exception if it is found (to avoid silent failures which might be noted to log but otherwise ignored).

    Parameters
    ----------
    outputtext : str
        String listing output text from an (AMBER) command which should be checked for errors.
    other_errors : list(str), default None
        If specified, provide strings for other errors which will be chcked for, such as "improper number of arguments", etc.
    ignore_errors: list(str), default None
        If specified, AMBER output lines containing errors but also containing any of the specified strings will be ignored (because, for example, AMBER issues an "ERROR" for non-integer charges in some cases when only a warning is needed). 

    Notes
    -----
    If error(s) are found, raise a RuntimeError and attept to print the appropriate errors from the processed text."""
    lines = outputtext.split('\n')
    error_lines = []
    for line in lines:
        if 'ERROR' in line.upper():
            error_lines.append( line )
        if not other_errors == None:
            for err in other_errors:
                if err.upper() in line.upper():
                    error_lines.append( line )

    if not ignore_errors == None and len(error_lines)>0:
        new_error_lines = []
        for ign in ignore_errors:
            ignore = False
            for err in error_lines:
                if ign in err:
                    ignore = True
            if not ignore:
                new_error_lines.append( err )
        error_lines = new_error_lines 

    if len(error_lines) > 0:
        print("Unexpected errors encountered running AMBER tool. Offending output:")
        for line in error_lines: print(line)
        raise(RuntimeError("Error encountered running AMBER tool. Exiting."))

    return


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


def run_antechamber(molecule_name, input_filename, charge_method="bcc", net_charge=None, gaff_mol2_filename=None, frcmod_filename=None):
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
    gaff_mol2_filename : str, optional, default=None
        Name of GAFF mol2 filename to output.  If None, uses local directory
        and molecule_name
    frcmod_filename : str, optional, default=None
        Name of GAFF frcmod filename to output.  If None, uses local directory
        and molecule_name

    Returns
    -------
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    """

    ext = utils.parse_ligand_filename(input_filename)[1]

    filetype = ext[1:]
    if filetype != "mol2":
        raise(ValueError("Must input mol2 filename"))


    if gaff_mol2_filename is None:
        gaff_mol2_filename = molecule_name + '.gaff.mol2'
    if frcmod_filename is None:
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
    check_for_errors( output  )

    return gaff_mol2_filename, frcmod_filename


def run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename, prmtop_filename=None, inpcrd_filename=None):
    """Run AmberTools tleap to create simulation files for AMBER

    Parameters
    ----------
    molecule_name : str
        The name of the molecule    
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    prmtop_filename : str, optional, default=None
        Amber prmtop file produced by tleap, defaults to molecule_name
    inpcrd_filename : str, optional, default=None
        Amber inpcrd file produced by tleap, defaults to molecule_name  

    Returns
    -------
    prmtop_filename : str
        Amber prmtop file produced by tleap
    inpcrd_filename : str
        Amber inpcrd file produced by tleap
    """
    if prmtop_filename is None:
        prmtop_filename = "%s.prmtop" % molecule_name
    if inpcrd_filename is None:
        inpcrd_filename = "%s.inpcrd" % molecule_name
    
    assert ' ' not in gaff_mol2_filename, "Error: tleap cannot process mol2 filenames containing spaces."
    assert ' ' not in frcmod_filename, "Error: tleap cannot process filenames containing spaces."
    assert ' ' not in prmtop_filename, "Error: tleap cannot process filenames containing spaces."
    assert ' ' not in inpcrd_filename, "Error: tleap cannot process filenames containing spaces."

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

    check_for_errors( output, other_errors = ['Improper number of arguments'] )

    file_handle.close()

    return prmtop_filename, inpcrd_filename

