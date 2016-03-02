import mdtraj as md
import tempfile
import logging
import os, os.path
import shutil
from distutils.spawn import find_executable
from mdtraj.utils.delay_import import import_
import mdtraj.utils

from openmoltools.utils import getoutput

logger = logging.getLogger(__name__)


# http://ambermd.org/tutorials/advanced/tutorial15/Tutorial2.xhtml
# Run tLEaP with input file:
# $ tleap -f commands.in

TLEAP_TEMPLATE = """
source leaprc.gaff
source oldff/leaprc.ff99SB
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


def build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename, water_model = 'TIP3P'):
    """Create a prmtop and inpcrd from a collection of mol2 and frcmod files
    as well as a single box PDB.  We have used this for setting up
    simulations of neat liquids or binary mixtures.

    Parameters
    ----------
    mol2_filenames : list(str)
        Filenames of GAFF flavored mol2 files.  Each must contain exactly
        ONE ligand.
    frcmod_filenames : str
        Filename of input GAFF frcmod filenames.
    box_filename : str
        Filename of PDB containing an arbitrary box of the mol2 molecules.
    prmtop_filename : str
        output prmtop filename.  Should have suffix .prmtop
    inpcrd_filename : str
        output inpcrd filename.  Should have suffix .inpcrd
    water_model : str, optional. Default: "TIP3P"
        String specifying water model to be used IF water is present as a component of the mixture. Valid options are currently "TIP3P", "SPC", or None. If None is specified, flexible GAFF-water will be used as for any other solute (old behavior).

    Returns
    -------
    tleap_commands : str
        The string of commands piped to tleap for building the prmtop
        and inpcrd files.  This will *already* have been run, but the
        output can be useful for debugging or archival purposes. However,
        this will reflect temporary file names for both input and output
        file as these are used to avoid tleap filename restrictions.

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
    if len(mol2_filenames) != len(frcmod_filenames):
        raise(ValueError("Must provide an equal number of frcmod and mol2 file names."))

    #Get number of files
    nfiles = len(mol2_filenames)

    #Check validity of water model options
    valid_water = ['TIP3P', 'SPC', None]
    if not water_model in valid_water:
        raise(ValueError("Must provide a valid water model."))

    #If we are requesting a different water model, check if there is water present
    if not water_model==None:
        parmed = import_("parmed")
        solventIsWater = []
        waterPresent = False
        for i in range(nfiles):
            mol = parmed.load_file( mol2_filenames[i] )
            #Check if it is water by checking GAFF atom names
            types = [ atom.type for atom in mol.atoms ]
            if 'oh' in types and types.count('ho')==2 and len(types)==3:
                solventIsWater.append(True)
                waterPresent = True
            else:
                solventIsWater.append(False)

        #In this case, if we have any water, we will now work on fewer .mol2 and .frcmod files and instead use the force field files for those. So, reduce nfiles and remove the files we don't need from the .mol2 and .frcmod filename lists
        #After doing so, go on to interpret the specified water model and compose the water model string needed for tleap
        if waterPresent:
            new_mol2_filenames = []
            new_frcmod_filenames = []
            water_mol2_filenames = []
            for i in range( nfiles ):
                if not solventIsWater[i]:
                    new_mol2_filenames.append( mol2_filenames[i] )
                    new_frcmod_filenames.append( frcmod_filenames[i] )
                else:
                    water_mol2_filenames.append( mol2_filenames[i] )
            nfiles = len(new_mol2_filenames)
            mol2_filenames = new_mol2_filenames
            frcmod_filenames = new_frcmod_filenames

            #Now interpret the specified water model and translate into AMBER nomenclature
            if water_model=='TIP3P':
                water_model = 'TP3'
            elif water_model =='SPC':
                water_model = 'SPC'
            else:
                raise(ValueError("Cannot translate specified water model into one of the available models."))


            #Compose string for loading specified water molecule
            water_string = '\n'
            water_names = [md.load(filename).top.residue(0).name for filename in water_mol2_filenames]
            for name in water_names:
                water_string += '%s = %s\n' % (name, water_model )
                #Also if not TIP3P, update to source correct frcmod file
                if water_model == 'SPC':
                    water_string += 'loadamberparams frcmod.spce\n'
                elif water_model =='TP3':
                    continue
                else:
                    raise(ValueError("Cannot identify water frcmod file to be loaded."))

            #Rename water atoms in box file to match what is expected by AMBER
            packmol = import_("openmoltools.packmol")
            packmol.rename_water_atoms(box_filename)
    else:
        waterPresent = False

    #Make temporary, hardcoded filenames for mol2 and frcmod input to avoid tleap filename restrictions
    tmp_mol2_filenames = [ 'in%d.mol2' % n for n in range(nfiles) ]
    tmp_frcmod_filenames = [ 'in%d.frcmod' % n for n in range(nfiles) ]

    #Make temporary, hardcoded filenames for output files to avoid tleap filename restrictions
    tmp_prmtop_filename = 'out.prmtop'
    tmp_inpcrd_filename = 'out.inpcrd'
    tmp_box_filename = 'tbox.pdb'

    #Build absolute paths of input files so we can use context and temporary directory
    infiles = mol2_filenames + frcmod_filenames + [box_filename]
    infiles = [ os.path.abspath(filenm) for filenm in infiles ]

    #Build absolute paths of output files so we can copy them back
    prmtop_filename = os.path.abspath( prmtop_filename )
    inpcrd_filename = os.path.abspath( inpcrd_filename )

    #Use temporary directory and do the setup
    with mdtraj.utils.enter_temp_directory():

        #Copy input files to temporary file names in target directory
        for (infile, outfile) in zip( infiles, tmp_mol2_filenames+tmp_frcmod_filenames+[tmp_box_filename] ):
            shutil.copy( infile, outfile)
            logger.debug('Copying input file %s to %s...\n' % (infile, outfile))


        all_names = [md.load(filename).top.residue(0).name for filename in tmp_mol2_filenames]

        mol2_section = "\n".join("%s = loadmol2 %s" % (all_names[k], filename) for k, filename in enumerate(tmp_mol2_filenames))
        #If non-GAFF water is present, load desired parameters for that water as well.
        if waterPresent:
            mol2_section += water_string
        amberparams_section = "\n".join("loadamberparams %s" % (filename) for k, filename in enumerate(tmp_frcmod_filenames))

        tleap_commands = TLEAP_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section, box_filename=tmp_box_filename, prmtop_filename=tmp_prmtop_filename, inpcrd_filename=tmp_inpcrd_filename)
        print(tleap_commands)

        file_handle = open('tleap_commands', 'w')
        file_handle.writelines(tleap_commands)
        file_handle.close()

        logger.debug('Running tleap in temporary directory.')
        cmd = "tleap -f %s " % file_handle.name
        logger.debug(cmd)

        output = getoutput(cmd)
        logger.debug(output)
        check_for_errors( output, other_errors = ['Improper number of arguments'], ignore_errors = ['unperturbed charge of the unit', 'ignoring the error'] )

        #Copy stuff back to right filenames
        for (tfile, finalfile) in zip( [tmp_prmtop_filename, tmp_inpcrd_filename], [prmtop_filename, inpcrd_filename] ):
            shutil.copy( tfile, finalfile)

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


def run_antechamber(molecule_name, input_filename, charge_method="bcc", net_charge=None, gaff_mol2_filename=None, frcmod_filename=None,
    input_format='mol2', resname=False, log_debug_output=False):
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
    input_format : str, optional, default='mol2'
        Format specifier for input file to pass to antechamber.
    resname : bool, optional, default=False
        Set the residue name used within output files to molecule_name
    log_debug_output : bool, optional, default=False
        If true, will send output of tleap to logger.

    Returns
    -------
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    """
    utils = import_("openmoltools.utils")
    ext = utils.parse_ligand_filename(input_filename)[1]

    if gaff_mol2_filename is None:
        gaff_mol2_filename = molecule_name + '.gaff.mol2'
    if frcmod_filename is None:
        frcmod_filename = molecule_name + '.frcmod'

    #Build absolute paths for input and output files
    gaff_mol2_filename = os.path.abspath( gaff_mol2_filename )
    frcmod_filename = os.path.abspath( frcmod_filename )
    input_filename = os.path.abspath( input_filename )

    def read_file_contents(filename):
        infile = open(filename, 'r')
        contents = infile.read()
        infile.close()
        return contents

    #Use temporary directory context to do this to avoid issues with spaces in filenames, etc.
    with mdtraj.utils.enter_temp_directory():
        local_input_filename = 'in.' + input_format
        shutil.copy( input_filename, local_input_filename )

        # Run antechamber.
        cmd = "antechamber -i %(local_input_filename)s -fi %(input_format)s -o out.mol2 -fo mol2 -s 2" % vars()
        if charge_method is not None:
            cmd += ' -c %s' % charge_method
        if net_charge is not None:
            cmd += ' -nc %d' % net_charge
        if resname:
            cmd += ' -rn %s' % molecule_name

        if log_debug_output: logger.debug(cmd)
        output = getoutput(cmd)
        if not os.path.exists('out.mol2'):
            msg  = "antechamber failed to produce output mol2 file\n"
            msg += "command: %s\n" % cmd
            msg += "output:\n"
            msg += 8 * "----------" + '\n'
            msg += output
            msg += 8 * "----------" + '\n'
            msg += "input mol2:\n"
            msg += 8 * "----------" + '\n'
            msg += read_file_contents(local_input_filename)
            msg += 8 * "----------" + '\n'
            raise Exception(msg)
        if log_debug_output: logger.debug(output)

        # Run parmchk.
        cmd = "parmchk2 -i out.mol2 -f mol2 -o out.frcmod"
        if log_debug_output: logger.debug(cmd)
        output = getoutput(cmd)
        if not os.path.exists('out.frcmod'):
            msg  = "parmchk2 failed to produce output frcmod file\n"
            msg += "command: %s\n" % cmd
            msg += "output:\n"
            msg += 8 * "----------" + '\n'
            msg += output
            msg += 8 * "----------" + '\n'
            msg += "input mol2:\n"
            msg += 8 * "----------" + '\n'
            msg += read_file_contents('out.mol2')
            msg += 8 * "----------" + '\n'
            raise Exception(msg)
        if log_debug_output: logger.debug(output)
        check_for_errors(output)

        #Copy back
        shutil.copy( 'out.mol2', gaff_mol2_filename )
        shutil.copy( 'out.frcmod', frcmod_filename )

    return gaff_mol2_filename, frcmod_filename


def run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename, prmtop_filename=None, inpcrd_filename=None, log_debug_output=False):
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
    log_debug_output : bool, optional, default=False
        If true, will send output of tleap to logger.

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

    #Get absolute paths for input/output
    gaff_mol2_filename = os.path.abspath( gaff_mol2_filename )
    frcmod_filename = os.path.abspath( frcmod_filename )
    prmtop_filename = os.path.abspath( prmtop_filename )
    inpcrd_filename = os.path.abspath( inpcrd_filename )

    #Work in a temporary directory, on hard coded filenames, to avoid any issues AMBER may have with spaces and other special characters in filenames
    with mdtraj.utils.enter_temp_directory():
        shutil.copy( gaff_mol2_filename, 'file.mol2' )
        shutil.copy( frcmod_filename, 'file.frcmod' )

        tleap_input = """
    source oldff/leaprc.ff99SB
    source leaprc.gaff
    LIG = loadmol2 file.mol2
    check LIG
    loadamberparams file.frcmod
    saveamberparm LIG out.prmtop out.inpcrd
    quit

"""

        file_handle = open('tleap_commands', 'w')
        file_handle.writelines(tleap_input)
        file_handle.close()

        cmd = "tleap -f %s " % file_handle.name
        if log_debug_output: logger.debug(cmd)

        output = getoutput(cmd)
        if log_debug_output: logger.debug(output)

        check_for_errors( output, other_errors = ['Improper number of arguments'] )

        #Copy back target files
        shutil.copy( 'out.prmtop', prmtop_filename )
        shutil.copy( 'out.inpcrd', inpcrd_filename )

    return prmtop_filename, inpcrd_filename
