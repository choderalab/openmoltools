import os
import shutil
import logging
import mdtraj.utils
from distutils.spawn import find_executable
import parmed

from openmoltools.utils import getoutput

logger = logging.getLogger(__name__)

GROMACS_PATH = find_executable('gmx')


def stripcomments(line):
    """From a GROMACS topology formatted line, return (line, comments) with whitespace and comments stripped. Comments are given with ;.

    Parameters
    ----------
    line : str
        GROMACS line to be stripped

    Returns
    -------
    line : str
        GROMACS line with comments and whitespace (leading/trailing) stripped
    """
    # strip comments
    index = line.find(';')
    comments =''
    if index > -1:
       comments = line[index:]
       line = line[0:index]
    # strip whitespace
    line = line.strip()
    comments = comments.strip()
    # return stripped line
    return line,comments


def extract_section(lines, section):
    """Identify lines associate with a GROMACS topology section.

    Parameters
    ----------
    lines : list (str)
        the lines in the file (or portion of a file) to search
    section : str
        The section name to locate, without formatters, i.e. "atoms"

    Returns
    -------
    status : bool
        Whether or not section is found. True if found, False if not. 
    indices :  list (int)
        Line indices within lines belonging to section excluding the header and counting from zero

    """

    indices = list()

    status = False

    nlines = len(lines)
    for start_index in range(nlines):
       # get line
       line,comments = stripcomments(lines[start_index])
       # split into elements
       elements = line.split()
       # see if keyword is matched
       if (len(elements) == 3):
          if (elements[0]=='[') and (elements[1]==section) and (elements[2]==']'):
             # increment counter to start of section data and abort search
             start_index += 1
             status = True
             break

    # print if section not found (for debugging)
    if (start_index == nlines):
       print( "Section %(section)s not found." % vars())
       status = False

    # Locate end of section.
    fnd_end = False
    for idx in range(start_index, nlines):
       # get line
       line , comments = stripcomments(lines[idx])
       # split into elements
       elements = line.split()
       # see if keyword is matched
       if (len(elements) == 3):
          if (elements[0]=='['):
             end_index = end_index + 1
             fnd_end = True
             break
       if elements != []:
           end_index = idx + 1

    if not fnd_end:
        end_index=end_index + 1

    # compute indices of lines in section. The indices are related to the section without the header
    indices = range(start_index, end_index-1)

    # return these indices
    return status, indices

def change_molecules_section( input_topology, output_topology, molecule_names, molecule_numbers):
    """Create a GROMACS topology file where the  molecule numbers are replaced by new molecule numbers in the gromacs [ molecules ] section.
        
    Parameters
    ----------
    input_topology : str
        The input gromacs topology file name
    output_topology : str
        Output topology file names to be written/created
    molecule_names : list (str)
        Molecule names to be searched in the gromacs [ molecules ] section
    molecule_numbers : list (int)
        The new molecule numbers to be used in [ molecules ] section
        
    Returns
    -------
    nothing is returned
        
    Notes
    -----
    This function reads in a gromacs topology file and changes the number of atoms related to the passed molecule name list. 
    If in the topology file one molecule name is not present in the passed molecule name list an exception is raised.
    Currently assumes the components are single-residue, single-molecule (i.e. the molecule names and residue names are equivalent).
    """
    
    #The molecule name list and the molecule number list must have the same size otherwise an exception is raised
    assert len(molecule_names) == len(molecule_numbers), "The molecule name list and the molecule name number must have the same size"
    
    #Check for non negative integer number of molecules
    check_nni = all(item >=0  and isinstance(item, int) for item in molecule_numbers)

    if not check_nni:
        raise ValueError("The molecule number list must contain only non-negative integer value")

    #Read in the topology file as ParmEd object
    top = parmed.load_file( input_topology )

    #Split the topology file to its component molecules
    components = top.split()

    #Get objects and quantities
    molecules = []
    numbers = []
    current_names = []
    for c in components:
        molecules.append( c[0] )
        numbers.append( c[1] )    
        current_names.append( c[0].residues[0].name )
    
    #Check length
    assert len(molecules) == len(molecule_numbers), "The number of molecules in the topology file is not equal to the number of molecule numbers/molecules provided."


    #Cross-check current molecule names with expected molecule names
    for idx in range( len(molecule_names) ):
        if not molecule_names[idx] == current_names[idx]:
            raise ValueError("Molecule name provided (%s) does not match molecule/residue name in topology file (%s)." % (molecule_names[idx], current_names[idx] ) )

    #Re-compose topology file
    newtop = molecules[0] * molecule_numbers[0]
    for idx in range( 1, len(molecule_names) ):
        newtop += molecules[idx] * molecule_numbers[idx]
    

    #Write topology file
    newtop.write( output_topology ) 


def do_solvate( top_filename, gro_filename, top_solv_filename, gro_solv_filename, box_dim, box_type, water_model, water_top, FF = 'amber99sb-ildn.ff' ):

    """ This function creates water solvated molecule coordinate files and its corresponding topology
        
        PARAMETERS:
            top_filename: str
                          Topology path/filename
            gro_filename: str
                          Coordinates path/filename
            top_solv_filename: str
                          Topology path/filename (solvated system)
            gro_solv_filename: str
                          Coordinates path/filename (solvated system)
            box_dim: float
                          cubic box dimension (nm); will be passed to GROMACS wiin .2f format
            box_type: str
                          box type (string passed to gmx solvate)
            water_model: str
                          Water model string to tell gmx solvate to use when solvating, i.e. "spc216"
            water_top: str
                          Water include file to ensure is present in topology file, i.e. "tip3p.itp"
            FF : str, optional, default = 'amber99sb-ildn.ff'
                          String specifying base force field directory for include files (i.e. 'amber99sb-ildn.ff').  

        NOTES:
        -----
        Primarily tested on 3 point water models. May need adjustment for other models.
"""

    #Setting up proper environment variables
    os.environ['GMX_MAXBACKUP'] = '-1' # Avoids unnecessary GROMACS backup files
    os.environ['GMX_NO_QUOTES'] = '1' # Supresses end-of-file quotes (gcq)

    #Get absolute paths for input/output
    top_filename = os.path.abspath( top_filename )
    gro_filename = os.path.abspath( gro_filename )
    top_solv_filename = os.path.abspath( top_solv_filename )
    gro_solv_filename = os.path.abspath( gro_solv_filename )

    with mdtraj.utils.enter_temp_directory(): #Work on hard coded filenames in temporary directory

        shutil.copy( gro_filename, 'in.gro' )
        shutil.copy( top_filename, 'out.top' )

        #string with the Gromacs 5.0.4 box generating commands
        cmdbox = 'gmx editconf -f in.gro -o out.gro -c -d %.2f -bt %s' % ( box_dim, box_type)
        output = getoutput(cmdbox)
        logger.debug(output)
        check_for_errors(output)

        #string with the Gromacs 5.0.4 solvation tool (it is not genbox anymore)
        cmdsolv = 'gmx solvate -cp out.gro -cs %s -o out.gro -p out.top' % (water_model)
        output = getoutput(cmdsolv)
        logger.debug(output)
        check_for_errors(output)

        #Insert Force Field specifications
        ensure_forcefield( 'out.top', 'out.top', FF = FF)

        #Insert line for water topology portion of the code
        try:
            file = open('out.top','r')
            text = file.readlines()
            file.close()
        except:
            raise NameError('The file out.top is missing' )

        #Insert water model
        wateritp = os.path.join(FF, water_top ) # e.g water_top = 'tip3p.itp'
        index = 0
        while '[ system ]' not in text[index]:
            index += 1
        text[index] = '#include "%s"\n\n' % wateritp + text[index]

        #Write the topology file
        try:
            file = open('out.top','w+')
            file.writelines( text )
            file.close()
        except:
            raise NameError('The file %s is missing' % 'out.top')

        #Check if file exist and is not empty;
        if os.stat( 'out.gro' ) == 0 or os.stat( 'out.top' ).st_size == 0:
            raise(ValueError("Solvent insertion failed"))

        #Copy back files
        shutil.copy( 'out.gro', gro_solv_filename )
        shutil.copy( 'out.top', top_solv_filename )

    return

def ensure_forcefield( intop, outtop, FF = 'ffamber99sb-ildn.ff'):
    """Open a topology file, and check to ensure that includes the desired forcefield itp file. If not, remove any [ defaults ] section (which will be provided by the FF) and include the forcefield itp. Useful when working with files set up by acpypi -- these need to have a water model included in order to work, and most water models require a force field included in order for them to work.
            
            ARGUMENTS:
            - intop: Input topology
            - outtop: Output topology
        OPTIONAL:
        - FF: String corresponding to desired force field; default ffamber99sb.-ildn.ff
        
        Limitations:
        - If you use this on a topology file that already includes a DIFFERENT forcefield, the result will be a topology file including two forcefields.
    """
    
    file = open(intop, 'r')
    text= file.readlines()
    file.close()
    
    FFstring = FF+'/forcefield.itp'
    
    #Check if force field is included somewhere
    found = False
    for line in text:
        if FFstring in line:
            found = True
    #If not, add it after any comments at the top
    if not found:
        idx = 0
        while text[idx].find(';')==0:
            idx+=1
        text[idx] = '\n#include "%s"\n\n' % FFstring + text[idx]
    
    #Remove any defaults section
    found = False
    foundidx = None
    endidx = None
    for (idx, line) in enumerate(text):
        if '[ defaults ]' in line:
            foundidx = idx
            found = True
        #If we've already found the defaults section, find location of start of next section
        #Assumes next section can be found by looking for a bracket at the start of a line
        elif found and '[' in line:
            #Only store if we didn't already store
            if endidx == None:
                endidx = idx
    #Now remove defaults section
    
    if found:
        text = text[0:foundidx] + text[endidx:]
    
    
    #Write results
    file = open( outtop, 'w')
    file.writelines(text)
    file.close()


    
def check_for_errors( outputtext, other_errors = None, ignore_errors = None ):
    """Check GROMACS package output for the string 'ERROR' (upper or lowercase) and (optionally) specified other strings and raise an exception if it is found (to avoid silent failures which might be noted to log but otherwise ignored).

    Parameters
    ----------
    outputtext : str
        String listing output text from an (GROMACS) command which should be checked for errors.
    other_errors : list(str), default None
        If specified, provide strings for other errors which will be chcked for, such as "improper number of arguments", etc.
    ignore_errors: list(str), default None
        If specified, GROMACS output lines containing errors but also containing any of the specified strings will be ignored

    Notes
    -----
    If error(s) are found, raise a RuntimeError and attept to print the appropriate errors from the processed text. This only currently prints useful output if the word ERROR occurs on the same line as the cause of the error."""
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
        print("Unexpected errors encountered running GROMACS tool. Offending output:")
        for line in error_lines: print(line)
        raise(RuntimeError("Error encountered running GROMACS tool. Exiting."))

    return

def merge_topologies( input_topologies, output_topology, system_name, molecule_names = None, molecule_numbers = None ):
    """Merge GROMACS topology files specified in a list of input topologies and write the result into a final topology file specified. Optionally specify a list of molecule names to be used in the [ moleculetype ] and [ molecules ] sections, overriding what is already present. If molecule_names and/or molecule_numbers are specified, the input topologies are expected to be single-molecule, single-residue topology files.

    Parameters
    ----------
    input_topologies : list (str)
        A list of input topology files to be merged
    output_topology : str
        Output topology file to be written/created
    system_name : str
        Name to be used in final [ system ] section
    molecule_names : list (str), optional
        Molecule names to use in [ moleculetype ] and [ system ] sections in resulting topology file
    molecule_numbers : list (int), optional
        Molecule numbers to use in [ system ] section, overriding what is present in the existing [ system ] sections.

    Returns
    -------
    status : bool
        True if successful

    Notes
    -----
    This simply takes the contents of provided topology files and consolidates them to make a single resulting topology file. Existing [ moleculetype ] definitions are preserved and molecules are kept separate.
    Merging should work fine on general topology files. However, if molecule_names and/or molecule_numbers are provided the input topologies must contain single-molecule, single-residue topologies. Molecule_numbers are used to replicate the input topologies via ParmEd, and molecule_names are used to change the residue names from the input topologies.
    """

    #PRELIMINARIES
    N_tops = len( input_topologies )

    #Check for obvious input problems - do we have the right number of everything, do all the input files exist
    if molecule_numbers != None:
        assert len( molecule_numbers ) == N_tops, "Must provide same number of molecule numbers as topology files."
    for filenm in input_topologies:
        assert( os.path.isfile( filenm )), "Error: Can't find input file %s provided to merge_topologies." % filenm

    #WORK ON TOPOLOGIES
    tops = []
    for filenm in input_topologies:
        top = parmed.gromacs.GromacsTopologyFile( filenm )
        tops.append( top )

    #List numbers of each molecule if not provided
    if molecule_numbers == None:
        molecule_numbers = [ 1] * N_tops

    #Check that we've been provided with the correct number of molecule_names if any
    if molecule_names != None:
        total_molecules = 0
        for topnr in range(N_tops): 
            total_molecules += len( tops[topnr].residues ) 
        assert total_molecules == len( molecule_names ), "Must provide a number of molecule names equal to your total number of residues, but you have %s and %s, respectively." % ( len( molecule_names), total_molecules )

        #Rename residues
        ctr = 0
        for topnr in range(N_tops):
            for resnr in range(len(tops[topnr].residues)):
                tops[topnr].residues[resnr].name = molecule_names[ ctr ]
                ctr += 1
 
    #Construct final topology
    final = tops[0] * molecule_numbers[ 0 ] 
    for topnr in range( 1, N_tops ):
        final += tops[ topnr ] * molecule_numbers[ topnr ] 

    #Set system name
    final.title = system_name

    #Write topology
    parmed.gromacs.GromacsTopologyFile.write( final, output_topology ) 

    return True

 
