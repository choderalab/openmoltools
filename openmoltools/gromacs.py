import os
import shutil
import logging
try:
    from subprocess import getoutput  # If python 3
except ImportError:
    from commands import getoutput  # If python 2

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")


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


def merge_topologies( input_topologies, output_topology, system_name, molecule_names = None, molecule_numbers = None ):
    """Merge GROMACS topology files specified in a list of input topologies and write the result into a final topology file specified. Currently these are required to be SINGLE-MOLECULE topology files. Optionally specify a list of molecule names to be used in the [ moleculetype ] and [ molecules ] sections, overriding what is already present.

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
    This is not necessarily a particularly robust function, as it is expected to be replaced by functionality within ParmEd so we have not paid a great deal of attention to ensuring all possible topology file formats are handled.
    Provided topology files are required to be SINGLE-MOLECULE topology files.
    Include files are not handled carefully. Initial include files at the top of topology files will be removed. Include files occurring after a molecule definition may be retained, but there is no checking whether these are duplicated or placed properly in the final topology file.
    Currently will NOT correctly handle a case where there might be multiple copies of a single section (i.e. multiple dihedrals sections) within a single molecule. Second and additional such sections will be ignored.
    """

    #PRELIMINARIES
    N_tops = len( input_topologies ) 

    #Check for obvious input problems - do we have the right number of everything, do all the input files exist
    if molecule_names != None:
        assert len( molecule_names ) == N_tops, "Must provide same number of molecule names as topology files." 
    if molecule_numbers != None:
        assert len( molecule_numbers ) == N_tops, "Must provide same number of molecule numbers as topology files."
    for filenm in input_topologies:
        assert( os.path.isfile( filenm )), "Error: Can't find input file %s provided to merge_topologies." % filenm


    #Read in topology files
    topology_text = []
    for filenm in input_topologies:
        file = open(filenm, 'r')
        topology_text.append(file.readlines() )
        file.close()

    #We require single-molecule topology files. Check that these are.
    for topnr in range(N_tops):
        thistop = topology_text[ topnr ]
        found_moltype = False
        for entry in thistop:
            line, comments = stripcomments( entry )
            elements = line.split()
            #If this is a section name, check and see if it's a moleculetype 
            if len(elements) == 3:
                if (elements[0]=='[') and (elements[2]==']'):
                    secname = elements[1]
                    if secname=='moleculetype':
                        if not found_moltype: found_moltype = True
                        elif found_moltype:
                            raise ValueError("merge_topologies requires single-molecule topology files. Multiple moelculetype definitions are found in in topology %s. Halting to avoid creating erroneous topologies." % (input_topologies[topnr]) )


    #Store a list of sections processed so we make sure there are no sections we DON'T process (throw an exception if there is one we don't handle)
    sections_processed = [] 

    #Handle overall sections - atomtypes, defaults, moleculetype, (NOT sections for specific molecules - just those for the final whole system)
    #Note that `moleculetype` is only processed here for the purposes of checking exclusions; the individual moleculetypes will be handled later.
    section_names = [ 'defaults', 'atomtypes', 'moleculetype', 'molecules'  ]
    sections_processed += section_names + ['system']
    section_contents = {}
    #For every section, store the lines to a list if they are not already present (and formatted identically). If they are already present, skip. No error checking is done (i.e. duplicate parameters which are not formatted identically will be included, and conflicting parameter definitions will also be included), EXCEPT that defaults and moleculetype sections are checked to ensure that combination rules/fudgeLJ/fudgeQQ and exclusions match. An exception is raised if not.
    for sec in section_names:
        section_contents[sec] = []
        for topnr in range(N_tops):
            thistop = topology_text[topnr]
            status, indices = extract_section( thistop, sec )
            if status:
                for index in indices:
                    line, comments = stripcomments( thistop[index] )
                    if len(line) > 0:
                        #If it's the [ molecules ] section, append even if duplicate since we'll change the names later
                        if sec=='molecules':
                            section_contents[sec].append( thistop[index] )
                        #But if it's defaults and it's already there, check if it's OK
                        elif sec=='defaults':
                            tmp = line.split()
                            identical = None
                            for entry in section_contents[sec]:
                                line2, comments2 = stripcomments( entry )
                                if len(line2) > 2:
                                    tmp2 = line2.split()
                                    for idx in range(len(tmp)):
                                        if tmp[idx] != tmp2[idx]:
                                            identical = False
                                            raise ValueError('Non-equivalent defaults entries in topology files; unsure how to proceed. Offending entries are %s and %s.' % (line, line2) )
                                        else:
                                            identical = True
                            if not identical:
                                section_contents[sec].append( thistop[index] )
                            
                        #If it's moleculetype and it's already there, check that the exclusions are OK
                        elif sec=='moleculetype':
                            tmp = line.split()
                            identical = None
                            for entry in section_contents[sec]:
                                line2, comments2 = stripcomments( entry )
                                if len(line2) > 1:
                                    tmp2 = line2.split()
                                    if tmp[1] != tmp2[1]:
                                        identical = False
                                        raise ValueError('Non-equivalent number of exclusions in molecule definitions; unsure how to proceed. Offending entries are %s and %s." % (line, line2) )')
                                    else:
                                        identical = True
                            if not identical:
                                section_contents[sec].append( thistop[index] )

                        #For everything else, if there is stuff here, store it if not already present
                        elif thistop[index] not in section_contents[sec]:
                            section_contents[sec].append( thistop[index] )
                    #If it's just comments, store it if it's not already there
                    elif len(comments) > 0:
                        if thistop[index] not in section_contents[sec]:
                            section_contents[sec].append( thistop[index] )

    #Now we've handled all the generic sections. Next handle all the per-molecule sections, tracking which topology they came from
    molecule_sections = [ 'moleculetype', 'atoms', 'bonds', 'pairs', 'angles', 'dihedrals' ]
    sections_processed += molecule_sections
    topology_sections = {}
    for sec in molecule_sections:
        topology_sections[sec] = {}
        for topnr in range(N_tops):
            topology_sections[sec][topnr] = []
            thistop = topology_text[ topnr ]
            status, indices = extract_section( thistop, sec )
            if status:
                for idx in indices[1:]:
                    line, comments = stripcomments( thistop[idx] )
                    #If this is a non-commented line in a moleculetype section and we were provided with names, we need to change the name
                    if sec=='moleculetype' and molecule_names != None and len( line ) >1:
                        newline = '%s         %s\n' % ( molecule_names[topnr], line.split()[1] )
                        topology_sections[sec][topnr].append( newline ) 
                    else:
                        topology_sections[sec][topnr].append( thistop[idx] )


    #Check sections_processed against the actual sections present in the topology file and throw an exception if there are sections present in the topology file which we don't have here
    for topnr in range(N_tops):
        thistop = topology_text[ topnr ]
        for entry in thistop:
            line, comments = stripcomments( entry )
            elements = line.split()
            #If this is a section name, make sure it's one we know...
            if len(elements) == 3:
                if (elements[0]=='[') and (elements[2]==']'):
                    secname = elements[1]
                    if not secname in sections_processed:
                        raise ValueError("Unprocessed section in topology %s; section name is %s. Halting to avoid creating erroneous topologies." % (input_topologies[topnr], secname) ) 

  
     
    #Construct final molecules section
    #If names are provided, construct section at least partially from scratch
    if molecule_names != None:
        if molecule_numbers != None: 
            #If molecule numbers are also provided, build fully from scratch
            section_contents['molecules'] = []
            for topnr in range(N_tops):
                section_contents['molecules'].append( '%s    %s\n' % (molecule_names[topnr], molecule_numbers[topnr] ) )
                
        else: #If they're not provided, use existing numbers
            new_contents = []
            topnr = 0
            for line in section_contents['molecules']:
                entry, comments = stripcomments(line)
                if len(entry) > 1:
                    new_contents.append( '%s    %s\n' % (molecule_names[topnr], entry.split()[1] ) )  
                    topnr += 1 
            section_contents['molecules'] = new_contents 

    #If names are not provided, use what we already have
    else:
        if molecule_numbers != None:
            #If molecule numbers are provided, build with existing names and new numbers
            new_contents = []
            topnr = 0
            for line in section_contents['molecules']:
                entry, comments = stripcomments(line)
                if len(entry) > 1:
                    new_contents.append( '%s    %s\n' % ( entry.split()[0], molecule_numbers[topnr] ) )
                    topnr +=1
            section_contents['molecules'] = new_contents

        #Otherwise we do nothing - we're just using existing section

    #Check that names in molecules are unique
    used_names = []
    for line in section_contents['molecules'][1:]:
        name = line.split()[0]
        if name in used_names:
            raise ValueError("Duplicate name in final molecules section, which will result in an incorrect topology file. Duplicate is %s. Halting." % name )
        used_names.append( name )

    #Now build final topology file
    topology_lines = []
    #First handle stuff which occurs only once at the top of the topology file
    for sec in [ 'defaults', 'atomtypes']:
        #Only store these if the sections aren't empty - i.e. we might not have a defaults or atomtypes section.
        if len( section_contents[sec] ) > 0:
            topology_lines.append( '[ %s ]\n' % sec )
            for line in section_contents[sec]:
                topology_lines.append( line ) 
            topology_lines.append('\n')

    #Now handle stuff for the individual molecules present
    for topnr in range(N_tops):
        for secname in molecule_sections:
            topology_lines.append( '[ %s ]\n' % secname )
            topology_lines += topology_sections[ secname ] [ topnr ]          
            #Section extraction strips whitespace; rebuild newline at end
            topology_lines.append('\n')

    #Now add system name and molecules sections
    topology_lines.append( '[ system ]\n')
    topology_lines.append( '%s\n' % system_name )
    topology_lines.append('\n') 
    topology_lines.append( '[ molecules ]\n' )
    for line in section_contents['molecules']:
        topology_lines.append( line )

    #Write final topology file
    file = open( output_topology, 'w')
    file.writelines( topology_lines )
    file.close()
    
    return True

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
    
    """
    
    #The molecule name list and the molecule number list must have the same size otherwise an exception is raised
    if (len(molecule_names) != len(molecule_numbers)):
        raise ValueError("The molecule name list and the molecule name number must have the same size")
    
    #Check for non negative integer number of molecules
    check_nni = all(item >=0  and isinstance(item, int) for item in molecule_numbers)

    if not check_nni:
        raise ValueError("The molecule number list contains a non negative integer value")

    #Read in the topology file
    try:
        file = open(input_topology, 'r')
        text = file.readlines()
        file.close()
    except IOError:
        raise NameError(input_topology)


    #Extract the [ molecules ] section from the gromacs topolgy file
    check , indices = extract_section(text, "molecules")

    if(check == False):
        raise ValueError("In the selected gromacs topology file is missing the [ molecules ] section")

    counter = 0


    for idx in indices:

        #Read in the molecule name and the molecule number present in the topology file
        line = text[idx]		
        line = line.split()
   
        #Skip comments
        if line[0] == ';':
            continue
        
        #The current name of the molecule in the topology file
        mol_name = line[0]

        #Check if the molecule name in the gromacs topology is present in the passed molecule name list otherwise an exception is raised
        try:
            index_name = molecule_names.index(mol_name)
        except ValueError:
            msg = "In the selected gromacs topology file the molecule name '%s' does not match any molecule names in the passed molecule name list: %s" % (mol_name , molecule_names)
            raise ValueError(msg)
        #Create a valid string to re-write back in the output topology file
        line[1] = str(molecule_numbers[index_name])
        strg = ' '.join(line)+'\n'
        text[idx] = strg
        counter=counter+1

    if (counter != len(molecule_names)):
        raise ValueError("The passed molecule list contains more molecule names than the gromacs topology file")
    
    try:
        file = open( output_topology, 'w')
        file.writelines( text )
        file.close()
    except IOError:
        msg = "It was not possible to write the output topology file: %s" % output_topology
        raise NameError(msg)


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
                          cubic box dimension (nm); will be passed to GROMACS with one digit of precision (%3.1f)
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

    #Setting up proper environment variable (avoid unnecessary GROMCAS backup files)
    os.environ['GMX_MAXBACKUP'] = '-1'


    #copies topology file to new directory
    shutil.copyfile(top_filename, top_solv_filename) 

    #string with the Gromacs 5.0.4 box generating commands
    cmdbox = 'gmx editconf -f %s -o %s -c -d %3.1f -bt %s' % (gro_filename, gro_solv_filename, box_dim, box_type)
    output = getoutput(cmdbox)
    logger.debug(output)

    #string with the Gromacs 5.0.4 solvation tool (it is not genbox anymore)
    cmdsolv = 'gmx solvate -cp %s -cs %s -o %s -p %s' % (gro_solv_filename, water_model, gro_solv_filename, top_solv_filename)
    output = getoutput(cmdsolv)
    logger.debug(output)

    #Insert Force Field specifications
    ensure_forcefield( top_solv_filename, top_solv_filename, FF = FF)

    #Insert line for water topology portion of the code
    try:
        file = open(top_solv_filename,'r')
        text = file.readlines()
        file.close()
    except:
        raise NameError('The file %s is missing' % top_solv_filename)

    #Insert water model
    wateritp = os.path.join(FF, water_top ) # e.g water_top = 'tip3p.itp'
    index = 0
    while '[ system ]' not in text[index]:
        index += 1
    text[index] = '#include "%s"\n\n' % wateritp + text[index]

    #Write the topology file
    try:
        file = open(top_solv_filename,'w+')
        file.writelines( text )
        file.close()
    except:
        raise NameError('The file %s is missing' % top_solv_filename)

    #Check if file exist and is not empty;
    if os.stat( gro_solv_filename ) == 0 or os.stat( top_solv_filename ).st_size == 0:
        raise(ValueError("Solvent insertion failed"))

    return

