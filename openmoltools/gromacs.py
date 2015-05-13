import os

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
        Line indices within lines belonging to section

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
    for end_index in range(start_index, nlines):
       # get line
       line,comments = stripcomments(lines[end_index])
       # split into elements
       elements = line.split()
       # see if keyword is matched
       if (len(elements) == 3):
          if (elements[0]=='['):
             fnd_end = True
             break
    if not fnd_end: end_index = nlines

    # compute indices of lines in section
    indices = range(start_index, end_index)

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
    Currently will NOT correctly handle a case where there might be multiple copies of a single section (i.e. multiple dihedrals sections) within a single molecule. Second and additional such sections will be ignored.
    """

    #PRELIMINARIES
    N_tops = len( input_topologies ) 

    #Check for obvious input problems - do we have the right number of everything, do all the input files exist
    if molecule_names <> None:
        assert len( molecule_names ) == N_tops, "Must provide same number of molecule names as topology files." 
    if molecule_numbers <> None:
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
            for index in indices:
                line, comments = stripcomments( thistop[index] )
                #If there is stuff here, store it if not already present
                if len(line) > 0:
                    if thistop[index] not in section_contents[sec]:
                        section_contents[sec].append( thistop[index] )
                    #If it's the [ molecules ] section, append even if duplicate since we'll change the names later
                    elif sec=='molecules':
                        section_contents[sec].append( thistop[index] )
                    #But if it's defaults and it's already there, check if it's OK
                    elif sec=='defaults':
                        tmp = line.split()
                        for entry in section_contents[sec]:
                            line2, comments2 = stripcomments( entry )
                            if len(line2) > 2:
                                tmp2 = line2.split()
                                if tmp2 <> tmp:
                                    raise ValueError('Non-equivalent defaults entries in topology files; unsure how to proceed. Offending entries are %s and %s.' % (line, line2) )
                    #If it's moleculetype and it's already there, check that the exclusions are OK
                    elif sec=='moleculetype':
                        tmp = line.split()
                        for entry in section_contents[sec]:
                            line2, comments2 = stripcomments( entry )
                            if len(line2) > 1:
                                if tmp[1] <> line2.split()[1]:
                                    raise ValueError('Non-equivalent number of exclusions in molecule definitions; unsure how to proceed. Offending entries are %s and %s." % (line, line2) )')
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
                    if sec=='moleculetype' and molecule_names <> None and len( line ) >1:
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
    if molecule_names <> None:
        if molecule_numbers <> None: 
            #If molecule numbers are also provided, build fully from scratch
            section_contents['molecules'] = []
            for topnr in range(N_tops):
                section_contents['molecules'].append( '%s    %s\n' % (molecule_names[topnr], molecule_numbers[topnr] ) )
                
        else: #If they're not provided, use existing numbers
            new_contents = []
            for line in section_contents['molecules']:
                entry, comments = stripcomments(line)
                if len(entry) > 1:
                    new_contents.append( '%s    %s\n' % (molecule_names[topnr], entry.split()[1] ) )         
            section_contents['molecules'] = new_contents 

    #If names are not provided, use what we already have
    else:
        if molecule_numbers <> None:
            #If molecule numbers are provided, build with existing names and new numbers
            new_contents = []
            for line in section_contents['molecules']:
                entry, comments = stripcomments(line)
                if len(entry) > 1:
                    new_contents.append( '%s    %s\n' % ( entry.split()[0], molecule_numbers[topnr] ) )
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
        topology_lines.append( '[ %s ]\n' % sec )
        for line in section_contents[sec]:
            topology_lines.append( line ) 
        topology_lines.append('\n')

    #Now handle stuff for the individual molecules present
    for topnr in range(N_tops):
        for secname in molecule_sections:
            topology_lines.append( '[ %s ]\n' % secname )
            topology_lines += topology_sections[ secname ] [ topnr ]          


    #Now add system name and molecules sections
    topology_lines.append( '[ system ]\n')
    topology_lines.append( '%s\n' % system_name )
    topology_lines.append('\n') 
    topology_lines.append( '[ molecules ]\n' )
    for line in section_contents['molecules']:
        topology_lines.append( line )
    topology_lines.append('\n') 

    #Write final topology file
    file = open( output_topology, 'w')
    file.writelines( topology_lines )
    file.close()
    
    return True
    
