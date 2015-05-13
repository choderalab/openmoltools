def read_next_topo_section(textarray):
   """Process a text array containing contents of a GROMACS topology file, and return a 
    text array consisting of the header and contents of the next [ named ] section in the topology file, as well as a second text array consisting of anything prior to the next named section.

    Parameters
    ----------
    textarray : list(str)
        Lines from a GROMACS topology file to be processed.

    Returns
    -------
    section : list(str)
        Contents of next section, including line with the header for the section
    extra : list(str)
        Contents of the input textarray prior to the start of the next section

    Notes
    -----
    Reads the [ named ] section in the GROMACS topology file from the passed textarray 
    and returns an array containing only that section (stopping before the subsequent 
    section). len(section) will be the number of lines in that section. Also returns a 
    second array containing that portion of textarray prior to the beginning of the next
    section.
    """

   #Pattern to match for start
   secstart=re.compile(r'\s*\[.+\]')

   #I think this is actually robust for comments: Comments shouldn't be recognized
   #by the search string above (nor should any bracketed object not preceded by)
   #spaces or nothing and so commented out sections will be treated as 'headers' and left untouched.

   #Find matches
   section=[]
   foundstart=False
   startline=0
   ctr=0
   for line in textarray:
     #search for match
     m=secstart.match(line)
     #If there is a match
     if m:
       #If this is the first match
       if not foundstart:
          #We have found the start of the section we want
          foundstart=True
          #Save this line
          section.append(line)
          startline=ctr
          continue #Go back to the top of the loop (the next line)
       #If this is not the first match
       elif foundstart:
          #Don't keep going into the next section
          break
     #Otherwise if this is not a match, so it doesn't start or end a section.
     #Only save if we already got to the section.
     elif foundstart:
       section.append(line)
     ctr+=1
   #End loop over array; we've now got the section.

   #Return the extra portion also
   extra=textarray[0:startline]
   return section,extra


def merge_topologies( input_topologies, output_topology, system_name, molecule_names = None, molecule_numbers = None ):
    """Merge GROMACS topology files specified in a list of input topologies and write the result into a final topology file specified. Optionally specify a list of molecule names to be used in the [ moleculetype ] and [ molecules ] sections, overriding what is already present.

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

    #
     
