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
