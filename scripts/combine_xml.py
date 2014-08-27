#!/usr/bin/env python

from lxml import etree as ET
import os, sys
from collections import OrderedDict

"""     
    Combine several OpenMM XML files into one. 
    Author: Lee-Ping Wang

    Assumes that residue names and atom types 
    are non-overlapping, so make sure of that first! 

    Place all XML files into a "Molecules" 
    subfolder first in no particular order.
"""

# Using this parser removes abnormalities with indentation.
parser = ET.XMLParser(remove_blank_text=True)

# The new XML object.
newxml = None

# A dictionary of second-level elements in the new XML file (AtomTypes, Residues, etc.).
new_elems = {}

# Dictionary of lists of existing third-level elements in each second-level element.
# For speed, third-level elements are stored as string representations of the attribute dictionary.
# Used to check for existing elements before adding them.
have_attribs = OrderedDict()

# Dictionary of lists of existing atom-class combinations in second-level elements
# corresponding to class-based forces (HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce).
# Used to check for attempts to override existing elements that have the same atom class combinations.
have_classes = OrderedDict()

for fnm in os.popen('find Molecules -name \*.xml | xargs ls -v').readlines():
    # Loop over XML files in Molecules subdirectory.
    fnm = fnm.strip()
    if newxml == None:
        # If the first XML file, parse and record existing elements.
        newxml = ET.parse(fnm, parser)
        for mainlvl in newxml.getroot():
            new_elems[mainlvl.tag] = mainlvl
            for sublvl in mainlvl:
                substr = str(OrderedDict([(k, sublvl.attrib[k]) for k in sorted(sublvl.attrib.keys())]))
                have_attribs.setdefault(mainlvl.tag, []).append(substr)
                have_classes.setdefault(mainlvl.tag, []).append([sublvl.attrib[i] for i in sorted([k for k in sublvl.attrib.keys() if 'class' in k])])
    else:
        # Otherwise, parse and add new elements if any.
        addxml = ET.parse(fnm, parser)
        root = addxml.getroot()
        for mainlvl in root:
            for sublvl in mainlvl:
                substr = str(OrderedDict([(k, sublvl.attrib[k]) for k in sorted(sublvl.attrib.keys())]))
                if mainlvl.tag in ["AtomTypes", "Residues", "NonbondedForce"]:
                    if substr in have_attribs[mainlvl.tag]:
                        # We don't expect elements in AtomTypes, Residues, NonbondedForce to be duplicated.
                        raise RuntimeError('Found a duplicate element in AtomTypes / Residues / NonbondedForce!')
                    new_elems[mainlvl.tag].append(sublvl)
                    have_attribs[mainlvl.tag].append(substr)
                else:
                    if substr in have_attribs[mainlvl.tag]:
                        # On the other hand, the Force seconds contain lots of duplicate elements.
                        continue
                    else:
                        # Notify the user whenever new elements are added.
                        classes = [sublvl.attrib[i] for i in sorted([k for k in sublvl.attrib.keys() if 'class' in k])]
                        if classes in have_classes.values():
                            print sublvl.attrib, "\x1b[1;91mskipping because exist\x1b[0m", have_classes[classes]
                        else:
                            print sublvl.attrib, "\x1b[1;92madding\x1b[0m"
                            new_elems[mainlvl.tag].append(sublvl)
                            have_attribs[mainlvl.tag].append(substr)
                            have_classes[mainlvl.tag].append(classes)
                           
# Write final output. 
with open("combined.xml", "w") as f:
    newxml.write(f, pretty_print=True)

