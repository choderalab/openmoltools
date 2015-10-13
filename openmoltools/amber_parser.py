#!/usr/bin/env python
import sys
import math
import simtk.openmm.app.element as element
import simtk.unit as unit
import subprocess
import datetime
from six.moves import cStringIO
import mdtraj as md

import logging
logger = logging.getLogger(__name__)


def fix(atomClass):
    if atomClass == 'X':
        return ''
    return atomClass

elements = {}
for elem in element.Element._elements_by_symbol.values():
    num = elem.atomic_number
    if num not in elements or elem.mass < elements[num].mass:
        elements[num] = elem

OTHER = 0
ATOMS = 1
CONNECT = 2
CONNECTIVITY = 3
RESIDUECONNECT = 4
section = OTHER

charge14scale = 1.0 / 1.2
epsilon14scale = 0.5

skipResidues = ['CIO', 'IB']  # "Generic" ions defined by Amber, which are identical to other real ions
skipClasses = ['OW', 'HW']  # Skip water atoms, since we define these in separate files


class AmberParser(object):

    def __init__(self, override_mol2_residue_name=None):
        """Create an AmberParser object for converting amber force field files to XML format.
        
        Parameters
        ----------
        override_mol2_residue_name : str, default=None
            If given, use this name to override mol2 residue names.
            Useful to ensure that multiple ligands have unique residue
            names, as required by the OpenMM ffXML parser.
        """
        
        self.override_mol2_residue_name = override_mol2_residue_name
        self.current_mol2 = 0
        
        self.residueAtoms = {}
        self.residueBonds = {}
        self.residueConnections = {}

        self.types = []
        self.type_names = []
        self.masses = {}
        self.resAtomTypes = {}
        self.vdwEquivalents = {}
        self.vdw = {}
        self.charge = {}
        self.bonds = []
        self.angles = []
        self.torsions = []
        self.impropers = []

        self.set_provenance()

    def addAtom(self, residue, atomName, atomClass, element, charge, use_numeric_types=True):
        """Add an atom to the database of FF data.

        Notes
        -----

        use_numeric_types was not originally present in the OpenMM AMBER
        parsers.  It was added so that we can have atom types of the form
        "RES-X", where RES is the name of the molecule or residue and X
        is the atom numbering within that molecule.  use_numeric_types is
        set to False when processing mol2 files--e.g. for ligands.
        """
        if residue is None:
            return
        type_id = len(self.types)
        self.residueAtoms[residue].append([atomName, type_id])
        self.types.append((atomClass, element, charge))
        if use_numeric_types:
            self.type_names.append("%d" % (type_id))
        else:
            self.type_names.append("%s-%s" % (residue, atomName))

    def addBond(self, residue, atom1, atom2):
        """Add a bond to the database of FF data."""
        if residue is None:
            return
        self.residueBonds[residue].append((atom1, atom2))

    def addExternalBond(self, residue, atom):
        """Add an external bond to the database of FF data."""
        if residue is None:
            return
        if atom != -1:
            self.residueConnections[residue] += [atom]

    def process_mol2_file(self, inputfile):
        """Process an AMBER GAFF-compatible mol2 file.

        Parameters
        ----------
        inputfile : str
            filename of an .mol2 file

        Notes
        -----

        Antechamber is known to produce NONSTANDARD mol2 files.  This function
        is designed to work with those nonstandard mol2 files, not
        Tripos standard mol2 files.  We are forced to live with the poor
        decisions of our predecessors...

        """
        atoms, bonds = md.formats.mol2.mol2_to_dataframes(inputfile)
        
        if self.override_mol2_residue_name is None:
            residue_name = atoms.resName[1]  # To Do: Add check for consistency
        else:
            residue_name = self.override_mol2_residue_name
        
        # Give each mol2 file a unique numbering to avoid conflicts.
        residue_name = "%s-%d" % (residue_name, self.current_mol2)
        self.current_mol2 += 1

        self.residueAtoms[residue_name] = []
        self.residueBonds[residue_name] = []
        self.residueConnections[residue_name] = []

        for (i0, i1, name, x, y, z, atype, code, resname, charge) in atoms.itertuples(index=True):
            # i0 and i1 are zero-based and one-based indices, respectively
            full_name = residue_name + "_" + name
            element_symbol = md.formats.mol2.gaff_elements[atype]
            e = element.Element.getBySymbol(element_symbol)
            self.addAtom(residue_name, name, atype, e, charge, use_numeric_types=False)  # use_numeric_types set to false to use string-based atom names, rather than numbers
            self.vdwEquivalents[full_name] = atype

        for (id0, id1, bond_type) in bonds.itertuples(False):
            i = id0 - 1  # Subtract 1 for zero based indexing in OpenMM???
            j = id1 - 1  # Subtract 1 for zero based indexing in OpenMM???
            self.addBond(residue_name, i, j)

    def process_library_file(self, inputfile):
        """Process an AMBER .lib file.

        Parameters
        ----------
        inputfile : str
            filename of an .lib file

        """
        for line in open(inputfile):
            if line.startswith('!entry'):
                fields = line.split('.')
                residue = fields[1]
                if residue in skipResidues:
                    residue = None
                    continue
                key = fields[3].split()[0]
                if key == 'atoms':
                    section = ATOMS
                    self.residueAtoms[residue] = []
                    self.residueBonds[residue] = []
                    self.residueConnections[residue] = []
                elif key == 'connect':
                    section = CONNECT
                elif key == 'connectivity':
                    section = CONNECTIVITY
                elif key == 'residueconnect':
                    section = RESIDUECONNECT
                else:
                    section = OTHER
            elif section == ATOMS:
                fields = line.split()
                atomName = fields[0][1:-1]
                atomClass = fields[1][1:-1]
                if fields[6] == '-1':
                    # Workaround for bug in some Amber files.
                    if atomClass[0] == 'C':
                        elem = elements[6]
                    elif atomClass[0] == 'H':
                        elem = elements[1]
                    else:
                        raise ValueError('Illegal atomic number: ' + line)
                else:
                    elem = elements[int(fields[6])]
                self.charge = float(fields[7])
                self.addAtom(residue, atomName, atomClass, elem, self.charge)
            elif section == CONNECT:
                self.addExternalBond(residue, int(line) - 1)
            elif section == CONNECTIVITY:
                fields = line.split()
                self.addBond(residue, int(fields[0]) - 1, int(fields[1]) - 1)
            elif section == RESIDUECONNECT:
                # Some Amber files have errors in them, incorrectly listing atoms that should not be
                # connected in the first two positions.  We therefore rely on the "connect" section for
                # those, using this block only for other external connections.
                for atom in [int(x) - 1 for x in line.split()[2:]]:
                    self.addExternalBond(residue, atom)

    def process_dat_file(self, inputfile):
        """Process an AMBER .dat file.

        Parameters
        ----------
        inputfile : str
            filename of an .dat file

        """
        block = 0
        continueTorsion = False
        for line in open(inputfile):
            line = line.strip()
            if block == 0:     # Title
                block += 1
            elif block == 1:   # Mass
                fields = line.split()
                if len(fields) == 0:
                    block += 1
                else:
                    self.masses[fields[0]] = float(fields[1])
            elif block == 2:   # Hydrophilic atoms
                block += 1
            elif block == 3:   # Bonds
                if len(line) == 0:
                    block += 1
                else:
                    fields = line[5:].split()
                    self.bonds.append((line[:2].strip(), line[3:5].strip(), fields[0], fields[1]))
            elif block == 4:   # Angles
                if len(line) == 0:
                    block += 1
                else:
                    fields = line[8:].split()
                    self.angles.append((line[:2].strip(), line[3:5].strip(), line[6:8].strip(), fields[0], fields[1]))
            elif block == 5:   # Torsions
                if len(line) == 0:
                    block += 1
                else:
                    fields = line[11:].split()
                    periodicity = int(float(fields[3]))
                    if continueTorsion:
                        self.torsions[-1] += [float(fields[1]) / float(fields[0]), fields[2], abs(periodicity)]
                    else:
                        self.torsions.append([line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip(), float(fields[1]) / float(fields[0]), fields[2], abs(periodicity)])
                    continueTorsion = (periodicity < 0)
            elif block == 6:   # Improper torsions
                if len(line) == 0:
                    block += 1
                else:
                    fields = line[11:].split()
                    self.impropers.append((line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip(), fields[0], fields[1], fields[2]))
            elif block == 7:   # 10-12 hbond potential
                if len(line) == 0:
                    block += 1
            elif block == 8:   # VDW equivalents
                if len(line) == 0:
                    block += 1
                else:
                    fields = line.split()
                    for atom in fields[1:]:
                        self.vdwEquivalents[atom] = fields[0]
            elif block == 9:   # VDW type
                block += 1
                self.vdwType = line.split()[1]
                if self.vdwType not in ['RE', 'AC']:
                    raise ValueError('Nonbonded type (KINDNB) must be RE or AC')
            elif block == 10:   # VDW parameters
                if len(line) == 0:
                    block += 1
                else:
                    fields = line.split()
                    self.vdw[fields[0]] = (fields[1], fields[2])

    def process_frc_file(self, inputfile):
        """Process an AMBER .frc file.

        Parameters
        ----------
        inputfile : str
            filename of an .frc file

        """
        block = ''
        continueTorsion = False
        first = True
        for line in open(inputfile):
            line = line.strip()
            if len(line) == 0 or first:
                block = None
                first = False
            elif block is None:
                block = line
            elif block.startswith('MASS'):
                fields = line.split()
                self.masses[fields[0]] = float(fields[1])
            elif block.startswith('BOND'):
                fields = line[5:].split()
                self.bonds.append((line[:2].strip(), line[3:5].strip(), fields[0], fields[1]))
            elif block.startswith('ANGL'):
                fields = line[8:].split()
                self.angles.append((line[:2].strip(), line[3:5].strip(), line[6:8].strip(), fields[0], fields[1]))
            elif block.startswith('DIHE'):
                fields = line[11:].split()
                periodicity = int(float(fields[3]))
                if continueTorsion:
                    self.torsions[-1] += [float(fields[1]) / float(fields[0]), fields[2], abs(periodicity)]
                else:
                    self.torsions.append([line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip(), float(fields[1]) / float(fields[0]), fields[2], abs(periodicity)])
                continueTorsion = (periodicity < 0)
            elif block.startswith('IMPR'):
                fields = line[11:].split()
                self.impropers.append((line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip(), fields[0], fields[1], fields[2]))
            elif block.startswith('NONB'):
                fields = line.split()
                self.vdw[fields[0]] = (fields[1], fields[2])

    def generate_xml(self):
        """Return the processed forcefield files as an XML stream.

        Returns
        -------
        stream : cStringIO
            The text of the output XML forcefield data.

        Notes
        -----

        The stream can be written to disk via:

        outfile = open("my_forcefield.xml", 'w')
        outfile.write(stream.read())
        outfile.close()

        """
        stream = cStringIO()
        write_stream = lambda x: stream.write(x + "\n")
        write_stream(self.provenance)
        write_stream("<ForceField>")
        write_stream(" <AtomTypes>")
        for index, type in enumerate(self.types):
            write_stream("""  <Type name="%s" class="%s" element="%s" mass="%s"/>""" % (self.type_names[index], type[0], type[1].symbol, type[1].mass.value_in_unit(unit.amu)))
        write_stream(" </AtomTypes>")
        write_stream(" <Residues>")
        for res in sorted(self.residueAtoms):
            write_stream("""  <Residue name="%s">""" % res)
            for atom in self.residueAtoms[res]:
                atom_name, type_id = tuple(atom)
                atom_type = self.type_names[type_id]
                write_stream("   <Atom name=\"%s\" type=\"%s\"/>" % (atom_name, atom_type))
            if res in self.residueBonds:
                for bond in self.residueBonds[res]:
                    write_stream("""   <Bond from="%d" to="%d"/>""" % bond)
            if res in self.residueConnections:
                for bond in self.residueConnections[res]:
                    write_stream("""   <ExternalBond from="%d"/>""" % bond)
            write_stream("  </Residue>")
        write_stream(" </Residues>")
        write_stream(" <HarmonicBondForce>")
        processed = set()
        for bond in self.bonds:
            signature = (bond[0], bond[1])
            if signature in processed:
                continue
            if any([c in skipClasses for c in signature]):
                continue
            processed.add(signature)
            length = float(bond[3]) * 0.1
            k = float(bond[2]) * 2 * 100 * 4.184
            write_stream("""  <Bond class1="%s" class2="%s" length="%s" k="%s"/>""" % (bond[0], bond[1], str(length), str(k)))
        write_stream(" </HarmonicBondForce>")
        write_stream(" <HarmonicAngleForce>")
        processed = set()
        for angle in self.angles:
            signature = (angle[0], angle[1], angle[2])
            if signature in processed:
                continue
            if any([c in skipClasses for c in signature]):
                continue
            processed.add(signature)
            theta = float(angle[4]) * math.pi / 180.0
            k = float(angle[3]) * 2 * 4.184
            write_stream("""  <Angle class1="%s" class2="%s" class3="%s" angle="%s" k="%s"/>""" % (angle[0], angle[1], angle[2], str(theta), str(k)))
        write_stream(" </HarmonicAngleForce>")
        write_stream(" <PeriodicTorsionForce>")
        processed = set()
        for tor in reversed(self.torsions):
            signature = (fix(tor[0]), fix(tor[1]), fix(tor[2]), fix(tor[3]))
            if signature in processed:
                continue
            if any([c in skipClasses for c in signature]):
                continue
            processed.add(signature)
            tag = "  <Proper class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\"" % signature
            i = 4
            while i < len(tor):
                index = i / 3
                periodicity = int(float(tor[i + 2]))
                phase = float(tor[i + 1]) * math.pi / 180.0
                k = tor[i] * 4.184
                tag += " periodicity%d=\"%d\" phase%d=\"%s\" k%d=\"%s\"" % (index, periodicity, index, str(phase), index, str(k))
                i += 3
            tag += "/>"
            write_stream(tag)
        processed = set()
        for tor in reversed(self.impropers):
            signature = (fix(tor[2]), fix(tor[0]), fix(tor[1]), fix(tor[3]))
            if signature in processed:
                continue
            if any([c in skipClasses for c in signature]):
                continue
            processed.add(signature)
            tag = "  <Improper class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\"" % signature
            i = 4
            while i < len(tor):
                index = i / 3
                periodicity = int(float(tor[i + 2]))
                phase = float(tor[i + 1]) * math.pi / 180.0
                k = float(tor[i]) * 4.184
                tag += " periodicity%d=\"%d\" phase%d=\"%s\" k%d=\"%s\"" % (index, periodicity, index, str(phase), index, str(k))
                i += 3
            tag += "/>"
            write_stream(tag)
        write_stream(" </PeriodicTorsionForce>")
        write_stream(""" <NonbondedForce coulomb14scale="%g" lj14scale="%s">""" % (charge14scale, epsilon14scale))
        sigmaScale = 0.1 * 2.0 / (2.0 ** (1.0 / 6.0))
        for index, type in enumerate(self.types):
            atomClass = type[0]
            q = type[2]
            if atomClass in self.vdwEquivalents:
                atomClass = self.vdwEquivalents[atomClass]
            if atomClass in self.vdw:
                params = [float(x) for x in self.vdw[atomClass]]
                if self.vdwType == 'RE':
                    sigma = params[0] * sigmaScale
                    epsilon = params[1] * 4.184
                else:
                    sigma = (params[0] / params[1]) ** (1.0 / 6.0)
                    epsilon = 4.184 * params[1] * params[1] / (4 * params[0])
            else:
                sigma = 1.0
                epsilon = 0
            if q != 0 or epsilon != 0:
                write_stream("""  <Atom type="%s" charge="%s" sigma="%s" epsilon="%s"/>""" % (self.type_names[index], q, sigma, epsilon))
        write_stream(" </NonbondedForce>")
        write_stream("</ForceField>")
        stream.seek(0)

        return stream

    def parse_filenames(self, filenames):
        """Process a list of filenames according to their filetype suffixes

        Parameters
        ----------
        filenames : list (of strings)
            List of filenames of type (lib, off, dat, or mol2)

        Notes
        -----
        When parameterizing small molecules, the correct order of inputs is:

        $AMBER_LIB_PATH/gaff.dat ligand_name.mol2 ligand_name.frcmod

        """
        for inputfile in filenames:
            if inputfile.endswith('.lib') or inputfile.endswith('.off'):
                self.process_library_file(inputfile)
            elif inputfile.endswith('.dat'):
                self.process_dat_file(inputfile)
            elif inputfile.endswith("mol2"):
                self.process_mol2_file(inputfile)
            else:
                self.process_frc_file(inputfile)

        self.reduce_atomtypes()

    def reduce_atomtypes(self, symmetrize_protons=False):
        """Reduce the list of atom self.types.

        Parameters
        ----------
        symmetrize_protons : bool, default=False
            if True, multiple hydrogens bound to the same heavy atom
            should all use the same type.

        Notes
        -----

        The default behavior of symmetrize_protons differs from the
        original OpenMM version of this script.  For arbitrary small
        molecules, we can not assume symmetric protons.
        """

        removeType = [False] * len(self.types)
        for res in self.residueAtoms:
            if res not in self.residueBonds:
                continue
            atomBonds = [[] for atom in self.residueAtoms[res]]
            for bond in self.residueBonds[res]:
                atomBonds[bond[0]].append(bond[1])
                atomBonds[bond[1]].append(bond[0])
            if symmetrize_protons is True:
                for index, atom in enumerate(self.residueAtoms[res]):
                    hydrogens = [x for x in atomBonds[index] if self.types[self.residueAtoms[res][x][1]][1] == element.hydrogen]
                    for h in hydrogens[1:]:
                        removeType[self.residueAtoms[res][h][1]] = True
                        self.residueAtoms[res][h][1] = self.residueAtoms[res][hydrogens[0]][1]
        newTypes = []
        replaceWithType = [0] * len(self.types)
        for i in range(len(self.types)):
            if not removeType[i]:
                newTypes.append(self.types[i])
            replaceWithType[i] = len(newTypes) - 1
        self.types = newTypes
        for res in self.residueAtoms:
            for atom in self.residueAtoms[res]:
                atom[1] = replaceWithType[atom[1]]

    def set_provenance(self):
        """Set the provenance attribute with information about the current python session."""
        self.provenance = []
        line = """<!-- %s -->\n""" % "Time and parameters of origin:"
        self.provenance.append(line)
        now = datetime.datetime.now()
        line = """<!-- %s -->\n""" % str(now)
        self.provenance.append(line)
        cmd_string = subprocess.list2cmdline(sys.argv[1:])
        cmd_string = cmd_string.replace("-", " ")  # Replace XML specific characters that can break some XML parsers
        cmd_string = cmd_string.replace(">", " ")  #
        cmd_string = cmd_string.replace("<", " ")  #
        line = """<!-- %s -->\n""" % cmd_string
        self.provenance.append(line)
        self.provenance = "".join(self.provenance)
