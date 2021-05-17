#!/usr/bin/env python
import sys
import math
from simtk.openmm.app import Element
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
for elem in Element._elements_by_symbol.values():
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

# Manually create the hydrogen element
hydrogen = Element(1, "hydrogen", "H", 1.007947*unit.daltons)


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
            e = Element.getBySymbol(element_symbol)
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
            # Use to detect blank lines
            line_length = len(line.strip())
            if block == 0:     # Title
                block += 1
            elif block == 1:   # Mass
                if line_length == 0:
                    block += 1
                else:
                    params = self._parse_dat_atom_symbols_and_masses(line)
                    self.masses[params['kndsym']] = float(params['amass'])
            elif block == 2:   # Hydrophilic atoms
                block += 1
            elif block == 3:   # Bonds
                if line_length == 0:
                    block += 1
                else:
                    params = self._parse_dat_bond_length_parameters(line)
                    self.bonds.append([params['ibt'], params['jbt'], params['rk'], params['req']])
            elif block == 4:   # Angles
                if line_length == 0:
                    block += 1
                else:
                    params = self._parse_dat_bond_angle_parameters(line)
                    self.angles.append([params['itt'], params['jtt'], params['ktt'], params['tk'], params['teq']])
            elif block == 5:   # Torsions
                if line_length == 0:
                    block += 1
                else:
                    params = self._parse_dat_dihedral_parameters(line)
                    # Periodicity parameter pn is an int stored as a float,
                    # and a negative sign indicates additional dihedral terms are added on the next line
                    pn = int(float(params['pn']))
                    pk_over_idivf = float(params['pk']) / float(params['idivf'])
                    if continueTorsion:
                        self.torsions[-1] += [pk_over_idivf, params['phase'], abs(pn)]
                    else:
                        self.torsions.append([params['ipt'], params['jpt'], params['kpt'], params['lpt'], pk_over_idivf, params['phase'], abs(pn)])
                    continueTorsion = (pn < 0)
            elif block == 6:   # Improper torsions
                if line_length == 0:
                    block += 1
                else:
                    params = self._parse_dat_improper_dihedral_parameters(line)
                    self.impropers.append((params['ipt'], params['jpt'], params['kpt'], params['lpt'],params['pk'], params['phase'], params['pn']))
            elif block == 7:   # 10-12 hbond potential
                if line_length == 0:
                    block += 1
            elif block == 8:   # VDW equivalents
                if line_length == 0:
                    block += 1
                else:
                    symbols = self._parse_dat_6_12_equivalence_symbols(line)
                    for atom in symbols['ieqv']:
                        self.vdwEquivalents[atom] = symbols['iorg']
            elif block == 9:   # VDW type
                block += 1
                spec = self._parse_dat_6_12_potential_kind(line)
                self.vdwType = spec['kindnb'].upper()
                if self.vdwType not in ['RE', 'AC']:
                    raise ValueError('Nonbonded type (KINDNB) must be RE or AC')
            elif block == 10:   # VDW parameters
                if line_length == 0:
                    block += 1
                else:
                    params = self._parse_dat_6_12_nb_parameters(line, spec['kindnb'])
                    if self.vdwType == "RE":
                        self.vdw[params['ltynb']] = (params['r'], params['edep'])
                    elif self.vdwType == "AC":
                        self.vdw[params['ltynb']] = (params['a'], params['c'])

    @staticmethod
    def _parse_dat_atom_symbols_and_masses(line):
        """
        Parse a line in a parm.dat file using atom symbol and mass specification.

        Parameters
        ----------
        line : str
            Single line string containing atom symbol and mass parameters in a parm.dat file.

        Returns
        -------
        dict containing
        kndsym : str
        amass : str
        atpol : str
        line :str

        Notes
        -----
        Original format specification http://ambermd.org/formats.html#parm.dat

        - 2 -      ***** INPUT FOR ATOM SYMBOLS AND MASSES *****

                    KNDSYM , AMASS, ATPOL

                        FORMAT(A2,2X,F10.2x,f10.2)

         KNDSYM     The unique atom symbol used in the system.

         AMASS      Atomic mass of the center having the symbol "KNDSYM".

         ATPOL      The atomic polarizability for each atom (in A**3)
                    This is the type of polarizability used in sander
                    and gibbs. No parameters are supplied for this since
                    the feature is still in development (Amber 4.1).

              NOTE: All the unique atomic symbols and their masses must
                    be read.  The input is terminated by a blank card.

        Examples
        --------
        c3 12.01         0.878               Sp3 C
        ca 12.01         0.360               Sp2 C in pure aromatic systems

        """
        kndsym = line[0:2].strip()
        amass = line[4:14].strip()

        # prevent potential IndexError from line being too short
        try:
            apol = line[14:24].split()[0].strip()
        except IndexError:
            apol = line[14:-1].split()[0].strip()
        return locals()



    @staticmethod
    def _parse_dat_bond_length_parameters(line):
        """
        Parse a line in a parm.dat file using bond length format specification.

        Parameters
        ----------
        line : str
            Single line string containing bond length parameters in a parm.dat file.

        Returns
        -------
        dict containing
        ibt : str
        jbt : str
        rk : str
        req : str
        line : str

        Notes
        -----
        Original format specification http://ambermd.org/formats.html#parm.dat

        - 4 -      ***** INPUT FOR BOND LENGTH PARAMETERS *****

                    IBT , JBT , RK , REQ

                        FORMAT(A2,1X,A2,2F10.2)

         IBT,JBT    Atom symbols for the two bonded atoms.

         RK         The harmonic force constant for the bond "IBT"-"JBT".
                    The unit is kcal/mol/(A**2).

         REQ        The equilibrium bond length for the above bond in Angstroms

                    The input is terminated by a blank card.

        Examples
        --------

        n -os  395.0    1.4103       SOURCE4      30	0.0112
        no-s4  143.0    1.9960       SOURCE3       3	0.0313
        no-s6  149.6    1.9760       SOURCE3       3	0.0520
        """

        ibt = line[0:2].strip()
        jbt = line[3:5].strip()
        rk = line[5:15].strip()
        # prevent potential IndexError from line being too short
        try:
            req = line[15:25].split()[0].strip()
        except IndexError:
            req = line[15:-1].split()[0].strip()

        return locals()

    @staticmethod
    def _parse_dat_bond_angle_parameters(line):
        """
        Parse a line in a parm.dat file using bond angle format specification.
        Parameters
        ----------
        line : str
            Single line string containing bond angle parameters in a parm.dat file.

        Returns
        -------
        dict containing
        itt : str
        jtt : str
        ktt : str
        tk : str
        teq : str
        line : str

        Notes
        -----
        Original format specification http://ambermd.org/formats.html#parm.dat

        - 5 -      ***** INPUT FOR BOND ANGLE PARAMETERS *****


                    ITT , JTT , KTT , TK , TEQ

                        FORMAT(A2,1X,A2,1X,A2,2F10.2)

         ITT,...    The atom symbols for the atoms making an angle.

         TK         The harmonic force constants for the angle "ITT"-"JTT"-
                    "KTT" in units of kcal/mol/(rad**2) (radians are the
                    traditional unit for angle parameters in force fields).

         TEQ        The equilibrium bond angle for the above angle in degrees.

                    The input is terminated by a blank card.

        Examples
        --------
        n3-c3-n3   69.61      109.59	SOURCE4           27	1.8125
        n3-c3-nc   68.79      113.29	SOURCE3            1	0.0000
        n3-c3-nd   68.79      113.29	SOURCE3            1	same_as_n3-c3-nc
        c1-sh-hs   48.23       95.99	calculated_based_on_C#C-SH       0
        """
        itt = line[0:2].strip()
        jtt = line[3:5].strip()
        ktt = line[6:8].strip()
        tk = line[8:18].strip()

        # prevent potential IndexError from line being too short
        try:
            teq = line[18:28].split()[0].strip()
        except IndexError:
            teq = line[18:-1].split()[0].strip()
        return locals()

    @staticmethod
    def _parse_dat_dihedral_parameters(line):
        """
        Parse a line in a parm.dat file using dihedral format specification.
        Parameters
        ----------
        line : str
            Single line string containing dihedral parameters in a parm.dat file.

        Returns
        -------
        dict containing
        ipt : str
        jpt : str
        kpt : str
        lpt : str
        idivf : str
        pk : str
        phase : str
        pn : str
        line : str

        Notes
        -----
        Original format specification http://ambermd.org/formats.html#parm.dat

        - 6 -      ***** INPUT FOR DIHEDRAL PARAMETERS *****

                    IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN

                        FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)

         IPT, ...   The atom symbols for the atoms forming a dihedral
                    angle.  If IPT .eq. 'X ' .and. LPT .eq. 'X ' then
                    any dihedrals in the system involving the atoms "JPT" and
                    and "KPT" are assigned the same parameters.  This is
                    called the general dihedral type and is of the form
                    "X "-"JPT"-"KPT"-"X ".

         IDIVF      The factor by which the torsional barrier is divided.
                    Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
                    details. Basically, the actual torsional potential is

                           (PK/IDIVF) * (1 + cos(PN*phi - PHASE))

         PK         The barrier height divided by a factor of 2.

         PHASE      The phase shift angle in the torsional function.

                    The unit is degrees.

         PN         The periodicity of the torsional barrier.
                    NOTE: If PN .lt. 0.0 then the torsional potential
                          is assumed to have more than one term, and the
                          values of the rest of the terms are read from the
                          next cards until a positive PN is encountered.  The
                          negative value of pn is used only for identifying
                          the existence of the next term and only the
                          absolute value of PN is kept.

            The input is terminated by a blank card.

        Examples
        --------
        X -c -cy-X    6    0.000       180.000           2.000      JCC, 7, (1986), 230
        X -c -ca-X    4    4.000       180.000           2.000      optimized by Junmei Wang, Jan-2013
        X -c -cc-X    4   11.500       180.000           2.000      statistic value


        """

        ipt = line[0:2].strip()
        jpt = line[3:5].strip()
        kpt = line[6:8].strip()
        lpt = line[9:11].strip()
        idivf = line[11:15].strip()
        pk = line[15:30].strip()
        phase = line[30:45].strip()
        # prevent potential IndexError from line being too short
        try:
            pn = line[45:60].split()[0].strip()
        except IndexError:
            pn = line[45:-1].split()[0].strip()

        return locals()

    @staticmethod
    def _parse_dat_improper_dihedral_parameters(line):
        """
        Parse a line in a parm.dat file using improper dihedral format specification.
        Parameters
        ----------
        line : str
            Single line string containing dihedral parameters in a parm.dat file.

        Returns
        -------
        dict containing
        ipt : str
        jpt : str
        kpt : str
        lpt : str
        idivf : str
        pk : str
        phase : str
        pn : str
        line : str

        Notes
        -----
        Original format specification http://ambermd.org/formats.html#parm.dat

         - 7 -      ***** INPUT FOR IMPROPER DIHEDRAL PARAMETERS *****

                    IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN

                        FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)

                    The input is the same as in for the dihedrals except that
                    the torsional barrier height is NOT divided by the factor
                    idivf.  The improper torsions are defined between any four
                    atoms not bonded (in a successive fashion) with each other
                    as in the case of "regular" or "proper" dihedrals.  Improper
                    dihedrals are used to keep certain groups planar and to
                    prevent the racemization of certain centers in the united
                    atom model.  Consult the above reference for details.

             Important note: all general type improper dihedrals
                             (e.g. x -x -ct-hc) should appear before all
                             specifics (ct-ct-ct-hc) in the parm list.
                             Otherwise the generals will override the
                             specific with no warning.

             The input is terminated by a blank card.

        Examples
        --------
        X -o -c -o          1.1          180.          2.           JCC,7,(1986),230
        X -X -c -o          10.5         180.          2.           JCC,7,(1986),230


        """

        ipt = line[0:2].strip()
        jpt = line[3:5].strip()
        kpt = line[6:8].strip()
        lpt = line[9:11].strip()
        idivf = line[11:15].strip()
        pk = line[15:30].strip()
        phase = line[30:45].strip()
        # prevent potential IndexError from line being too short
        try:
            pn = line[45:60].split()[0].strip()
        except IndexError:
            pn = line[45:-1].split()[0].strip()

        return locals()

    @staticmethod
    def _parse_dat_6_12_equivalence_symbols(line):
        """
        Parse a line in a parm.dat file using equivalencing symbols for non-bonded 6-12 potential specification.

        Parameters
        ----------
        line : str
            Single line string containing equivalencing symbols for 6-12 non-bonded parameters in a parm.dat file.

        Returns
        -------
        dict containing
        iorg : str
        ieqv : list of str
        line : str

        Notes
        -----
        Original format specification http://ambermd.org/formats.html#parm.dat

        - 9 -      ***** INPUT FOR EQUIVALENCING ATOM SYMBOLS FOR
                          THE NON-BONDED 6-12 POTENTIAL PARAMETERS *****

                          IORG , IEQV(I) , I = 1 , 19

                              FORMAT(20(A2,2X))

         IORG        The atom symbols to which other atom symbols are to be
                     equivalenced in generating the 6-12 potential parameters.

         IEQV(I)     The atoms symbols which are to be equivalenced to the
                     atom symbol "IORG".  If more than 19 atom symbols have
                     to be equivalenced to a given atom symbol they can be
                     included as extra cards.

                     It is advisable not to equivalence any hydrogen bond
                     atom type atoms with any other atom types.

          NOTE: The input is terminated by a blank card.

        """
        ieqv = list()
        iorg = line[0:2].strip()
        # continue adding names till line runs out,
        # or reaches 19 which is the maximum according to format
        try:
            for n in range(1,20):
                ieqv.append(line[4*n:4*n+2].split()[0].strip())
        except IndexError:
            pass

        return dict(ieqv=ieqv, iorg=iorg, line=line)

    @staticmethod
    def _parse_dat_6_12_potential_kind(line):
        """
        Parse a line in a parm.dat file using input for non-bonded 6-12 potential specification.

        Parameters
        ----------
        line : str
            Single line string containing the input format for 6-12 non-bonded parameters in a parm.dat file.

        Returns
        -------
        dict containing
        label : str
        kindb : str
        line : str

        Notes
        -----
        Original format specification http://ambermd.org/formats.html#parm.dat

        - 10 -      ***** INPUT FOR THE 6-12 POTENTIAL PARAMETERS *****

                     LABEL , KINDNB

                         FORMAT(A4,6X,A2)

         LABEL       The name of the non-bonded input parameter to be
                     used.  It has to be matched with "NAMNB" read through
                     unit 5.  The program searches the file to load the
                     the required non-bonded parameters.  If that name is
                     not found the run will be terminated.

         KINDNB      Flag for the type of 6-12 parameters.

          'SK'       Slater-Kirkwood parameters are input.
                     see "caution" below.

          'RE'       van der Waals radius and the potential well depth
                     parameters are read.

          'AC'       The 6-12 potential coefficients are read.

             NOTE: All the non equivalenced atoms' parameters have to
                   be given.

           The input is terminated when label .eq. 'END'

        Examples
        --------
        MOD4      RE

        """

        label = line[0:4].strip()
        kindnb = line[10:12].strip()

        if kindnb not in ["SK", "RE", "AC"]:
            raise ValueError("Unsupported 6-12 potential format {kindb}".format(**locals()))

        return locals()

    @staticmethod
    def _parse_dat_6_12_nb_parameters(line, kindnb):
        """
        Parse a line in a parm.dat file using RE format for 6-12 potential specification.

        Parameters
        ----------
        line : str
            Single line string containing equivalencing symbols for 6-12 non-bonded parameters in a parm.dat file.

        kindnb : str
            The kind of format for the nonbonded parameter line ("SK", "RE", or "AC")

        Returns
        -------
        dict containing
        lytnb : str
        line : str
        kindnb : str

        and for SK

        pol : str
        xneff : str
        rmin :str

        or for RE

        r : str
        edep : str

        or for AC
        a : str
        c : str


        Notes
        -----

        This code assumes the format is   FORMAT(2X,A2,6X,2F10.6) for 10B and 10C

        Original format specification http://ambermd.org/formats.html#parm.dat
         - 10A -     ***** ONLY IF KINDNB .EQ. 'SK' *****

                     LTYNB , POL , XNEFF , RMIN

                         FORMAT(2X,A2,6X,3F10.6)

         LTYNB       Atom symbol.

         POL         Atomic polarizability for the atom centers having the
                     the above symbol.

         XNEFF       Effective number of electrons on the atom centers having
                     the above symbol.

         RMIN        van der Waals radius of the atom center having the above
                     symbol.


        - 10B -     ***** ONLY IF KINDNB .EQ. 'RE' *****

                     LTYNB , R , EDEP

         LTYNB       Atom symbol.

         R           The van der Waals radius of the atoms having the symbol
                     "LTYNB"  (Angstoms)

         EDEP        The 6-12 potential well depth. (kcal/mol)

        ------------------------------------------------------------------------

         - 10C -     ***** ONLY IF KINDNB .EQ. 'AC' *****

                     LTYNB , A , C

         LTYNB       Atom symbol.

         A           The coefficient of the 12th power term (A/r**12).

         C           The coefficient of the 6th power term (-C/r**6).

        Examples
        --------
          c1          1.9080  0.2100             cp C DLM 11/2007 well depth from OPLS replacing 0.0860
          c2          1.9080  0.0860             sp2 atom in the middle of C=CD-CD=C

        """
        ltynb = line[2:4].strip()
        if kindnb.upper() == "SK":
            pol = line[10:20].strip()
            xneff = line[20:30].strip()

            # prevent IndexError from line being too short
            try:
                rmin = line[30:40].split()[0].strip()
            except IndexError:
                rmin = line[30:-1].split()[0].strip()
        elif kindnb.upper() == "RE":
            r = line[10:20].strip()
            # prevent IndexError from line being too short
            try:
                edep = line[20:30].split()[0].strip()
            except IndexError:
                edep = line[20:-1].split()[0].strip()
        elif kindnb.upper() == "AC":
            a = line[10:20].strip()
            # prevent IndexError from line being too short
            try:
                c = line[20:30].split()[0].strip()
            except IndexError:
                c = line[20:-1].split()[0].strip()
        else:
            raise ValueError("Unsupported NB format {nbformat}".format(**locals()))

        return locals()

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
                    hydrogens = [x for x in atomBonds[index] if self.types[self.residueAtoms[res][x][1]][1] == hydrogen]
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
