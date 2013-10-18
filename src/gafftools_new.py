import mdtraj
import numpy as np
import pandas as pd
import cStringIO
import simtk.openmm.app.element
import simtk.unit as unit
import itertools

GAFF_NONBONDED_PATH = "/home/kyleb/src/kyleabeauchamp/gaff2xml/xml/gaff_nonbonded.dat"
GAFF_ELEMENTS_PATH = "/home/kyleb/src/kyleabeauchamp/gaff2xml/xml/gaff_elements.csv"

def parse_mol2_sections(x):    
    if x.startswith('@<TRIPOS>'):
        parse_mol2_sections.key = x
    return parse_mol2_sections.key

def mol2_to_dataframes(filename):
    """Convert a GAFF mol2 file to a pair of pandas dataframes."""
    with open(filename) as f:
        data = dict((key, list(grp)) for key, grp in itertools.groupby(f, parse_mol2_sections))

    csv = cStringIO.StringIO()
    csv.writelines(data["@<TRIPOS>BOND\n"][1:])
    csv.reset()
    bonds_frame = pd.read_table(csv, delim_whitespace=True, names=["bond_id", "id0","id1", "bond_type"], index_col=0, header=None, sep="\s*",)

    csv = cStringIO.StringIO()
    csv.writelines(data["@<TRIPOS>ATOM\n"][1:])
    csv.reset()
    atoms_frame = pd.read_csv(csv, delim_whitespace=True, names=["serial", "name", "x", "y" ,"z", "atype", "code", "resName", "charge"], header=None, usecols=range(1, 10))

    return atoms_frame, bonds_frame

def load_gaff_nonbonded():
    nonbonded = pd.read_fwf(GAFF_NONBONDED_PATH, widths=(4,16, 8), names=["atype", "sigma", "epsilon"], index_col=0)
    return nonbonded

def lookup_vdw(atype):
    return nonbonded.ix[atype]

nonbonded = load_gaff_nonbonded()
gaff_elements = pd.read_csv(GAFF_ELEMENTS_PATH, index_col=0).element

class Mol2Parser(object):
    def __init__(self, filename):
        self.atoms, self.bonds = mol2_to_dataframes(filename)

    def to_mdtraj(self):
        atoms, bonds = self.atoms, self.bonds
        atoms_mdtraj = atoms[["name", "resName"]]
        atoms_mdtraj["serial"] = atoms.index
        atoms_mdtraj["element"] = atoms.atype.map(gaff_elements)
        atoms_mdtraj["resSeq"] = np.ones(len(atoms))
        atoms_mdtraj["chainID"] = np.ones(len(atoms))

        bonds_mdtraj = bonds[["id0", "id1"]].values
        offset = bonds_mdtraj.min()
        bonds_mdtraj -= offset

        top = mdtraj.Topology.from_dataframe(atoms_mdtraj, bonds_mdtraj)
        xyzlist = np.array([atoms[["x","y","z"]].values])
        traj = mdtraj.Trajectory(xyzlist, top)
        return traj

    def to_openmm(self):
        traj = self.to_mdtraj()
        top_mm = traj.top.to_openmm()
        xyz_mm = (traj.xyz[0] / 10.0).tolist()  # Convert from angstrom to nanometer

        return top_mm, xyz_mm


def generate_gaff_xml(atoms, bonds):
    residue_name = atoms.resName[1]  # To Do: Add check for consistency
    xml_text = cStringIO.StringIO()
    xml_text.write("<ForceField>\n")

    xml_text.write("<AtomTypes>\n")
    for (i, name, x, y, z, atype, code, resname, charge) in atoms.itertuples(False):
        sigma, epsilon = lookup_vdw(atype)
        full_name = residue_name + "_" + name
        element = gaff_elements[atype]
        omm_element = simtk.openmm.app.element.get_by_symbol(element)
        mass = omm_element.mass / unit.daltons
        
        line = """<Type name="%s" class="%s" element="%s" mass="%f"/>\n""" % (full_name, atype, element, mass)
        xml_text.write(line)
    xml_text.write("</AtomTypes>\n")

    xml_text.write("<Residues>\n")
    xml_text.write("""<Residue name="%s">\n""" % residue_name)
    for (i, name, x, y, z, atype, code, resname, charge) in atoms.itertuples(False):
        sigma, epsilon = lookup_vdw(atype)
        full_name = residue_name + "_" + name
        element = gaff_elements[atype]
        omm_element = simtk.openmm.app.element.get_by_symbol(element)
        mass = omm_element.mass / unit.daltons        
        line = """   <Atom name="%s" type="%s"/>\n""" % (name, full_name)
        xml_text.write(line)
    
    for (id0, id1, bond_type) in bonds.itertuples(False):
        i = id0 - 1  # Subtract 1 for zero based indexing in OpenMM???  
        j = id1 - 1  # Subtract 1 for zero based indexing in OpenMM???  
        xml_text.write("""<Bond from="%d" to="%d"/>\n""" % (i, j))
        
    xml_text.write("</Residue>\n")
    xml_text.write("</Residues>\n")


    xml_text.write("""<NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">\n""")
    for (i, name, x, y, z, atype, code, resname, charge) in atoms.itertuples(False):
        sigma, epsilon = lookup_vdw(atype)
        full_name = residue_name + "_" + name    
        line = """<Atom type="%s" charge="%f" sigma="%f" epsilon="%f"/>\n""" % (full_name, charge, sigma, epsilon)
        xml_text.write(line)
    xml_text.write("</NonbondedForce>\n")
    
    xml_text.write("</ForceField>\n")
    xml_text.reset()
    return xml_text

_FRCMOD_HEADERS = ["MASS", "BOND", "ANGL", "DIHE", "IMPR", "NONB"]


def parse_frcmod_sections(x):
    if x[0:4] in _FRCMOD_HEADERS:
        parse_frcmod_sections.key = x[0:4]
    return parse_frcmod_sections.key

parse_frcmod_sections.key = None

parsers = {
"BOND": lambda stream: pd.read_fwf(stream, colspecs=((0,2),(3,5), (6, 8), (11, 17), (22, 29)), header=None, names=["type0", "type1", "k", "r0"]), 
"ANGL": lambda stream: pd.read_fwf(stream, colspecs=((0,2),(3,5), (6, 8), (11, 17), (22, 29)), header=None, names=["type0", "type1", "type2", "k", "theta0"]),
"DIHE": lambda stream: pd.read_fwf(stream, colspecs=((0,2),(3,5), (6, 8), (9, 11), (18, 24), (32, 39), (46, 51)), names=["type0", "type1", "type2", "type3", "k", "phase", "periodicity"], header=None),
"IMPR": lambda stream: pd.read_fwf(stream, colspecs=((0,2),(3,5), (6, 8), (9, 11), (18, 24), (32, 39), (46, 51)), names=["type0", "type1", "type2", "type3", "k", "phase", "periodicity"], header=None)
}

def bonds_to_xml(frame, xml_text):
    pass

def angles_to_xml(frame, xml_text):
    xml_text.write("<Angle/>\n")
    for (type0, type1, type2, k, theta) in frame.itertuples(False):
          line = """<Angle class1="%s" class2="%s" class3="%s" angle="%f" k="%f"/>\n""" % (type0, type1, type2, theta * np.pi / 180., k * 4.184)  # Radians and KJ / mol
          xml_text.write(line)
    xml_text.write("</Angle/>\n")
    
def dihedrals_to_xml(frame, xml_text):
    xml_text.write("<<PeriodicTorsionForce>/>\n")
    processed = []
    for (type0, type1, type2, type3, k, phase, periodicity) in frame.itertuples(False):
        signature = (type0, type1, type2, type3)
        if signature in processed:
            continue
        processed.add(signature)

    for (type0, type1, type2, type3, k, phase, periodicity) in frame.itertuples(False):
        signature = (type0, type1, type2, type3)
        if signature in processed:
            continue
        processed.add(signature)

    tag = "  <Proper class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\"" % signature
    i = 4
    while i < len(tor):
        index = i/3
        periodicity = int(float(tor[i+2]))
        phase = float(tor[i+1])*math.pi/180.0
        k = tor[i]*4.184
        tag += " periodicity%d=\"%d\" phase%d=\"%s\" k%d=\"%s\"" % (index, periodicity, index, str(phase), index, str(k))
        i += 3
    tag += "/>"

        #line = """<Proper class1="%s" class2="%s" class3="%s" class4="%s" angle="%f" k="%f"/>\n""" % (type0, type1, type2, type3,)
        #xml_text.write(line)
    xml_text.write("</<PeriodicTorsionForce>/>\n")

def impropers_to_xml(frame, xml_text):
    xml_text.write("<<PeriodicTorsionForce>/>\n")
    for (type0, type1, type2, type3, k, phase, periodicity) in frame.itertuples(False):
        pass
        #line = """<Improper class1="%s" class2="%s" class3="%s" class4="%s" angle="%f" k="%f"/>\n""" % (type0, type1, type2, type3, theta, k)
        #xml_text.write(line)
    xml_text.write("</<PeriodicTorsionForce>/>\n")

xml_mungers = {
"BOND": bonds_to_xml,
"ANGL": angles_to_xml,
"DIHE": dihedrals_to_xml,
"IMPR": impropers_to_xml,
}

def frcmod_to_xml(filename):
    with open(filename) as f:
        data = dict((key, list(grp)) for key, grp in itertools.groupby(f, parse_frcmod_sections))

    section_frames = {}
    for header, lines in data.iteritems():
        if header in parsers.keys() and len(lines) > 2:
            print(header)
            stream = cStringIO.StringIO()
            stream.writelines(lines[1:-1])  # Remove first (header) and last (blank) lines
            stream.reset()
            frame = parsers[header](stream)
            section_frames[header] = frame
            print(frame)

    xml_text = cStringIO.StringIO()
    xml_text.write("<ForceField>\n")

    for header, frame in section_frames.iteritems():
        print(header)
        print(frame)
        xml_mungers[header](frame, xml_text)

    xml_text.write("</ForceField>\n")
    return xml_text
