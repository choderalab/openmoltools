import mdtraj
import numpy as np
import pandas as pd
import cStringIO
import simtk.openmm.app.element
import simtk.unit as unit
import itertools

GAFF_NONBONDED_PATH = "/home/kyleb/src/kyleabeauchamp/gaff2xml/xml/gaff_nonbonded.dat"
GAFF_ELEMENTS_PATH = "/home/kyleb/src/kyleabeauchamp/gaff2xml/xml/gaff_elements.csv"

gaff_elements = {'br': 'Br',
 'c': 'C',
 'c1': 'C',
 'c2': 'C',
 'c3': 'C',
 'ca': 'C',
 'cc': 'C',
 'cd': 'C',
 'ce': 'C',
 'cf': 'C',
 'cg': 'C',
 'ch': 'C',
 'cl': 'Cl',
 'cp': 'C',
 'cq': 'C',
 'cu': 'C',
 'cv': 'C',
 'cx': 'C',
 'cy': 'C',
 'cz': 'C',
 'f': 'F',
 'h1': 'H',
 'h2': 'H',
 'h3': 'H',
 'h4': 'H',
 'h5': 'H',
 'ha': 'H',
 'hc': 'H',
 'hn': 'H',
 'ho': 'H',
 'hp': 'H',
 'hs': 'H',
 'hw': 'H',
 'hx': 'H',
 'i': 'I',
 'n': 'N',
 'n1': 'N',
 'n2': 'N',
 'n3': 'N',
 'n4': 'N',
 'na': 'N',
 'nb': 'N',
 'nc': 'N',
 'nd': 'N',
 'ne': 'N',
 'nf': 'N',
 'nh': 'N',
 'no': 'N',
 'o': 'O',
 'oh': 'O',
 'os': 'O',
 'ow': 'O',
 'p2': 'P',
 'p3': 'P',
 'p4': 'P',
 'p5': 'P',
 'pb': 'P',
 'px': 'P',
 'py': 'P',
 's': 'S',
 's2': 'S',
 's4': 'S',
 's6': 'S',
 'sh': 'S',
 'ss': 'S',
 'sx': 'S',
 'sy': 'S'}

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

