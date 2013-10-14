import mdtraj
import numpy as np
import pandas as pd
import cStringIO
import pybel
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
        xyz_mm = traj.xyz[0].tolist()

        return top_mm, xyz_mm
