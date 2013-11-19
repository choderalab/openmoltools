import mdtraj
import numpy as np
import pandas as pd
import cStringIO
import simtk.openmm.app.element
import simtk.unit as unit
import itertools

import logging
logger = logging.getLogger(__name__)

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
        xyzlist /= 10.0  # Convert from angstrom to nanometer
        traj = mdtraj.Trajectory(xyzlist, top)
        return traj

    def to_openmm(self):
        traj = self.to_mdtraj()
        top_mm = traj.top.to_openmm()
        xyz_mm = traj.xyz[0].tolist()

        return top_mm, xyz_mm


