import pandas as pd
import cStringIO
import pybel
import time

nonbonded = load_gaff_nonbonded()

def correct_types(pybel_mol2, pandas_mol2):
    for k, atom_type in enumerate(pandas_mol2.type):
        pybel_mol2.OBMol.GetAtom(k + 1).SetType(atom_type)
        
        #pybel_mol2.atoms[k].OBAtom.SetType(atom_type)

def read_mol2_pandas(filename):
    handle = open(filename, 'r')
    csv_data = cStringIO.StringIO()
    appending = False

    for line in handle.readlines():
    
        if appending == True and """@<TRIPOS>""" in line:
            appending = False
            break

        if appending == True:
            csv_data.write(line)

        if """@<TRIPOS>ATOM""" in line:
            appending = True
            
    csv_data.reset()
    pandas_mol2 = pd.read_csv(csv_data, delim_whitespace=True, names=["id", "name", "x", "y" ,"z", "type", "code", "resname", "charge"], index_col=0, header=None, usecols=range(1, 10))
    return pandas_mol2

def load_antechamber_mol2(filename):
    pandas_mol2 = read_mol2_pandas(filename)
    pybel_mol2 = pybel.readfile("mol2", filename).next()
    correct_types(pybel_mol2, pandas_mol2)
    return pybel_mol2
        
def load_gaff_nonbonded():
    nonbonded = pd.read_fwf("./gaff_nonbonded.dat", widths=(4,16, 8), names=["type", "sigma", "epsilon"], index_col=0)
    return nonbonded


def write_gaff_xml(mol2_filename):
    pybel_mol2 = load_antechamber_mol2(mol2_filename)
    
"""
  <Atom type="1960" charge="1.0" sigma="0.526699322165" epsilon="0.00071128"/>
 </NonbondedForce>
</ForceField>
"""
