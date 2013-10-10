import pandas as pd
import cStringIO
import pybel
import simtk.openmm.app.element
import simtk.unit as unit

GAFF_NONBONDED_PATH = "/home/kyleb/src/kyleabeauchamp/gaff2xml/xml/gaff_nonbonded.dat"
GAFF_ELEMENTS_PATH = "/home/kyleb/src/kyleabeauchamp/gaff2xml/xml/gaff_elements.csv"

def correct_types(pybel_mol2, pandas_mol2):
    """Fix mol2 elements and types that are broken by OpenBabel due to non-standard GAFF mol2 files.
    
    Notes
    -----
    
    It appears necessary to do both GetType() and SetType() for each atom
    in order to "update" the internal state of the OBAtom objects.  Likely
    indicates bugs in pybel / python-openbabel.
    """
    for k, atom_type in enumerate(pandas_mol2.atype):
        old_type = pybel_mol2.atoms[k].OBAtom.GetType()  # Do not delete this line
        pybel_mol2.atoms[k].OBAtom.SetType(atom_type)
        new_type = pybel_mol2.atoms[k].OBAtom.GetType()  # Do not delete this line
        print("Changing type of atom %d from %s to %s." % (k, old_type, new_type))

    for k, atom_type in enumerate(pandas_mol2.atype):
        new_type = pybel_mol2.atoms[k].OBAtom.GetType()
        element = gaff_elements[new_type]
        omm_element = simtk.openmm.app.element.get_by_symbol(element)
        atomic_number = omm_element.atomic_number
        pybel_mol2.atoms[k].OBAtom.SetAtomicNum(atomic_number)

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
    pandas_mol2 = pd.read_csv(csv_data, delim_whitespace=True, names=["id", "name", "x", "y" ,"z", "atype", "code", "resname", "charge"], index_col=0, header=None, usecols=range(1, 10))
    return pandas_mol2

def load_antechamber_mol2(filename):
    pandas_mol2 = read_mol2_pandas(filename)
    pybel_mol2 = pybel.readfile("mol2", filename).next()
    correct_types(pybel_mol2, pandas_mol2)
    return pybel_mol2, pandas_mol2
        
def load_gaff_nonbonded():
    nonbonded = pd.read_fwf(GAFF_NONBONDED_PATH, widths=(4,16, 8), names=["atype", "sigma", "epsilon"], index_col=0)
    return nonbonded


def generate_gaff_xml(pandas_mol2, pybel_mol2):
    residue_name = pandas_mol2.resname[1]  # To Do: Add check for consistency
    xml_text = cStringIO.StringIO()
    xml_text.write("<ForceField>\n")

    xml_text.write("<AtomTypes>\n")
    for (i, name, x, y, z, atype, code, resname, charge) in pandas_mol2.itertuples():
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
    for (i, name, x, y, z, atype, code, resname, charge) in pandas_mol2.itertuples():
        sigma, epsilon = lookup_vdw(atype)
        full_name = residue_name + "_" + name
        element = gaff_elements[atype]
        omm_element = simtk.openmm.app.element.get_by_symbol(element)
        mass = omm_element.mass / unit.daltons        
        line = """   <Atom name="%s" type="%s"/>\n""" % (name, full_name)
        xml_text.write(line)
    
    for k in range(pybel_mol2.OBMol.NumBonds()):
        b = pybel_mol2.OBMol.GetBond(k)
        i = b.GetBeginAtomIdx() - 1  # Subtract 1 for zero based indexing in OpenMM???  
        j = b.GetEndAtomIdx() - 1  # Subtract 1 for zero based indexing in OpenMM???  
        xml_text.write("""<Bond from="%d" to="%d"/>\n""" % (i, j))
        
    xml_text.write("</Residue>\n")
    xml_text.write("</Residues>\n")


    xml_text.write("""<NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">\n""")
    for (i, name, x, y, z, atype, code, resname, charge) in pandas_mol2.itertuples():
        sigma, epsilon = lookup_vdw(atype)
        full_name = residue_name + "_" + name    
        line = """<Atom type="%s" charge="%f" sigma="%f" epsilon="%f"/>\n""" % (full_name, charge, sigma, epsilon)
        xml_text.write(line)
    xml_text.write("</NonbondedForce>\n")
    
    xml_text.write("</ForceField>\n")
    xml_text.reset()
    return xml_text

    
def lookup_vdw(atype):
    return nonbonded.ix[atype]

nonbonded = load_gaff_nonbonded()
gaff_elements = pd.read_csv(GAFF_ELEMENTS_PATH, index_col=0).element
