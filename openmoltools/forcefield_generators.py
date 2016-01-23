#!/usr/bin/env python
"""
OpenMM ForceField residue template generators.

"""

from __future__ import absolute_import

import numpy as np
import os, os.path
from simtk.openmm.app import ForceField
from openmoltools.amber import run_antechamber
from openmoltools.openeye import get_charges
from simtk.openmm.app.element import Element
import parmed

def PerceiveBondOrdersExplicitHydrogens(mol):
    """
    Perceive bond order from topology and element identities rather than geometry.

    Parameters
    ----------
    mol : OEChem.OEMol
        Molecule whose bonds are to be perceived.
        All bonds should be single-bonds.

    """
    from openeye import oechem
    # Expand a tree of potential bond order assignments
    from openeye.oechem import OEGetAtomicNum
    maximum_atomic_valences = {
        OEGetAtomicNum('H') : 1,
        OEGetAtomicNum('C') : 4,
        OEGetAtomicNum('N') : 3,
        OEGetAtomicNum('O') : 2,
        OEGetAtomicNum('P') : 5,
        OEGetAtomicNum('S') : 5,
        OEGetAtomicNum('F') : 1,
        OEGetAtomicNum('Cl') : 1,
        OEGetAtomicNum('Br') : 1,
    }
    # Determine number of bonds.
    nbonds = len(list(mol.GetBonds()))
    # Compute maximum atomic valences.
    max_atom_valences = [ maximum_atomic_valences[atom.GetAtomicNum()] for atom in mol.GetAtoms() ]
    # Determine the number of bonds impinging on each atom.
    nbonds_per_atom = [ len(list(atom.GetBonds())) for atom in mol.GetAtoms() ]
    # Compute maximum bond orders based on maximum atomic valence at either end of the bond.
    max_bond_orders = [ min(
        max_atom_valences[bond.GetBgnIdx()] - nbonds_per_atom[bond.GetBgnIdx()] + 1,
        max_atom_valences[bond.GetEndIdx()] - nbonds_per_atom[bond.GetEndIdx()] + 1)
            for bond in mol.GetBonds() ]
    # Compute maximum number of possibilities we have to search.
    max_possibilities = np.prod( np.array(max_bond_orders) )
    print('There are %d possible bond order assignments to search.' % max_possibilities)

    def select_bond_orders(index, max_bond_orders):
        bond_orders = np.zeros([nbonds], np.int32)
        for (bond_index, max_bond_order) in enumerate(max_bond_orders):
            bond_orders[bond_index] = index % max_bond_order + 1
            index = int(index/max_bond_order)
        return bond_orders

    # Create a list of bonds impinging on each atom.
    bond_lists = [ list() for atom in mol.GetAtoms() ]
    for (bond_index, bond) in enumerate(mol.GetBonds()):
        bond_lists[bond.GetBgnIdx()].append(bond_index)
        bond_lists[bond.GetEndIdx()].append(bond_index)

    # Enumerate possibilities.
    valid_bond_orders = list()
    for index in range(max_possibilities):
        # Set possible bond orders.
        bond_orders = select_bond_orders(index, max_bond_orders)
        # Compute atomic valences.
        atom_valences = [ np.sum(bond_orders[bond_list]) for bond_list in bond_lists ]
        # Check that maximum valences are satisfied.
        # TODO: We also need to handle charge state variants, which may not satisfy maximum valences.
        if np.all(np.array(atom_valences) == np.array(max_atom_valences)):
            valid_bond_orders.append(bond_orders)
    print('There are %d valid bond order arrangements' % len(valid_bond_orders))
    print(valid_bond_orders)

    # Reduce these to molecules.
    valid_molecules = list()
    for bond_orders in valid_bond_orders:
        newmol = oechem.OEGraphMol(mol)
        for (bond, order) in zip(newmol.GetBonds(), bond_orders):
            bond.SetOrder(int(order))
        oechem.OEClearAromaticFlags(newmol);
        oechem.OEAssignAromaticFlags(newmol)
        oechem.OEAssignFormalCharges(newmol);
        valid_molecules.append(newmol)

    # TODO: Filter molecules to return only unique ones.

    return valid_molecules

def generateTopologyFromOEMol(molecule):
    """
    Generate an OpenMM Topology object from an OEMol molecule.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule from which a Topology object is to be generated.

    Returns
    -------
    topology : simtk.openmm.app.Topology
        The Topology object generated from `molecule`.

    """
    from simtk.openmm.app import Topology
    topology = Topology()
    chain = topology.addChain()
    resname = molecule.GetTitle()
    residue = topology.addResidue(resname)

    # Create atoms.
    for atom in molecule.GetAtoms():
        name = atom.GetName()
        element = Element.getByAtomicNumber(atom.GetAtomicNum())
        atom = topology.addAtom(name, element, residue)

    # Index atoms.
    atoms = { atom : atom.name for atom in topology.atoms() }

    # Create bonds.
    for bond in molecule.GetBonds():
        topology.addBond(atoms[bond.GetBgn().GetName()], atoms[bond.GetEnd().GetName()])

    return topology

def generateResidueTemplate(molecule, residue_atoms=None):
    """
    Generate an residue template for simtk.openmm.app.ForceField using GAFF/AM1-BCC.

    This requires the OpenEye toolkit.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule to be parameterized.
        The molecule must have explicit hydrogens.
        Charge will be inferred from the net formal charge.
    residue_atomset : set of OEAtom, optional, default=None
        If not None, only the atoms in this set will be used to construct the residue template

    Returns
    -------
    template : simtk.openmm.app.forcefield._TemplateData
        Residue template for ForceField using atom types and parameters from `gaff.xml`.
    additional_parameters_ffxml : str
        Contents of ForceField `ffxml` file defining additional parameters from parmchk(2).

    """
    # Generate a unique residue template name to avoid namespace collisions.
    # TODO: Can we come up with a more intelligent name?
    from uuid import uuid4
    template_name = uuid4()

    # Compute net formal charge.
    oechem.OEAssignFormalCharges(molecule)
    charges = [ atom.GetFormalCharge() for atom in molecule.GetAtoms() ]
    net_charge = np.array(charges).sum()

    # Generate canonical AM1-BCC charges and a reference conformation.
    molecule = get_charges(molecule, strictStereo=False, keep_confs=1)

    # Create temporary directory for running antechamber.
    import tempfile
    tmpdir = tempfile.mkdtemp()
    print('temporary directory: %s' % tmpdir) # DEBUG
    input_mol2_filename = os.path.join(tmpdir, template_name + '.tripos.mol2')
    gaff_mol2_filename = os.path.join(tmpdir, template_name + '.gaff.mol2')
    frcmod_filename = os.path.join(tmpdir, template_name + '.frcmod')

    # Write Tripos mol2 file as antechamber input.
    ofs = oemolostream(input_mol2_filename)
    oechem.OEWriteMolecule(ofs, molecule)
    ofs.close()

    # Parameterize the molecule with antechamber.
    run_antechamber(template_name, input_mol2_filename, charge_method=None, net_charge=net_charge, gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename)

    # Read the resulting GAFF mol2 file as a ParmEd structure.
    #structure = parmed.load_file(gaff_mol2_filename)
    #structure_atoms = { atom.name : atom for atom in structure.atoms }
    ifs = oemolistream(gaff_mol2_filename)
    ifs.SetFlavor(oechem.OEFormat_MOL2H, oechem.OEIFlavor_MOL2)
    oechem.OEReadMolecule(ifs, molecule)
    ifs.close()

    # If residue_atoms = None, add all atoms to the residues
    if residue_atoms == None:
        residue_atoms = [ atom for atom in molecule.GetAtoms() ]

    # Modify partial charges so that charge on residue atoms is integral.
    residue_charge = 0.0
    for atom in residue_atoms:
        partial_charge += atom.GetPartialCharge()
    factor = net_charge / residue_charge
    for atom in residue_atoms:
        atom.SetPartialCharge( factor * atom.GetPartialCharge() )

    # Create residue template.
    template = ForceField._TemplateData(template_name)
    for atom in molecule.GetAtoms():
        atomname = atom.GetName()
        typename = atom.GetType()
        element = Element.getByAtomicNumber(atom.GetAtomicNum())
        charge = atom.GetPartialCharge()
        parameters = { 'charge' : charge }
        #typename = structure_atoms[atomname].type # assigned GAFF atom type
        #parameters = { 'charge' : structure_atoms[atomname].charge } # partial atomic charge
        atom_template = ForceField._TemplateAtomData(atomname, typename, element, parameters)
        template.atoms.append(atom_template)
    for bond in molecule.GetBonds():
        if (bond.GetBgn() in residue_atoms) and (bond.GetEnd() in residue_atoms):
            template.addBondByName(bond.GetBgn().GetName(), bond.GetEnd().GetName())
        elif (bond.GetBgn() in residue_atoms) and (bond.GetEnd() not in residue_atoms):
            template.addExternalBondByName(bond.GetBgn().GetName())
        elif (bond.GetBgn() not in residue_atoms) and (bond.GetEnd() in residue_atoms):
            template.addExternalBondByName(bond.GetEnd().GetName())

    # Generate ffxml file contents for parmchk-generated frcmod output.
    leaprc = StringIO("parm = loadamberparams %s" % frcmod_filename)
    params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    ffxml = ""
    params.write(ffxml)

    return template, ffxml

def gaffTemplateGenerator(forcefield, residue, structure=None):
    """
    OpenMM ForceField residue template generator for GAFF/AM1-BCC.

    NOTE: This implementation currently only handles small molecules, not polymeric residues.
    NOTE: We presume we have already loaded the `gaff.xml` force definitions into ForceField.

    Parameters
    ----------
    forcefield : simtk.openmm.app.ForceField
        The ForceField object to which residue templates and/or parameters are to be added.
    residue : simtk.openmm.app.Topology.Residue
        The residue topology for which a template is to be generated.

    Returns
    -------
    success : bool
        If the generator is able to successfully parameterize the residue, `True` is returned.
        If the generator cannot parameterize the residue, it should return `False` and not modify `forcefield`.

    """
    # Get a list of external bonds.
    external_bonds = [ bond for bond in residue.external_bonds() ]
    if len(external_bonds) > 0:
        # We can't parameterize residues with external bonds right now.
        return False

    # Generate a unique residue template name to avoid namespace collisions.
    # TODO: Can we come up with a more intelligent name?
    from uuid import uuid4
    template_name = uuid4()

    # Create ParmEd structure for residue.
    #structure = parmed.Structure()
    #for a in residue.atoms():
    #    if a.element is None:
    #        atom = parmed.ExtraPoint(name=a.name)
    #    else:
    #        atom = parmed.Atom(atomic_number=a.element.atomic_number, name=a.name, mass=a.element.mass)
    #    struct.add_atom(atom, residue.name, residue.index, 'A')
    #    atommap[a] = atom
    #for a1, a2 in topology.bonds():
    #    struct.bonds.append(Bond(atommap[a1], atommap[a2]))

    # Create an OpenEye molecule.
    from openeye import oechem
    mol = oechem.OEGraphMol()
    mol.SetTitle(template_name)
    oeatoms = dict()
    for atom in residue.atoms():
        oeatoms[atom] = mol.NewAtom(atom.element.atomic_number)
        oeatoms[atom].SetName(atom.name)
    for (atom1, atom2) in residue.internal_bonds():
        mol.NewBond(oeatoms[atom1], oeatoms[atom2])
    # Use atoms and bond graph to perceive bond orders and aromaticity.
    # TODO: Is this sufficient for determining bond orders and aromaticity correctly?
    # NOTE: No stereochemistry is defined here. This may cause problems with molecules with multiple chiral centers.
    oechem.OEFindRingAtomsAndBonds(mol)
    oechem.OEPerceiveBondOrders(mol)
    oechem.OEAssignFormalCharges(mol)
    oechem.OEAssignAromaticFlags(mol)
    #oechem.OEAssignHybridization(mol)
    #oechem.OEAssignAromaticFlags(mol)

    # Compute net formal charge.
    net_charge = 0.0
    for oeatom in oeatoms:
        net_charge += oeatom.GetFormalCharge()

    # Parameterize molecule using antechamber.
    import tempfile
    tmpdir = tempfile.mkdtemp()
    print('temporary directory: %s' % tmpdir) # DEBUG
    input_mol2_filename = os.path.join(tmpdir, template_name + '.tripos.mol2')
    gaff_mol2_filename = os.path.join(tmpdir, template_name + '.gaff.mol2')
    frcmod_filename = os.path.join(tmpdir, template_name + '.frcmod')

    # Render molecule
    from openeye import oedepict
    opts = oedepict.OE2DMolDisplayOptions()
    disp = oedepict.OE2DMolDisplay(mol, opts)
    pdf_filename = os.path.join(tmpdir, template_name + '.pdf')
    ofs = oeofstream(pdf_filename)
    oedepict.OERenderMolecule(ofs, 'pdf', disp)
    ofs.close()

    # Generate canonical AM1-BCC charges and a reference conformation.
    from openmoltools.openeye import get_charges
    mol = get_charges(mol, strictStereo=False, keep_confs=1)

    # Write the Tripos mol2 file.
    ofs = oemolostream(input_mol2_filename)
    oechem.OEWriteMolecule(ofs, mol)
    ofs.close()

    # Parameterize the molecule with antechamber.
    run_antechamber(template_name, input_mol2_filename, charge_method=None, net_charge=net_charge, gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename)

    # Read the resulting GAFF mol2 file.
    structure = parmed.load_file(gaff_mol2_filename)
    structure_atoms = { atom.name : atom for atom in structure.atoms }

    # Create residue template.
    template = ForceField._TemplateData(template_name)
    for atom in residue.atoms():
        typename = structure_atoms[atom.name].type # assigned GAFF atom type
        parameters = { 'charge' : structure_atoms[atom.name].charge } # partial atomic charge
        atom_template = ForceField._TemplateAtomData(atom.name, typename, atom.element, parameters)
        template.atoms.append(atom_template)
    for (atom1,atom2) in residue.internal_bonds():
        template.addBondByName(atom1.name, atom2.name)
    residue_atoms = [ atom for atom in residue.atoms() ]
    for (atom1,atom2) in residue.external_bonds():
        if atom1 in residue_atoms:
            template.addExternalBondByName(atom1.name)
        elif atom2 in residue_atoms:
            template.addExternalBondByName(atom2.name)

    # Register the template.
    forcefield.registerResidueTemplate(template)

    # Signal that we have successfully parameterized the residue.
    return True
