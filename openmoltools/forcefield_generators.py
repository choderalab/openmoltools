#!/usr/bin/env python
"""
OpenMM ForceField residue template generators.

"""

from __future__ import absolute_import

import numpy as np
import os, os.path, sys
from simtk.openmm.app import ForceField
from openmoltools.amber import run_antechamber
from openmoltools.openeye import get_charges
from simtk.openmm.app.element import Element
import parmed
if sys.version_info >= (3, 0):
    from io import StringIO
    from subprocess import getstatusoutput
else:
    from cStringIO import StringIO
    from commands import getstatusoutput

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
    # Create a Topology object with one Chain and one Residue.
    from simtk.openmm.app import Topology
    topology = Topology()
    chain = topology.addChain()
    resname = molecule.GetTitle()
    residue = topology.addResidue(resname, chain)

    # Create atoms in the residue.
    for atom in molecule.GetAtoms():
        name = atom.GetName()
        element = Element.getByAtomicNumber(atom.GetAtomicNum())
        atom = topology.addAtom(name, element, residue)

    # Create bonds.
    atoms = { atom.name : atom for atom in topology.atoms() }
    for bond in molecule.GetBonds():
        topology.addBond(atoms[bond.GetBgn().GetName()], atoms[bond.GetEnd().GetName()])

    return topology

def generateOEMolFromTopologyResidue(residue):
    """
    Generate an OpenEye OEMol molecule from an OpenMM Topology Residue.

    Parameters
    ----------
    residue : simtk.openmm.app.topology.Residue
        The topology Residue from which an OEMol is to be created.
        An Exception will be thrown if this residue has external bonds.

    Returns
    -------
    molecule : openeye.oechem.OEMol
        The OEMol molecule corresponding to the topology.
        Atom order will be preserved and bond orders assigned.

    The Antechamber `bondtype` program will be used to assign bond orders, and these
    will be converted back into OEMol bond type assignments.

    Note that there is no way to preserve stereochemistry since `Residue` does
    not note stereochemistry in any way.

    """
    # Raise an Exception if this residue has external bonds.
    if len(list(residue.external_bonds())) > 0:
        raise Exception("Cannot generate an OEMol from residue '%s' because it has external bonds." % residue.name)

    from openeye import oechem
    # Create OEMol where all atoms have bond order 1.
    molecule = oechem.OEMol()
    molecule.SetTitle(residue.name) # name molecule after first residue
    for atom in residue.atoms():
        oeatom = molecule.NewAtom(atom.element.atomic_number)
        oeatom.SetName(atom.name)
    oeatoms = { oeatom.GetName() : oeatom for oeatom in molecule.GetAtoms() }
    for (atom1, atom2) in residue.bonds():
        order = 1
        molecule.NewBond(oeatoms[atom1.name], oeatoms[atom2.name], order)

    # Write out a mol2 file without altering molecule.
    import tempfile
    tmpdir = tempfile.mkdtemp()
    mol2_input_filename = os.path.join(tmpdir,'molecule-before-bond-perception.mol2')
    ac_output_filename = os.path.join(tmpdir,'molecule-after-bond-perception.ac')
    ofs = oechem.oemolostream(mol2_input_filename)
    m2h = True
    substruct = False
    oechem.OEWriteMol2File(ofs, molecule, m2h, substruct)
    ofs.close()
    # Run Antechamber bondtype
    import subprocess
    #command = 'bondtype -i %s -o %s -f mol2 -j full' % (mol2_input_filename, ac_output_filename)
    command = 'antechamber -i %s -fi mol2 -o %s -fo ac -j 2' % (mol2_input_filename, ac_output_filename)
    [status, output] = getstatusoutput(command)

    # Define mapping from GAFF bond orders to OpenEye bond orders.
    order_map = { 1 : 1, 2 : 2, 3: 3, 7 : 1, 8: 2, 9 : 5, 10 : 5 }
    # Read bonds.
    infile = open(ac_output_filename)
    lines = infile.readlines()
    infile.close()
    antechamber_bond_types = list()
    for line in lines:
        elements = line.split()
        if elements[0] == 'BOND':
            antechamber_bond_types.append(int(elements[4]))
    for (bond, antechamber_bond_type) in zip(molecule.GetBonds(), antechamber_bond_types):
        bond.SetOrder(order_map[antechamber_bond_type])

    # Clean up.
    os.unlink(mol2_input_filename)
    os.unlink(ac_output_filename)
    os.rmdir(tmpdir)

    # Set aromaticity.
    # TODO: Is this necessary?
    oechem.OEClearAromaticFlags(molecule)
    oechem.OEAssignAromaticFlags(molecule, oechem.OEAroModelOpenEye)

    # Assign Tripos atom types.
    # TODO: Is this necessary?
    oechem.OETriposAtomTypeNames(molecule)
    oechem.OETriposBondTypeNames(molecule)

    # Assign geometry
    # TODO: Is this necessary?
    from openeye import oeomega
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(False)
    omega.SetStrictStereo(False)
    omega(molecule)

    return molecule

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

    Note that this method preserves stereochemistry during AM1-BCC charge parameterization.

    """
    # Generate a unique residue template name to avoid namespace collisions.
    # TODO: Can we come up with a more intelligent name?
    #from uuid import uuid4
    #template_name = str(uuid4())
    template_name = molecule.GetTitle()

    # Compute net formal charge.
    from openeye import oechem
    oechem.OEAssignFormalCharges(molecule)
    charges = [ atom.GetFormalCharge() for atom in molecule.GetAtoms() ]
    net_charge = np.array(charges).sum()

    # Generate canonical AM1-BCC charges and a reference conformation.
    molecule = get_charges(molecule, strictStereo=False, keep_confs=1)

    # Create temporary directory for running antechamber.
    import tempfile
    tmpdir = tempfile.mkdtemp()
    input_mol2_filename = os.path.join(tmpdir, template_name + '.tripos.mol2')
    gaff_mol2_filename = os.path.join(tmpdir, template_name + '.gaff.mol2')
    frcmod_filename = os.path.join(tmpdir, template_name + '.frcmod')

    # Write Tripos mol2 file as antechamber input.
    ofs = oechem.oemolostream(input_mol2_filename)
    oechem.OEWriteMolecule(ofs, molecule)
    ofs.close()

    # Parameterize the molecule with antechamber.
    run_antechamber(template_name, input_mol2_filename, charge_method=None, net_charge=net_charge, gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename)

    # Read the resulting GAFF mol2 file as a ParmEd structure.
    ifs = oechem.oemolistream(gaff_mol2_filename)
    ifs.SetFlavor(oechem.OEFormat_MOL2, oechem.OEIFlavor_MOL2_DEFAULT | oechem.OEIFlavor_MOL2_M2H | oechem.OEIFlavor_MOL2_Forcefield)
    m2h = True
    oechem.OEReadMolecule(ifs, molecule)
    ifs.close()

    # If residue_atoms = None, add all atoms to the residues
    if residue_atoms == None:
        residue_atoms = [ atom for atom in molecule.GetAtoms() ]

    # Modify partial charges so that charge on residue atoms is integral.
    residue_charge = 0.0
    sum_of_absolute_charge = 0.0
    for atom in residue_atoms:
        charge = atom.GetPartialCharge()
        residue_charge += charge
        sum_of_absolute_charge += abs(charge)
    excess_charge = residue_charge - net_charge
    if sum_of_absolute_charge == 0.0:
        sum_of_absolute_charge = 1.0
    for atom in residue_atoms:
        charge = atom.GetPartialCharge()
        atom.SetPartialCharge( charge + excess_charge * (abs(charge) / sum_of_absolute_charge) )

    # Create residue template.
    template = ForceField._TemplateData(template_name)
    for (index, atom) in enumerate(molecule.GetAtoms()):
        atomname = atom.GetName()
        typename = atom.GetType()
        element = Element.getByAtomicNumber(atom.GetAtomicNum())
        charge = atom.GetPartialCharge()
        parameters = { 'charge' : charge }
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
    ffxml = StringIO()
    params.write(ffxml)

    return template, ffxml.getvalue()

def createStructureFromResidue(residue):
    # Create ParmEd structure for residue.
    structure = parmed.Structure()
    for a in residue.atoms():
        if a.element is None:
            atom = parmed.ExtraPoint(name=a.name)
        else:
            atom = parmed.Atom(atomic_number=a.element.atomic_number, name=a.name, mass=a.element.mass)
        structure.add_atom(atom, residue.name, residue.index, 'A')
        atommap[a] = atom
    for a1, a2 in topology.bonds():
        structure.bonds.append(Bond(atommap[a1], atommap[a2]))

    return structure

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

    Note that there is no way to preserve stereochemistry since `Residue` does not specify stereochemistry in any way.
    Charge fitting is therefore performed on an indeterminate stereo form.

    """
    # Get a list of external bonds.
    external_bonds = [ bond for bond in residue.external_bonds() ]
    if len(external_bonds) > 0:
        # We can't parameterize residues with external bonds right now.
        return False

    # Generate an OpenEye OEMol molecule from the Topology Residue.
    molecule = generateOEMolFromTopologyResidue(residue)

    # Generate template and parameters.
    [template, ffxml] = generateResidueTemplate(molecule)

    # Register the template.
    forcefield.registerResidueTemplate(template)

    # Add the parameters.
    forcefield.loadFile(StringIO(ffxml))

    # Signal that we have successfully parameterized the residue.
    return True
