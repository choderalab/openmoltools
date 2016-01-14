"""
OpenMM ForceField residue template generators.

"""

import os, os.path
from simtk.openmm.app import ForceField
from openmoltools.amber import run_antechamber
import parmed
import oechem

def gaffTemplateGenerator(forcefield, residue):
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
    oechem.OEAssignHybridization(mol)
    oechem.OEAssignAromaticFlags(mol)

    # Compute net formal charge.
    oechem.OEAssignFormalCharges(mol)
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
    mol = get_charges(mol, strictStereo=False, keep_confs=1)

    # Write the Tripos mol2 file.
    #parmed.formats.mol2.write(structure, input_mol2_filename)
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


