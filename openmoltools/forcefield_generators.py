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
import json
import parmed
if sys.version_info >= (3, 0):
    from io import StringIO
    from subprocess import getstatusoutput, call
    def run_command(command):
        call(command.split())
else:
    from cStringIO import StringIO
    from commands import getstatusoutput
    def run_command(command):
        getstatusoutput(command)

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

def _ensureUniqueAtomNames(molecule):
    """
    Ensure all atom names are unique and not blank.
    If any atom names are degenerate or blank, Tripos atom names are assigned to all atoms.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule to be modified

    """
    from openeye import oechem
    atom_names = set()
    atom_names_are_unique = True
    for atom in molecule.GetAtoms():
        atom_name = atom.GetName()
        if (atom_name in atom_names) or (atom_name == ""):
            atom_names_are_unique = False
        atom_names.add(atom_name)
    if not atom_names_are_unique:
        oechem.OETriposAtomNames(molecule)

def generateOEMolFromTopologyResidue(residue, geometry=False, tripos_atom_names=False):
    """
    Generate an OpenEye OEMol molecule from an OpenMM Topology Residue.

    Parameters
    ----------
    residue : simtk.openmm.app.topology.Residue
        The topology Residue from which an OEMol is to be created.
        An Exception will be thrown if this residue has external bonds.
    geometry : bool, optional, default=False
        If True, will generate a single configuration with OEOmega.
        Note that stereochemistry will be *random*.
    tripos_atom_names : bool, optional, default=False
        If True, will generate and assign Tripos atom names.

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
        oeatom.AddData("topology_index", atom.index)
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
    run_command(command)

    # Define mapping from GAFF bond orders to OpenEye bond orders.
    order_map = { 1 : 1, 2 : 2, 3: 3, 7 : 1, 8 : 2, 9 : 5, 10 : 5 }
    # Read bonds.
    infile = open(ac_output_filename)
    lines = infile.readlines()
    infile.close()
    antechamber_bond_types = list()
    for line in lines:
        elements = line.split()
        if elements[0] == 'BOND':
            antechamber_bond_types.append(int(elements[4]))
    oechem.OEClearAromaticFlags(molecule)
    for (bond, antechamber_bond_type) in zip(molecule.GetBonds(), antechamber_bond_types):
        #bond.SetOrder(order_map[antechamber_bond_type])
        bond.SetIntType(order_map[antechamber_bond_type])
    oechem.OEFindRingAtomsAndBonds(molecule)
    oechem.OEKekulize(molecule)
    oechem.OEAssignFormalCharges(molecule)
    oechem.OEAssignAromaticFlags(molecule, oechem.OEAroModelOpenEye)

    # Clean up.
    os.unlink(mol2_input_filename)
    os.unlink(ac_output_filename)
    os.rmdir(tmpdir)

    # Generate Tripos atom names if requested.
    if tripos_atom_names:
        oechem.OETriposAtomNames(molecule)

    # Assign geometry
    if geometry:
        from openeye import oeomega
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(1)
        omega.SetIncludeInput(False)
        omega.SetStrictStereo(False)
        omega(molecule)

    return molecule


def _computeNetCharge(molecule):
    """
    Compute the net formal charge on the molecule.
    Formal charges are assigned by this function.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule for which a net formal charge is to be computed

    Returns
    -------
    net_charge : float
        The net formal charge on the molecule

    """
    from openeye import oechem
    oechem.OEAssignFormalCharges(molecule)
    charges = [ atom.GetFormalCharge() for atom in molecule.GetAtoms() ]
    net_charge = np.array(charges).sum()
    return net_charge


def _writeMolecule(molecule, output_filename, standardize=True):
    """
    Write the molecule to a file.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule to write (will be modified by writer).
    output_filename : str
        The filename of file to be written; type is autodetected by extension.
    standardize : bool, optional, default=True
        Standardize molecular properties such as atom names in the output file.

    """
    from openmoltools.openeye import molecule_to_mol2
    molecule_to_mol2(molecule, tripos_mol2_filename=output_filename, conformer=0, residue_name=molecule.GetTitle(), standardize=standardize)

def generateResidueTemplate(molecule, residue_atoms=None, normalize=True, gaff_version='gaff'):
    """
    Generate an residue template for simtk.openmm.app.ForceField using GAFF and AM1-BCC ELF10 charges.

    This requires the OpenEye toolkit.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule to be parameterized.
        The molecule must have explicit hydrogens.
        Net charge will be inferred from the net formal charge on each molecule.
        Partial charges will be determined automatically using oequacpac and canonical AM1-BCC charging rules.
    residue_atomset : set of OEAtom, optional, default=None
        If not None, only the atoms in this set will be used to construct the residue template
    normalize : bool, optional, default=True
        If True, normalize the molecule by checking aromaticity, adding
        explicit hydrogens, and renaming by IUPAC name.
    gaff_version : str, default = 'gaff'
        One of ['gaff', 'gaff2']; selects which atom types to use.


    Returns
    -------
    template : simtk.openmm.app.forcefield._TemplateData
        Residue template for ForceField using atom types and parameters from `gaff.xml` or `gaff2.xml`.
    additional_parameters_ffxml : str
        Contents of ForceField `ffxml` file defining additional parameters from parmchk(2).

    Notes
    -----
    The residue template will be named after the molecule title.
    This method preserves stereochemistry during AM1-BCC charge parameterization.
    Atom names in molecules will be assigned Tripos atom names if any are blank or not unique.

    """
    # Make a copy of the molecule so we don't modify the original
    molecule = molecule.CreateCopy()

    # Set the template name based on the molecule title plus a globally unique UUID.
    from uuid import uuid4
    template_name = molecule.GetTitle() + '-' + str(uuid4())

    # If any atom names are not unique, atom names
    _ensureUniqueAtomNames(molecule)

    # Compute net formal charge.
    net_charge = _computeNetCharge(molecule)

    # Generate canonical AM1-BCC ELF10 charges
    from openeye import oequacpac
    oequacpac.OEAssignCharges(molecule, oequacpac.OEAM1BCCELF10Charges())

    # Set title to something that antechamber can handle
    molecule.SetTitle('MOL')

    # Create temporary directory for running antechamber.
    import tempfile
    tmpdir = tempfile.mkdtemp()
    prefix = 'molecule'
    input_mol2_filename = os.path.join(tmpdir, prefix + '.tripos.mol2')
    gaff_mol2_filename = os.path.join(tmpdir, prefix + '.gaff.mol2')
    frcmod_filename = os.path.join(tmpdir, prefix + '.frcmod')

    # Write Tripos mol2 file as antechamber input.
    _writeMolecule(molecule, input_mol2_filename, standardize=normalize)

    # Parameterize the molecule with antechamber.
    run_antechamber(template_name, input_mol2_filename, charge_method=None, net_charge=net_charge, gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename, gaff_version=gaff_version)

    # Read the resulting GAFF mol2 file as a ParmEd structure
    from openeye import oechem
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

    # Generate additional parameters, if needed
    # TODO: Do we have to make sure that we don't duplicate existing parameters already loaded in the forcefield?
    from inspect import signature # use introspection to support multiple parmed versions
    leaprc = StringIO('parm = loadamberparams %s' % frcmod_filename)
    params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
    kwargs = {}
    if 'remediate_residues' in signature(parmed.openmm.OpenMMParameterSet.from_parameterset).parameters:
        kwargs['remediate_residues'] = False
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params, **kwargs)
    ffxml = StringIO()
    kwargs = {}
    if 'write_unused' in signature(params.write).parameters:
        kwargs['write_unused'] = True
    params.write(ffxml, **kwargs)
    additional_parameters_ffxml = ffxml.getvalue()

    return template, additional_parameters_ffxml

def generateForceFieldFromMolecules(molecules, ignoreFailures=False, generateUniqueNames=False, normalize=True, gaff_version='gaff'):
    """
    Generate ffxml file containing additional parameters and residue templates for simtk.openmm.app.ForceField using GAFF/AM1-BCC.

    This requires the OpenEye toolkit.

    Parameters
    ----------
    molecules : list of openeye.oechem.OEMol
        The molecules to be parameterized.
        All molecules must have explicit hydrogens.
        Net charge will be inferred from the net formal charge on each molecule.
        Partial charges will be determined automatically using oequacpac and canonical AM1-BCC charging rules.
    ignoreFailures: bool, optional, default=False
        Determines whether to add a failed molecule to the list of failed molecules (True),
        or raise an Exception (False).
    generateUniqueNames : bool, optional, default=False
        If True, will generate globally unique names for templates.
    normalize : bool, optional, default=True
        If True, normalize the molecule by checking aromaticity, adding
        explicit hydrogens, and renaming by IUPAC name.
    gaff_version : str, default = 'gaff'
        One of ['gaff', 'gaff2']; selects which atom types to use.

    Returns
    -------
    ffxml : str
        Contents of ForceField `ffxml` file defining additional parameters from parmchk(2) and residue templates.
    failed_molecule_list : list of openeye.oechem.OEMol
        List of the oemols that could not be parameterized. Only returned if ignoreFailures=True

    Notes
    -----
    This method preserves stereochemistry during AM1-BCC charge parameterization.
    Residue template names will be set from molecule names.
    Atom names in molecules will be assigned Tripos atom names if any are blank or not unique.

    """
    if not generateUniqueNames:
    # Check template names are unique.
        template_names = set()
        for molecule in molecules:
            template_name = molecule.GetTitle()
            if template_name == '<0>':
                raise Exception("Molecule '%s' has invalid name" % template_name)
            if template_name in template_names:
                raise Exception("Molecule '%s' has template name collision." % template_name)
            template_names.add(template_name)

    # Process molecules.
    import tempfile
    tmpdir = tempfile.mkdtemp()
    olddir = os.getcwd()
    os.chdir(tmpdir)
    leaprc = ""
    failed_molecule_list = []
    for (molecule_index, molecule) in enumerate(molecules):
        # Set the template name based on the molecule title.
        if generateUniqueNames:
            from uuid import uuid4
            template_name = molecule.GetTitle() + '-' + str(uuid4())
        else:
            template_name = molecule.GetTitle()

        # If any atom names are not unique, atom names
        _ensureUniqueAtomNames(molecule)

        # Compute net formal charge.
        net_charge = _computeNetCharge(molecule)

        # Generate canonical AM1-BCC charges and a reference conformation.
        if not ignoreFailures:
            molecule = get_charges(molecule, strictStereo=False, keep_confs=1, normalize=normalize)
        else:
            try:
                molecule = get_charges(molecule, strictStereo=False, keep_confs=1, normalize=normalize)
            except:
                failed_molecule_list.append(molecule)

        # Create a unique prefix.
        prefix = 'molecule%010d' % molecule_index

        # Create temporary directory for running antechamber.
        input_mol2_filename = prefix + '.tripos.mol2'
        gaff_mol2_filename  = prefix + '.gaff.mol2'
        frcmod_filename     = prefix + '.frcmod'

        # Write Tripos mol2 file as antechamber input.
        _writeMolecule(molecule, input_mol2_filename, standardize=normalize)

        # Parameterize the molecule with antechamber.
        run_antechamber(prefix, input_mol2_filename, charge_method=None, net_charge=net_charge, gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename, gaff_version=gaff_version)

        # Append to leaprc input for parmed.
        leaprc += '%s = loadmol2 %s\n' % (prefix, gaff_mol2_filename)
        leaprc += 'loadamberparams %s\n' % frcmod_filename

    # Generate ffxml file contents for parmchk-generated frcmod output.
    leaprc = StringIO(leaprc)
    params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    ffxml = StringIO()
    params.write(ffxml)

    # TODO: Clean up temporary directory.
    os.chdir(olddir)

    if ignoreFailures:
        return ffxml.getvalue(), failed_molecule_list
    else:
        return ffxml.getvalue()

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

class OEGAFFTemplateGenerator(object):
    """
    OpenMM ForceField residue template generator for GAFF/AM1-BCC using pre-cached OpenEye toolkit OEMols.

    Examples
    --------

    Create a template generator for GAFF for a single OEMol and register it with ForceField:

    >>> from openmoltools.forcefield_generators import OEGAFFTemplateGenerator
    >>> template_generator = OEGAFFTemplateGenerator(oemols=oemol)
    >>> from simtk.openmm.app import ForceField
    >>> forcefield = ForceField('gaff.xml', 'amber14-all.xml', 'tip3p.xml')
    >>> forcefield.registerTemplateGenerator(template_generator.generator)

    Create a template generator for GAFF2 for multiple OEMols:

    >>> template_generator = OEGAFFTemplateGenerator(oemols=[oemol1, oemol2], gaff_version='gaff2')

    You can also some OEMols later on after the generator has been registered:

    >>> forcefield.add_oemols(oemol)
    >>> forcefield.add_oemols([oemol1, oemol2])

    You can optionally create or use a tiny database cache of pre-parameterized molecules:

    >>> template_generator = OEGAFFTemplateGenerator(cache='gaff-molecules.json')

    Newly parameterized molecules will be written to the cache, saving time next time!

    """
    def __init__(self, oemols=None, cache=None, gaff_version='gaff'):
        """
        Create an OEGAFFTemplateGenerator with some OpenEye toolkit OEMols

        Requies the OpenEye Toolkit.

        Parameters
        ----------
        oemols : OEMol or list of OEMol, optional, default=None
            If specified, these molecules will be recognized and parameterized with antechamber as needed.
            The parameters will be cached in case they are encountered again the future.
        cache : str or TinyDB instance, optional, default=None
            Filename or TinyDB instance for global caching of parameters.
            If specified, parameterized molecules will be stored in a TinyDB instance.
            Note that no checking is done to determine this cache was created with the same GAFF version.
        gaff_version : str, default = 'gaff'
            One of ['gaff', 'gaff2']; selects which GAFF major version to use.

        .. todo :: Should we support SMILES instead of OEMols?

        Examples
        --------

        Create a template generator for GAFF for a single OEMol and register it with ForceField:

        >>> from openmoltools.forcefield_generators import OEGAFFTemplateGenerator
        >>> template_generator = OEGAFFTemplateGenerator(oemols=oemol)

        Create a template generator for GAFF2 for multiple OEMols:

        >>> template_generator = OEGAFFTemplateGenerator(oemols=[oemol1, oemol2], gaff_version='gaff2')

        You can optionally create or use a tiny database cache of pre-parameterized molecules:

        >>> template_generator = OEGAFFTemplateGenerator(cache='gaff-molecules.json')

        """
        from openeye import oechem

        self._gaff_version = gaff_version

        # Add oemols to the dictionary
        self._oemols = dict()
        self.add_oemols(oemols)

        # Open TinyDB instance
        from tinydb import TinyDB
        if isinstance(cache, TinyDB):
            self._db = cache
        elif isinstance(cache, str):
            # Open the database instance
            self._db = TinyDB(cache)
        else:
            self._db = None

        self._smiles_added_to_db = set() # set of SMILES added to the database this session

    # TODO: Replace this encoder/decoder logic when openmm objects are properly serializable
    class _JSONEncoder(json.JSONEncoder):
        def default(self, o):
            from simtk.openmm.app import ForceField, Element
            if isinstance(o, ForceField._TemplateData):
                s = {'_type' : '_TemplateData'}
                s.update(o.__dict__)
                return s
            elif isinstance(o, ForceField._TemplateAtomData):
                s = {'_type' : '_TemplateAtomData'}
                s.update(o.__dict__)
                return s
            elif isinstance(o, Element):
                return {'_type' : 'Element', 'atomic_number' : o.atomic_number}
            else:
                return super(OEGAFFTemplateGenerator._JSONEncoder, self).default(o)

    class _JSONDecoder(json.JSONDecoder):
        def __init__(self, *args, **kwargs):
            json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

        def object_hook(self, obj):
            from simtk.openmm.app import ForceField, Element
            if '_type' not in obj:
                return obj
            type = obj['_type']
            if type == '_TemplateData':
                template = ForceField._TemplateData.__new__(ForceField._TemplateData)
                del obj['_type']
                template.__dict__ = obj
                return template
            if type == '_TemplateAtomData':
                atom = ForceField._TemplateAtomData.__new__(ForceField._TemplateAtomData)
                del obj['_type']
                atom.__dict__ = obj
                return atom
            elif type == 'Element':
                return Element.getByAtomicNumber(obj['atomic_number'])
            return obj

    def add_oemols(self, oemols=None):
        """
        Add specified list of OEMol objects to cached molecules that will be recognized.

        Parameters
        ----------
        oemols : OEMol or list of OEMol, optional, default=None
            If specified, these molecules will be recognized and parameterized with antechamber as needed.
            The parameters will be cached in case they are encountered again the future.

        Examples
        --------
        Add some OEMols later on after the generator has been registered:

        >>> forcefield.add_oemols(oemol)
        >>> forcefield.add_oemols([oemol1, oemol2])

        """
        # Return if empty
        if not oemols:
            return

        # Ensure oemols is iterable
        try:
            iterator = iter(oemols)
        except TypeError as te:
            oemols = [ oemols ]

        # Create copies
        oemols = [ oemol.CreateCopy() for oemol in oemols ]

        # Cache OEMols
        from openeye import oechem
        self._oemols.update( { oechem.OEMolToSmiles(oemol) : oemol for oemol in oemols } )

    @staticmethod
    def _match_residue(residue, oemol_template):
        """Determine whether a residue matches an OEMol template and return a list of corresponding atoms.

        This implementation uses NetworkX for graph isomorphism determination.

        Parameters
        ----------
        residue : simtk.openmm.app.topology.Residue
            The residue to check
        oemol_template : openeye.oechem.OEMol
            The OEMol template to compare it to

        Returns
        -------
        matches : dict of int : int
            matches[residue_atom_index] is the corresponding OEMol template atom index
            or None if it does not match the template

        """
        import networkx as nx
        from openeye import oechem

        # Make a copy of the template
        oemol_template = oemol_template.CreateCopy()

        # Ensure atom names are unique
        oechem.OETriposAtomNames(oemol_template)

        # Build list of external bonds for residue
        number_of_external_bonds = { atom : 0 for atom in residue.atoms() }
        for bond in residue.external_bonds():
            if bond[0] in number_of_external_bonds: number_of_external_bonds[bond[0]] += 1
            if bond[1] in number_of_external_bonds: number_of_external_bonds[bond[1]] += 1

        # Residue graph
        residue_graph = nx.Graph()
        for atom in residue.atoms():
            residue_graph.add_node(atom, element=atom.element.atomic_number, number_of_external_bonds=number_of_external_bonds[atom])
        for bond in residue.internal_bonds():
            residue_graph.add_edge(bond[0], bond[1])

        # Template graph
        # TODO: We can support templates with "external" bonds or atoms using attached string data in future
        # See https://docs.eyesopen.com/toolkits/python/oechemtk/OEChemClasses/OEAtomBase.html
        template_graph = nx.Graph()
        for oeatom in oemol_template.GetAtoms():
            template_graph.add_node(oeatom.GetName(), element=oeatom.GetAtomicNum(), number_of_external_bonds=0)
        for oebond in oemol_template.GetBonds():
            template_graph.add_edge(oebond.GetBgn().GetName(), oebond.GetEnd().GetName())

        # Determine graph isomorphism
        from networkx.algorithms import isomorphism
        graph_matcher = isomorphism.GraphMatcher(residue_graph, template_graph)
        if graph_matcher.is_isomorphic() == False:
            return None

        # Translate to local residue atom indices
        atom_index_within_residue = { atom : index for (index, atom) in enumerate(residue.atoms()) }
        atom_index_within_template = { oeatom.GetName() : index for (index, oeatom) in enumerate(oemol_template.GetAtoms()) }
        matches = { atom_index_within_residue[residue_atom] : atom_index_within_template[template_atom] for (residue_atom, template_atom) in graph_matcher.mapping.items() }

        return matches

    def generator(self, forcefield, residue, structure=None):
        """
        Residue template generator method to register with simtk.openmm.app.ForceField

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
        from openeye import oechem

        # If a database is specified, check against molecules in the database
        if self._db is not None:
            for entry in self._db:
                # Skip any molecules we've added to the database this session
                if entry['smiles'] in self._smiles_added_to_db:
                    continue

                # See if the template matches
                oemol_template = oechem.OEMol()
                oechem.OESmilesToMol(oemol_template, entry['smiles'])
                oechem.OEAddExplicitHydrogens(oemol_template)
                if self._match_residue(residue, oemol_template):
                    # Register the template
                    template = self._JSONDecoder().decode(entry['template'])
                    forcefield.registerResidueTemplate(template)
                    # Add the parameters
                    # TODO: Do we have to worry about parameter collisions?
                    forcefield.loadFile(StringIO(entry['ffxml']))
                    # Signal success
                    return True

        # Check against the molecules we know about
        for smiles, oemol_template in self._oemols.items():
            # See if the template matches
            if self._match_residue(residue, oemol_template):
                # Generate template and parameters.
                [template, ffxml] = generateResidueTemplate(oemol_template, gaff_version=self._gaff_version)
                # Register the template
                forcefield.registerResidueTemplate(template)
                # Add the parameters
                # TODO: Do we have to worry about parameter collisions?
                forcefield.loadFile(StringIO(ffxml))
                # Add it to the database, if we have one defined
                if self._db is not None:
                    self._db.insert({'smiles' : smiles, 'template' : self._JSONEncoder().encode(template), 'ffxml' : ffxml})
                    self._smiles_added_to_db.add(smiles)

                # Signal success
                return True

        # Report that we have failed to parameterize the residue
        return False

class SystemGenerator(object):
    """
    Utility factory to generate OpenMM Systems from Topology objects.

    Parameters
    ----------
    forcefields_to_use : list of string
        List of the names of ffxml files that will be used in system creation.
    forcefield_kwargs : dict of arguments to createSystem, optional
        Allows specification of various aspects of system creation.
    use_gaff : bool, optional, default=True
        If True, will add the GAFF residue template generator.

    Examples
    --------
    >>> from simtk.openmm import app
    >>> forcefield_kwargs={ 'nonbondedMethod' : app.NoCutoff, 'implicitSolvent' : None, 'constraints' : None }
    >>> system_generator = SystemGenerator(['amber99sbildn.xml'], forcefield_kwargs=forcefield_kwargs)
    >>> from openmmtools.testsystems import AlanineDipeptideVacuum
    >>> testsystem = AlanineDipeptideVacuum()
    >>> system = system_generator.createSystem(testsystem.topology)
    """

    def __init__(self, forcefields_to_use, forcefield_kwargs=None, use_gaff=True):
        self._forcefield_xmls = forcefields_to_use
        self._forcefield_kwargs = forcefield_kwargs if forcefield_kwargs is not None else {}
        from simtk.openmm.app import ForceField
        self._forcefield = ForceField(*self._forcefield_xmls)
        if use_gaff:
            self._forcefield.registerTemplateGenerator(gaffTemplateGenerator)

    def getForceField(self):
        """
        Return the associated ForceField object.

        Returns
        -------
        forcefield : simtk.openmm.app.ForceField
            The current ForceField object.
        """
        return self._forcefield

    def createSystem(self, topology):
        """
        Build a system from specified topology object.

        Parameters
        ----------
        topology : simtk.openmm.app.Topology object
            The topology of the system to construct.

        Returns
        -------
        system : openmm.System
            A system object generated from the topology
        """
        system = self._forcefield.createSystem(topology, **self._forcefield_kwargs)
        return system

    @property
    def ffxmls(self):
        return self._forcefield_xmls

    @property
    def forcefield(self):
        return self._forcefield
