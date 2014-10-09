from utils import import_

def normalize_molecule(molecule):
    """Normalize a copy of the molecule by checking aromaticity, adding explicit hydrogens, and renaming by IUPAC name.

    ARGUMENTS
    molecule (OEMol) - the molecule to be normalized.

    EXAMPLES
    # read a partial molecule and normalize it
    molecule = readMolecule('molecule.sdf')
    normalizeMolecule(molecule)
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))   
    oeiupac = import_("openeye.oeiupac")
    if not oeiupac.OEIUPACIsLicensed(): raise(ImportError("Need License for OEOmega!"))    
    molcopy = oechem.OEMol(molecule)
   
    # Assign aromaticity.
    oechem.OEAssignAromaticFlags(molcopy, oechem.OEAroModelOpenEye)

    # Add hydrogens.
    oechem.OEAddExplicitHydrogens(molcopy)

    # Set title to IUPAC name.
    name = oeiupac.OECreateIUPACName(molcopy)
    molcopy.SetTitle(name)

    return molcopy


def iupac_to_oemol(iupac_name):
    """Create a OEMolBuilder from a iupac name.
    
    Parameters
    ----------
    iupac_name : str
        IUPAC name of desired molecule.

    Examples
    --------
    >>> molecule = iupac_to_oemol("trans-2-fluoro-3-methylpent-2-ene")
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oeiupac = import_("openeye.oeiupac")
    if not oeiupac.OEIUPACIsLicensed(): raise(ImportError("Need License for OEOmega!"))

    # Create an OEMol molecule from IUPAC name.
    molecule = oechem.OEMol()  # create a molecule

    # Populate the MoleCule from the IUPAC name
    if not oeiupac.OEParseIUPACName(molecule, iupac_name):
        raise ValueError("The supplied IUPAC name '%s' could not be parsed." % iupac_name)

    normalize_molecule(molecule)

    return molecule

def smiles_to_oemol(smiles):
    """Create a OEMolBuilder from a smiles string.
    
    Parameters
    ----------
    smiles : str
        SMILES representation of desired molecule.
    Examples
    --------
    >>> molecule = smiles_to_oemol("c1ccncc1")
    """        
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    
    molecule = oechem.OEMol()
    if not oechem.OEParseSmiles(molecule, smiles):
        raise ValueError("The supplied SMILES '%s' could not be parsed." % smiles)

    normalize_molecule(molecule)

    return molecule

def generate_conformers(molecule, max_conformers, strictStereo=True):
    """Generate conformations for the supplied molecule

    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers
    max_conformers : int
        Max number of conformers to generate
    strictStereo : bool, optional, default=True
        Adhere to strict specification of stereo isomer

    Examples
    --------
    >>> molecule = iupac_to_oemol("trans-2-fluoro-3-methylpent-2-ene")
    >>> molecule_w_conformers = generate_conformers(molecule, 5, strictStereo=True)
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oeomega = import_("openeye.oeomega")
    if not oeomega.OEOmegaIsLicensed(): raise(ImportError("Need License for OEOmega!"))
    molcopy = oechem.OEMol(molecule)
    omega = oeomega.OEOmega()
    omega.SetStrictStereo(strictStereo)

    #if strictTyping != None:
    #    omega.SetStrictAtomTypes(strictTyping)
   
    omega.SetIncludeInput(False)  # don't include input
    omega.SetMaxConfs(max_conformers)
    omega(molcopy)  # generate conformation

    return molcopy


def generate_am1bcc_charges(molecule):
    """Use OEQuacPac to generate am1bcc charges for a molecule."""
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oequacpac = import_("openeye.oequacpac")
    if not oequacpac.OEQuacPacIsLicensed(): raise(ImportError("Need License for oequacpac!"))

    charged_copy = oechem.OEMol(molecule)
    oequacpac.OEAssignPartialCharges(charged_copy, oequacpac.OECharges_AM1BCC)
    return charged_copy


def generate_conformer_charges_am1bcc(molecule_w_conformers):
    """Generate charges for each conformer in a given molecule

    Parameters
    ----------
    molecule_w_conformers: OEMol
        A molecule with several 3D conformers

    Examples
    --------
    >>> molecule = iupac_to_oemol("trans-2-fluoro-3-methylpent-2-ene")
    >>> molecule_w_conformers = generate_conformers(molecule, 5, strictStereo=True)
    >>> charged_conformers = generate_conformer_charges_am1bcc(molecule_w_conformers)
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))

    molcopy = oechem.OEMol(molecule_w_conformers)
    charged_conformers = list()
    for conformer_index, conformation in enumerate(molcopy.GetConfs()):
        charged_molecule = generate_am1bcc_charges(conformation)
        charged_conformers.append(charged_molecule)

    return charged_conformers


def get_oeszybki_minimizer(use_charges=True):
    """Set up a minimizer with OESzybki"""
    oeszybki = import_("openeye.oeszybki")
    if not oeszybki.OESzybkiIsLicensed(): raise(ImportError("Need License for oeszybki!"))

    szybki = oeszybki.OESzybki()

    szybki.SetRunType(oeszybki.OERunType_CartesiansOpt)  # Set the runtype to minimizer

    if use_charges:
        szybki.SetUseCurrentCharges(True)  # use partial charges for minimization
    else:
        szybki.SetUseCurrentCharges(False) # use mmff94 charges for minimization

    return szybki


def minimize_conformers(conformer_list, use_charges=True):
    """Minimize the conformations in a list of conformers.

    Parameters
    ----------
    conformer_list: list of OEMol
        A list of conformers that need to be minimized
    use_charges: bool, optional, default=True
        Use the molecules current charges for minimization

    Examples
    --------
    >>> molecule = iupac_to_oemol("trans-2-fluoro-3-methylpent-2-ene")
    >>> molecule_w_conformers = generate_conformers(molecule, 5, strictStereo=True)
    >>> charged_conformers = generate_conformer_charges_am1bcc(molecule_w_conformers)
    >>> minimized_conformers = minimize_conformers(charged_conformers, use_charges=True)
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for oechem!"))

    szybki = get_oeszybki_minimizer(use_charges)

    minimized_conformers = list()
    conformer_energies = list()
    for conformer in conformer_list:
        conformer_copy = oechem.OEMol(conformer)
        szybki_result = szybki(conformer_copy)
        minimized_conformers.append(conformer_copy)
        lastresult = None
        for lastresult in szybki_result:
            pass
        conformer_energies.append(lastresult.GetTotalEnergy()) # Unit is kcal / mol ?
    return minimized_conformers, conformer_energies

def absolute_charges(molecule):
    """Set partial charges to the absolute value."""
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for oechem!"))
    molcopy = oechem.OEMol(molecule)

    for atom in molcopy.GetAtoms():
            atom.SetPartialCharge(abs(atom.GetPartialCharge()))

    return molcopy

def generate_extended_conformers(molecule, max_conformers, strictStereo=True):
    """Generate extended conformations for the supplied molecule using conformers minimized by absolute charges.

    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers
    max_conformers : int
        Max number of conformers to generate
    strictStereo : bool, optional, default=True
        Adhere to strict specification of stereo isomer

    Examples
    --------
    >>> molecule = iupac_to_oemol("7-azaniumylheptanoate")
    >>> extended_conformers = generate_extended_conformers(molecule, 50, strictStereo=True)
    >>> z = zip(*extended_conformers)
    >>> for p in z: print p[1]- p[2] # debug change in energy after minimization with abs charge
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for oechem!"))

    molcopy = oechem.OEMol(molecule)

    molecule_w_conformers = generate_conformers(molcopy, max_conformers, strictStereo)
    charged_conformers = generate_conformer_charges_am1bcc(molecule_w_conformers)
    minimized_charged_conformers, min_energies = minimize_conformers(charged_conformers, use_charges=True)

    absolute_charged_conformers = list()
    print get_mol_atoms_type_charge(minimized_charged_conformers[0])  # debug
    for charged_conf in minimized_charged_conformers:
        absolute_charged_conformers.append(absolute_charges(charged_conf))

    print get_mol_atoms_type_charge(absolute_charged_conformers[0]) #debug
    extended_conformers, ext_energies = minimize_conformers(absolute_charged_conformers, use_charges=True)

    return extended_conformers, ext_energies, min_energies

def get_mol_atoms_type_charge(oemol):
    """
    Print type and partial charge for every atom in OEMol
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for oechem!"))
    molcopy = oechem.OEMol(oemol)
    molrepr = str()
    for at in molcopy.GetAtoms():
        tp = at.GetType()
        pc = at.GetPartialCharge()
        molrepr += "%s %f \n" % (tp,pc)
    return molrepr