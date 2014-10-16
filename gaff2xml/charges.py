from utils import import_

# Note: We recommend having every function return *copies* of input, to avoid headaches associated with in-place changes

def get_charges(molecule, max_confs=200):
    """Generate charges for an OpenEye OEMol molecule.

    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers.  If molecule.NumConfs() <= 1
        and max_confs > 1, then Omega will be used to generate max_confs
        conformations.  
    max_confs : int
        Max number of conformers to generate
    
    Returns
    -------
    charged_copy : OEMol
        A molecule with OpenEye's recommended AM1BCC charge selection scheme.

    Notes
    -----
    See http://docs.eyesopen.com/toolkits/quacpac/java/molchargetheory.html
    for details.
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oequacpac = import_("openeye.oequacpac")
    if not oequacpac.OEQuacPacIsLicensed(): raise(ImportError("Need License for oequacpac!"))
    
    molecule = normalize_molecule(molecule)
    
    if molecule.NumConfs() <= 1 and max_confs > 1:
        charged_copy = generate_conformers(molecule, max_confs=max_confs)  # Generate up to max_confs conformers
    else:
        charged_copy = molecule  # Just charge the input molecule
    
    oequacpac.OEAssignPartialCharges(charged_copy, oequacpac.OECharges_AM1BCCSymSPt)  # AM1BCCSymSPt recommended by Chris Bayly to KAB, Oct. 15 2014.

    return charged_copy


def normalize_molecule(molecule):
    """Normalize a copy of the molecule by checking aromaticity, adding explicit hydrogens, and renaming by IUPAC name.

    Parameters
    ----------
    molecule : OEMol
        the molecule to be normalized.
    
    Returns
    -------
    molcopy : OEMol
        A (copied) version of the normalized molecule

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
    
    # Check for any missing atom names, if found reassign all of them.
    if any([atom.GetName() == '' for atom in molcopy.GetAtoms()]):
        oechem.OETriposAtomNames(molcopy)
    
    return molcopy


def iupac_to_oemol(iupac_name):
    """Create a OEMolBuilder from a iupac name.
    
    Parameters
    ----------
    iupac_name : str
        IUPAC name of desired molecule.

    Returns
    -------
    molecule : OEMol
        A normalized molecule with desired iupac name

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

    molecule = normalize_molecule(molecule)

    return molecule

def smiles_to_oemol(smiles):
    """Create a OEMolBuilder from a smiles string.
    
    Parameters
    ----------
    smiles : str
        SMILES representation of desired molecule.

    Returns
    -------
    molecule : OEMol
        A normalized molecule with desired smiles string.
    
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

def generate_conformers(molecule, max_confs=100, strictStereo=True):
    """Generate conformations for the supplied molecule

    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers
    max_confs : int
        Max number of conformers to generate
    strictStereo : bool, optional, default=True
        Adhere to strict specification of stereo isomer

    Returns
    -------
    molcopy : OEMol
        A multi-conformer molecule with up to max_confs conformers.

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
   
    omega.SetIncludeInput(False)  # don't include input
    omega.SetMaxConfs(max_confs)
    omega(molcopy)  # generate conformation

    return molcopy


def get_names_to_charges(molecule):
    """Return a dictionary of atom names and partial charges, as well as a string representation.

    Parameters
    ----------
    molecule : OEMol
        Molecule for which to grab charges

    Returns
    -------
    data : dictionary
        A dictinoary whose (key, val) pairs are the atom names and partial
        charges, respectively.
    molrepr : str
        A string representation of data
    """
    
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for oechem!"))
    molcopy = oechem.OEMol(molecule)
    molrepr = ""
    data = {}
    for atom in molcopy.GetAtoms():
        name = atom.GetName()
        charge = atom.GetPartialCharge()
        data[name] = charge
        molrepr += "%s %f \n" % (name, charge)
    return data, molrepr
