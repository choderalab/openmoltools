from utils import import_
from copy import deepcopy

def normalize_molecule(molecule):
    """Normalize a copy of the molecule by checking aromaticity, adding explicit hydrogens, and renaming by IUPAC name.

    ARGUMENTS
    molecule (OEMol) - the molecule to be normalized.

    EXAMPLES
    # read a partial molecule and normalize it
    molecule = readMolecule('molecule.sdf')
    normalizeMolecule(molecule)
    """
    molcopy = deepcopy(molecule)
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))   
    oeiupac = import_("openeye.oeiupac")
    if not oeiupac.OEIUPACIsLicensed(): raise(ImportError("Need License for OEOmega!"))    
   
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
    molecule = iupac_to_oemol("trans-2-fluoro-3-methylpent-2-ene")
    molecule = generate_conformers(molecule, 5, strictStereo=True)
    """
    molcopy = deepcopy(molecule)
    oeomega = import_("openeye.oeomega")
    if not oeomega.OEOmegaIsLicensed(): raise(ImportError("Need License for OEOmega!"))
    
    omega = oeomega.OEOmega()
    omega.SetStrictStereo(strictStereo)

    #if strictTyping != None:
    #    omega.SetStrictAtomTypes(strictTyping)
   
    omega.SetIncludeInput(False)  # don't include input
    omega.SetMaxConfs(max_conformers)
    omega(molcopy)  # generate conformation

    return molcopy


def find_conformer_for_charges(molecule, verbose = False):
    """Pick an optimal extended conformer to be used for AM1BCC charging in antechamber.

    Parameters
    ----------

    molecule : OEMol
        molecule for which charges are to be assigned, 
    verbose : bool, optional, default=False
        If True, information about the current calculation is printed.

    Returns
    -------
     output_molecule : OEMol
        Optimal conformation for further preparation via antechamber.

    Notes
    -----
    Ask us about why we're doing this.
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oeszybki = import_("openeye.oeszybki")
    if not oeszybki.OESzybkiIsLicensed(): raise(ImportError("Need License for oeszybki!"))
    oequacpac = import_("openeye.oequacpac")
    if not oequacpac.OEQuacPacIsLicensed(): raise(ImportError("Need License for oequacpac!"))

    oequacpac
    
    #Check that molecule has atom names; if not we need to assign them
    assignNames = False
    for atom in molecule.GetAtoms():
        if atom.GetName()=='':
            assignNames = True #In this case we are missing an atom name and will need to assign
    if assignNames:
        if verbose: print("Assigning TRIPOS names to atoms")
        oechem.OETriposAtomNames(molecule)

    # Check input pameters.
    supported_charge_models  = ['am1bcc']
    if not (charge_model in supported_charge_models):
        raise "Charge model %(charge_model)s not in supported set of %(supported_charge_models)s" % vars()

    # Set up storage for partial charges.
    partial_charges = dict()
    for atom in molecule.GetAtoms():
        name = atom.GetName()
        partial_charges[name] = 0.0

    for conformer_index, conformation in enumerate(molecule.GetConfs()):
        if verbose: print("Assigning partial charges to conformer %d / %d" % (conformer_index, molecule.NumConfs()))
        # Assign partial charges to a copy of the molecule.
        if verbose: print("assignPartialCharges: determining partial charges...")
        charged_molecule = oechem.OEMol(conformation)   
        oequacpac.OEAssignPartialCharges(charged_molecule, oequacpac.OECharges_AM1BCC)         
      
        #   Set partial charges to absolute value.
        for atom in charged_molecule.GetAtoms():
            atom.SetPartialCharge(abs(atom.GetPartialCharge()))
        # Minimize in Cartesian space to splay out substructures.
        szybki = oeszybki.OESzybki() # create an instance of OESzybki
        szybki.SetRunType(oeszybki.OERunType_CartesiansOpt) # set minimization         
        szybki.SetUseCurrentCharges(True) # use charges for minimization
        results = szybki(charged_molecule)

        if verbose: print("assignPartialCharges: redetermining partial charges...")
        oequacpac.OEAssignPartialCharges(charged_molecule, oequacpac.OECharges_AM1BCC)         
         
    # Compute and store average partial charges in a copy of the original molecule.
    charged_molecule = oechem.OEMol(molecule)
    for atom in charged_molecule.GetAtoms():
        name = atom.GetName()
        #atom.SetPartialCharge(partial_charges[name] / nconformers)

    # Return the charged molecule
    return charged_molecule
