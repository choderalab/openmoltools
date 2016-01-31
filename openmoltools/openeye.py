import shutil
import os
import mdtraj as md
from openmoltools.utils import import_, enter_temp_directory, create_ffxml_file
from openmoltools.amber import run_antechamber
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")


# Note: We recommend having every function return *copies* of input, to avoid headaches associated with in-place changes

def get_charges(molecule, max_confs=800, strictStereo=True, normalize=True, keep_confs=None):
    """Generate charges for an OpenEye OEMol molecule.

    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers.
        Omega will be used to generate max_confs conformations.
    max_confs : int, optional, default=800
        Max number of conformers to generate
    strictStereo : bool, optional, default=True
        If False, permits smiles strings with unspecified stereochemistry.
        See https://docs.eyesopen.com/omega/usage.html
    normalize : bool, optional, default=True
        If True, normalize the molecule by checking aromaticity, adding
        explicit hydrogens, and renaming by IUPAC name.
    keep_confs : int, optional, default=None
        If None, apply the charges to the provided conformation and return
        this conformation, unless no conformation is present.
        Otherwise, return some or all of the generated
        conformations. If -1, all generated conformations are returned.
        Otherwise, keep_confs = N will return an OEMol with up to N
        generated conformations.  Multiple conformations are still used to
        *determine* the charges.

    Returns
    -------
    charged_copy : OEMol
        A molecule with OpenEye's recommended AM1BCC charge selection scheme.

    Notes
    -----
    Roughly follows
    http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    """

    # If there is no geometry, return at least one conformation.
    if molecule.GetConfs() == 0:
        keep_confs = 1

    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oequacpac = import_("openeye.oequacpac")
    if not oequacpac.OEQuacPacIsLicensed(): raise(ImportError("Need License for oequacpac!"))

    if normalize:
        molecule = normalize_molecule(molecule)
    else:
        molecule = oechem.OEMol(molecule)

    charged_copy = generate_conformers(molecule, max_confs=max_confs, strictStereo=strictStereo)  # Generate up to max_confs conformers

    status = oequacpac.OEAssignPartialCharges(charged_copy, oequacpac.OECharges_AM1BCCSym)  # AM1BCCSym recommended by Chris Bayly to KAB+JDC, Oct. 20 2014.

    if not status:
        raise(RuntimeError("OEAssignPartialCharges returned error code %d" % status))

    #Determine conformations to return
    if keep_confs == None:
        #If returning original conformation
        original = molecule.GetCoords()
        #Delete conformers over 1
        for k, conf in enumerate( charged_copy.GetConfs() ):
            if k > 0:
                charged_copy.DeleteConf(conf)
        #Copy coordinates to single conformer
        charged_copy.SetCoords( original )
    elif keep_confs > 0:
        #Otherwise if a number is provided, return this many confs if available
        for k, conf in enumerate( charged_copy.GetConfs() ):
            if k > keep_confs - 1:
                charged_copy.DeleteConf(conf)
    elif keep_confs == -1:
        #If we want all conformations, continue
        pass
    else:
        #Not a valid option to keep_confs
        raise(ValueError('Not a valid option to keep_confs in get_charges.'))

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
        A normalized molecule with desired iupac name.

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

    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))

    molecule = oechem.OEMol()
    if not oechem.OEParseSmiles(molecule, smiles):
        raise ValueError("The supplied SMILES '%s' could not be parsed." % smiles)

    molecule = normalize_molecule(molecule)

    return molecule

def generate_conformers(molecule, max_confs=800, strictStereo=True, ewindow=15.0, rms_threshold=1.0, strictTypes = True):
    """Generate conformations for the supplied molecule

    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers
    max_confs : int, optional, default=800
        Max number of conformers to generate.  If None, use default OE Value.
    strictStereo : bool, optional, default=True
        If False, permits smiles strings with unspecified stereochemistry.
    strictTypes : bool, optional, default=True
        If True, requires that Omega have exact MMFF types for atoms in molecule; otherwise, allows the closest atom type of the same element to be used.

    Returns
    -------
    molcopy : OEMol
        A multi-conformer molecule with up to max_confs conformers.

    Notes
    -----
    Roughly follows
    http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html

    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oeomega = import_("openeye.oeomega")
    if not oeomega.OEOmegaIsLicensed(): raise(ImportError("Need License for OEOmega!"))

    molcopy = oechem.OEMol(molecule)
    omega = oeomega.OEOmega()

    # These parameters were chosen to match http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    omega.SetMaxConfs(max_confs)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(False)

    omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
    omega.SetEnergyWindow(ewindow)
    omega.SetRMSThreshold(rms_threshold)  # Word to the wise: skipping this step can lead to significantly different charges!

    omega.SetStrictStereo(strictStereo)
    omega.SetStrictAtomTypes(strictTypes)

    omega.SetIncludeInput(False)  # don't include input
    if max_confs is not None:
        omega.SetMaxConfs(max_confs)

    status = omega(molcopy)  # generate conformation
    if not status:
        raise(RuntimeError("omega returned error code %d" % status))


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


def molecule_to_mol2(molecule, tripos_mol2_filename=None, conformer=0, residue_name="MOL"):
    """Convert OE molecule to tripos mol2 file.

    Parameters
    ----------
    molecule : openeye.oechem.OEGraphMol
        The molecule to be converted.
    tripos_mol2_filename : str, optional, default=None
        Output filename.  If None, will create a filename similar to
        name.tripos.mol2, where name is the name of the OE molecule.
    conformer : int, optional, default=0
        Save this frame
    residue_name : str, optional, default="MOL"
        OpenEye writes mol2 files with <0> as the residue / ligand name.
        This chokes many mol2 parsers, so we replace it with a string of
        your choosing.

    Returns
    -------
    tripos_mol2_filename : str
        Filename of output tripos mol2 file

    """

    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for oechem!"))

    # Get molecule name.
    molecule_name = molecule.GetTitle()
    logger.debug(molecule_name)

    # Write molecule as Tripos mol2.
    if tripos_mol2_filename is None:
        tripos_mol2_filename = molecule_name + '.tripos.mol2'

    ofs = oechem.oemolostream(tripos_mol2_filename)
    ofs.SetFormat(oechem.OEFormat_MOL2H)
    for k, mol in enumerate(molecule.GetConfs()):
        if k == conformer:
            oechem.OEWriteMolecule(ofs, mol)

    ofs.close()

    # Replace <0> substructure names with valid text.
    infile = open(tripos_mol2_filename, 'r')
    lines = infile.readlines()
    infile.close()
    newlines = [line.replace('<0>', residue_name) for line in lines]
    outfile = open(tripos_mol2_filename, 'w')
    outfile.writelines(newlines)
    outfile.close()

    return molecule_name, tripos_mol2_filename


def oemols_to_ffxml(molecules, base_molecule_name="lig"):
    """Generate an OpenMM ffxml object and MDTraj trajectories from multiple OEMols

    Parameters
    ----------
    molecules : list(OEMole)
        Molecules WITH CHARGES.  Each can have multiple conformations.
        WILL GIVE UNDEFINED RESULTS IF NOT CHARGED.
    base_molecule_name : str, optional, default='lig'
        Base name of molecule to use inside parameter files.

    Returns
    -------
    trajectories : list(mdtraj.Trajectory)
        List of MDTraj Trajectories for molecule.  May contain multiple frames
    ffxml : StringIO
        StringIO representation of ffxml file.

    Notes
    -----
    We allow multiple different molecules at once so that they can all be
    included in a single ffxml file, which is currently the only recommended
    way to simulate multiple GAFF molecules in a single simulation.  For most
    applications, you will have just a single molecule:
    e.g. molecules = [my_oemol]
    The resulting ffxml StringIO object can be directly input to OpenMM e.g.
    `forcefield = app.ForceField(ffxml)`

    This will generate a lot of temporary files, so you may want to use
    utils.enter_temp_directory() to avoid clutter.
    """
    all_trajectories = []
    gaff_mol2_filenames = []
    frcmod_filenames = []

    print(os.getcwd())
    for i, molecule in enumerate(molecules):
        trajectories = []
        for j in range(molecule.NumConfs()):
            molecule_name = "%s-%d-%d" % (base_molecule_name, i, j)
            mol2_filename = "./%s.mol2" % molecule_name
            _unused = molecule_to_mol2(molecule, mol2_filename, conformer=j)
            gaff_mol2_filename, frcmod_filename = run_antechamber(molecule_name, mol2_filename, charge_method=None)  # It's redundant to run antechamber on each frame, fix me later.

            traj = md.load(gaff_mol2_filename)
            trajectories.append(traj)

            if j == 0:  # Only need 1 frame of forcefield files
                gaff_mol2_filenames.append(gaff_mol2_filename)
                frcmod_filenames.append(frcmod_filename)

        # Create a trajectory with all frames of the current molecule
        traj = trajectories[0].join(trajectories[1:])
        all_trajectories.append(traj)

    ffxml = create_ffxml_file(gaff_mol2_filenames, frcmod_filenames, override_mol2_residue_name=base_molecule_name)

    return all_trajectories, ffxml


def smiles_to_antechamber(smiles_string, gaff_mol2_filename, frcmod_filename, residue_name="MOL", strictStereo=False):
    """Build a molecule from a smiles string and run antechamber,
    generating GAFF mol2 and frcmod files from a smiles string.  Charges
    will be generated using the OpenEye QuacPac AM1-BCC implementation.

    Parameters
    ----------
    smiles_string : str
        Smiles string of molecule to construct and charge
    gaff_mol2_filename : str
        Filename of mol2 file output of antechamber, with charges
        created from openeye
    frcmod_filename : str
        Filename of frcmod file output of antechamber.  Most likely
        this file will be almost empty, at least for typical molecules.
    residue_name : str, optional, default="MOL"
        OpenEye writes mol2 files with <0> as the residue / ligand name.
        This chokes many mol2 parsers, so we replace it with a string of
        your choosing.  This might be useful for downstream applications
        if the residue names are required to be unique.
    strictStereo : bool, optional, default=False
        If False, permits smiles strings with unspecified stereochemistry.
        See https://docs.eyesopen.com/omega/usage.html
    """
    oechem = import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for oechem!"))

    # Get the absolute path so we can find these filenames from inside a temporary directory.
    gaff_mol2_filename = os.path.abspath(gaff_mol2_filename)
    frcmod_filename = os.path.abspath(frcmod_filename)

    m = smiles_to_oemol(smiles_string)
    m = get_charges(m, strictStereo=strictStereo, keep_confs=1)

    with enter_temp_directory():  # Avoid dumping 50 antechamber files in local directory.
        _unused = molecule_to_mol2(m, "./tmp.mol2", residue_name=residue_name)
        net_charge = oechem.OENetCharge(m)
        tmp_gaff_mol2_filename, tmp_frcmod_filename = run_antechamber("tmp", "./tmp.mol2", charge_method=None, net_charge=net_charge)  # USE OE AM1BCC charges!
        shutil.copy(tmp_gaff_mol2_filename, gaff_mol2_filename)
        shutil.copy(tmp_frcmod_filename, frcmod_filename)
