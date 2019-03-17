import unittest
import pytest
import tempfile
import os, sys
import numpy as np
from openmoltools import utils, forcefield_generators
from simtk import openmm, unit
from simtk.openmm.app import ForceField, NoCutoff, CutoffPeriodic, OBC2
if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO

ISTRAVIS = os.environ.get('TRAVIS', None) == 'true'

################################################################################
# Suppress matplotlib logging
################################################################################

import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

################################################################################
# OpenEye imports
################################################################################

try:
    oechem = utils.import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oequacpac = utils.import_("openeye.oequacpac")
    if not oequacpac.OEQuacPacIsLicensed(): raise(ImportError("Need License for oequacpac!"))
    oeiupac = utils.import_("openeye.oeiupac")
    if not oeiupac.OEIUPACIsLicensed(): raise(ImportError("Need License for OEOmega!"))
    oeomega = utils.import_("openeye.oeomega")
    if not oeomega.OEOmegaIsLicensed(): raise(ImportError("Need License for OEOmega!"))
    HAVE_OE = True
    openeye_exception_message = str()
except Exception as e:
    HAVE_OE = False
    openeye_exception_message = str(e)

################################################################################
# Utility methods
################################################################################

IUPAC_molecule_names = ['naproxen', 'aspirin', 'imatinib', 'bosutinib', 'dibenzyl ketone']
def createOEMolFromIUPAC(iupac_name='bosutinib'):
    from openeye import oechem, oeiupac, oeomega

    # Create molecule.
    mol = oechem.OEMol()
    oeiupac.OEParseIUPACName(mol, iupac_name)
    mol.SetTitle(iupac_name)

    # Assign aromaticity and hydrogens.
    oechem.OEAssignAromaticFlags(mol, oechem.OEAroModelOpenEye)
    oechem.OEAddExplicitHydrogens(mol)

    # Create atom names.
    oechem.OETriposAtomNames(mol)

    # Assign geometry
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(False)
    omega.SetStrictStereo(True)
    omega(mol)

    return mol

def extractPositionsFromOEMOL(molecule):
    positions = unit.Quantity(np.zeros([molecule.NumAtoms(), 3], np.float32), unit.angstroms)
    coords = molecule.GetCoords()
    for index in range(molecule.NumAtoms()):
        positions[index,:] = unit.Quantity(coords[index], unit.angstroms)
    return positions

def check_potential_is_finite(system, positions):
    """
    Check that the potential energy is finite.

    Parameters
    ----------
    system : simtk.openmm.Syste
        System object
    positions : simtk.unit.Quantity of (natoms,3) with units compatible with angstroms
        positions

    Returns
    -------
    potential : simtk.unit.Quantity
        The potential energy.

    Raises an Exception if the potential energy is not finite.

    """
    integrator = openmm.VerletIntegrator(1 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    potential = context.getState(getEnergy=True).getPotentialEnergy()
    del context, integrator

    if np.isnan(potential / unit.kilojoules_per_mole):
        raise Exception("Potential energy is infinite.")

    return potential

#@unittest.skipIf(os.getenv('AMBERHOME') is None, 'Cannot test w/out Amber')
def testWriteXMLParametersGAFF():
    """ Test writing XML parameters loaded from Amber GAFF parameter files """

    # Generate ffxml file contents for parmchk-generated frcmod output.
    gaff_dat_filename = utils.get_data_filename("parameters/gaff.dat")
    leaprc = StringIO("parm = loadamberparams %s" % gaff_dat_filename)
    import parmed
    params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    citations = """\
Wang, J., Wang, W., Kollman P. A.; Case, D. A. "Automatic atom type and bond type perception in molecular mechanical calculations". Journal of Molecular Graphics and Modelling , 25, 2006, 247260.
Wang, J., Wolf, R. M.; Caldwell, J. W.;Kollman, P. A.; Case, D. A. "Development and testing of a general AMBER force field". Journal of Computational Chemistry, 25, 2004, 1157-1174.
"""
    ffxml = str()
    provenance=dict(OriginalFile='gaff.dat', Reference=citations)
    outfile = open('gaff.xml', 'w')
    params.write(outfile, provenance=provenance)
    outfile.close()

@pytest.mark.skipif(not HAVE_OE, reason="Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_generate_ffxml_from_molecules():
    """
    Test generation of single ffxml file from a list of molecules
    """
    # Create a test set of molecules.
    molecules = [ createOEMolFromIUPAC(name) for name in IUPAC_molecule_names ]
    # Create an ffxml file.
    from openmoltools.forcefield_generators import generateForceFieldFromMolecules
    ffxml = generateForceFieldFromMolecules(molecules)
    # Create a ForceField.
    gaff_xml_filename = utils.get_data_filename("parameters/gaff.xml")
    forcefield = ForceField(gaff_xml_filename)
    try:
        forcefield.loadFile(StringIO(ffxml))
    except Exception as e:
        msg  = str(e)
        msg += "ffxml contents:\n"
        for (index, line) in enumerate(ffxml.split('\n')):
            msg += 'line %8d : %s\n' % (index, line)
        raise Exception(msg)

    # Parameterize the molecules.
    from openmoltools.forcefield_generators import generateTopologyFromOEMol
    for molecule in molecules:
        # Create topology from molecule.
        topology = generateTopologyFromOEMol(molecule)
        # Create system with forcefield.
        system = forcefield.createSystem(topology)
        # Check potential is finite.
        positions = extractPositionsFromOEMOL(molecule)
        check_potential_is_finite(system, positions)

@pytest.mark.skipif(not HAVE_OE, reason="Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_generate_gaff2_ffxml_from_molecules():
    """
    Test generation of single ffxml file from a list of molecules, using the gaff2 option.
    """
    # Create a test set of molecules.
    molecules = [ createOEMolFromIUPAC(name) for name in IUPAC_molecule_names ]
    # Create an ffxml file.
    from openmoltools.forcefield_generators import generateForceFieldFromMolecules
    ffxml = generateForceFieldFromMolecules(molecules, gaff_version='gaff2')
    # Create a ForceField.
    gaff_xml_filename = utils.get_data_filename("parameters/gaff2.xml")
    forcefield = ForceField(gaff_xml_filename)
    try:
        forcefield.loadFile(StringIO(ffxml))
    except Exception as e:
        msg  = str(e)
        msg += "ffxml contents:\n"
        for (index, line) in enumerate(ffxml.split('\n')):
            msg += 'line %8d : %s\n' % (index, line)
        raise Exception(msg)

    # Parameterize the molecules.
    from openmoltools.forcefield_generators import generateTopologyFromOEMol
    for molecule in molecules:
        # Create topology from molecule.
        topology = generateTopologyFromOEMol(molecule)
        # Create system with forcefield.
        system = forcefield.createSystem(topology)
        # Check potential is finite.
        positions = extractPositionsFromOEMOL(molecule)
        check_potential_is_finite(system, positions)

@pytest.mark.skipif(not HAVE_OE, reason="Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_topology_molecules_round_trip():
    """
    Test round-trips between OEMol and Topology
    """
    # Create a test set of molecules.
    molecules = [ createOEMolFromIUPAC(name) for name in IUPAC_molecule_names ]
    # Test round-trips.
    from openmoltools.forcefield_generators import generateTopologyFromOEMol, generateOEMolFromTopologyResidue
    for molecule in molecules:
        # Create topology from molecule.
        topology = generateTopologyFromOEMol(molecule)
        # Create molecule from topology.
        residues = [residue for residue in topology.residues()]
        molecule2 = generateOEMolFromTopologyResidue(residues[0])
        # Create topology form molecule.
        topology2 = generateTopologyFromOEMol(molecule2)
        # Create molecule from topology with geometry.
        residues2 = [residue for residue in topology2.residues()]
        molecule3 = generateOEMolFromTopologyResidue(residues2[0], geometry=True)
        # Create molecule from topology with Tripos atom names
        molecule4 = generateOEMolFromTopologyResidue(residues2[0], tripos_atom_names=True)

class TestForceFieldGenerators(unittest.TestCase):
    @pytest.mark.skipif(not HAVE_OE, reason="Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
    def test_generate_Topology_and_OEMol(self):
        """
        Test round-trip from OEMol >> Topology >> OEMol
        """
        from openmoltools.forcefield_generators import generateTopologyFromOEMol, generateOEMolFromTopologyResidue
        from openeye import oechem, oeiupac
        for molecule_name in IUPAC_molecule_names:
            molecule1 = createOEMolFromIUPAC(molecule_name)

            # Generate Topology from OEMol
            topology = generateTopologyFromOEMol(molecule1)
            # Check resulting Topology.
            residues = [ residue for residue in topology.residues() ]
            self.assertEqual(len(residues), 1)
            self.assertEqual(residues[0].name, molecule1.GetTitle())
            for (top_atom, mol_atom) in zip(topology.atoms(), molecule1.GetAtoms()):
                self.assertEqual(top_atom.name, mol_atom.GetName())
            for (top_bond, mol_bond) in zip(topology.bonds(), molecule1.GetBonds()):
                self.assertEqual(top_bond[0].name, mol_bond.GetBgn().GetName())
                self.assertEqual(top_bond[1].name, mol_bond.GetEnd().GetName())

            # Generate OEMol from Topology
            molecule2 = generateOEMolFromTopologyResidue(residues[0])
            # Check resulting molecule.
            self.assertEqual(molecule1.GetTitle(), molecule2.GetTitle())
            for (atom1, atom2) in zip(molecule1.GetAtoms(), molecule2.GetAtoms()):
                self.assertEqual(atom1.GetName(), atom2.GetName())
                self.assertEqual(atom1.GetAtomicNum(), atom2.GetAtomicNum())
            for (bond1, bond2) in zip(molecule1.GetBonds(), molecule2.GetBonds()):
                self.assertEqual(bond1.GetBgn().GetName(), bond2.GetBgn().GetName())
                self.assertEqual(bond1.GetEnd().GetName(), bond2.GetEnd().GetName())

@pytest.mark.skipif(not HAVE_OE, reason="Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_generateResidueTemplate():
    """
    Test GAFF residue template generation from OEMol molecules.
    """
    from openeye import oechem, oeiupac

    from pkg_resources import resource_filename
    gaff_xml_filename = utils.get_data_filename("parameters/gaff.xml")

    # Test independent ForceField instances.
    for molecule_name in IUPAC_molecule_names:
        mol = createOEMolFromIUPAC(molecule_name)
        # Generate an ffxml residue template.
        from openmoltools.forcefield_generators import generateResidueTemplate
        [template, ffxml] = generateResidueTemplate(mol)
        # Create a ForceField object.
        forcefield = ForceField(gaff_xml_filename)
        # Add the additional parameters and template to the forcefield.
        forcefield.registerResidueTemplate(template)
        forcefield.loadFile(StringIO(ffxml))
        # Create a Topology from the molecule.
        from openmoltools.forcefield_generators import generateTopologyFromOEMol
        topology = generateTopologyFromOEMol(mol)
        # Parameterize system.
        system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)
        # Check potential is finite.
        positions = extractPositionsFromOEMOL(mol)
        check_potential_is_finite(system, positions)

    # Test adding multiple molecules to a single ForceField instance.
    forcefield = ForceField(gaff_xml_filename)
    for molecule_name in IUPAC_molecule_names:
        mol = createOEMolFromIUPAC(molecule_name)
        # Generate an ffxml residue template.
        from openmoltools.forcefield_generators import generateResidueTemplate
        [template, ffxml] = generateResidueTemplate(mol)
        # Add the additional parameters and template to the forcefield.
        forcefield.registerResidueTemplate(template)
        forcefield.loadFile(StringIO(ffxml))
        # Create a Topology from the molecule.
        from openmoltools.forcefield_generators import generateTopologyFromOEMol
        topology = generateTopologyFromOEMol(mol)
        # Parameterize system.
        system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)
        # Check potential is finite.
        positions = extractPositionsFromOEMOL(mol)
        check_potential_is_finite(system, positions)


@pytest.mark.skipif(not HAVE_OE, reason="Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_generateResidueTemplate_gaff2():
    """
    Test GAFF2 residue template generation from OEMol molecules.
    """
    from openeye import oechem, oeiupac

    from pkg_resources import resource_filename
    gaff_xml_filename = utils.get_data_filename("parameters/gaff2.xml")

    # Test independent ForceField instances.
    for molecule_name in IUPAC_molecule_names:
        mol = createOEMolFromIUPAC(molecule_name)
        # Generate an ffxml residue template.
        from openmoltools.forcefield_generators import generateResidueTemplate
        [template, ffxml] = generateResidueTemplate(mol,gaff_version='gaff2')
        # Create a ForceField object.
        forcefield = ForceField(gaff_xml_filename)
        # Add the additional parameters and template to the forcefield.
        forcefield.registerResidueTemplate(template)
        forcefield.loadFile(StringIO(ffxml))
        # Create a Topology from the molecule.
        from openmoltools.forcefield_generators import generateTopologyFromOEMol
        topology = generateTopologyFromOEMol(mol)
        # Parameterize system.
        system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)
        # Check potential is finite.
        positions = extractPositionsFromOEMOL(mol)
        check_potential_is_finite(system, positions)

    # Test adding multiple molecules to a single ForceField instance.
    forcefield = ForceField(gaff_xml_filename)
    for molecule_name in IUPAC_molecule_names:
        mol = createOEMolFromIUPAC(molecule_name)
        # Generate an ffxml residue template.
        from openmoltools.forcefield_generators import generateResidueTemplate
        [template, ffxml] = generateResidueTemplate(mol, gaff_version='gaff2')
        # Add the additional parameters and template to the forcefield.
        forcefield.registerResidueTemplate(template)
        forcefield.loadFile(StringIO(ffxml))
        # Create a Topology from the molecule.
        from openmoltools.forcefield_generators import generateTopologyFromOEMol
        topology = generateTopologyFromOEMol(mol)
        # Parameterize system.
        system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)
        # Check potential is finite.
        positions = extractPositionsFromOEMOL(mol)
        check_potential_is_finite(system, positions)


def check_energy_components_vs_prmtop(prmtop=None, inpcrd=None, system=None, MAX_ALLOWED_DEVIATION=5.0):
    """
    """
    import parmed as pmd
    structure = pmd.load_file(prmtop, inpcrd)
    prmtop_components = dict(pmd.openmm.energy_decomposition_system(structure, structure.createSystem(nonbondedMethod=NoCutoff)))
    system_components = dict(pmd.openmm.energy_decomposition_system(structure, system))

    msg  = "\n"
    msg += "Energy components:\n"
    test_pass = True
    msg += "%20s %12s %12s : %12s\n" % ('component', 'prmtop (kcal/mol)', 'system (kcal/mol)', 'deviation')
    for key in prmtop_components:
        e1 = prmtop_components[key]
        e2 = system_components[key]
        deviation = abs(e1-e2)
        if (deviation > MAX_ALLOWED_DEVIATION):
            test_pass = False
        msg += "%20s %20.6f %20.6f : %20.6f\n" % (key, e1, e2, deviation)

    if not test_pass:
        msg += "Maximum allowed deviation (%f) exceeded.\n" % MAX_ALLOWED_DEVIATION
        #raise Exception(msg) # TODO: Re-enable when we have force tag merging sorted out in simtk.openmm.app.ForceField
        print(msg) # DEBUG

@pytest.mark.skipif(not HAVE_OE, reason="Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_gaffResidueTemplateGenerator():
    """
    Test the GAFF residue template generator.
    """

    #
    # Test where we generate parameters for only a ligand.
    #

    # Load the PDB file.
    from simtk.openmm.app import PDBFile
    pdb_filename = utils.get_data_filename("chemicals/imatinib/imatinib.pdb")
    pdb = PDBFile(pdb_filename)
    # Create a ForceField object.
    gaff_xml_filename = utils.get_data_filename("parameters/gaff.xml")
    forcefield = ForceField(gaff_xml_filename)
    # Add the residue template generator.
    from openmoltools.forcefield_generators import gaffTemplateGenerator
    forcefield.registerTemplateGenerator(gaffTemplateGenerator)
    # Parameterize system.
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    # Check potential is finite.
    check_potential_is_finite(system, pdb.positions)
    # Check energy matches prmtop route.
    check_energy_components_vs_prmtop(
        prmtop=utils.get_data_filename('chemicals/imatinib/imatinib.prmtop'),
        inpcrd=utils.get_data_filename('chemicals/imatinib/imatinib.inpcrd'),
        system=system)

    #
    # Test where we generate parameters for only a ligand in a protein.
    #

    # Load the PDB file.
    from simtk.openmm.app import PDBFile
    pdb_filename = utils.get_data_filename("chemicals/proteins/T4-lysozyme-L99A-p-xylene-implicit.pdb")
    pdb = PDBFile(pdb_filename)
    # Create a ForceField object.
    gaff_xml_filename = utils.get_data_filename("parameters/gaff.xml")
    forcefield = ForceField('amber99sb.xml', gaff_xml_filename)
    # Add the residue template generator.
    from openmoltools.forcefield_generators import gaffTemplateGenerator
    forcefield.registerTemplateGenerator(gaffTemplateGenerator)
    # Parameterize system.
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    # Check potential is finite.
    check_potential_is_finite(system, pdb.positions)

@pytest.mark.skipif(not HAVE_OE, reason="Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_atom_topology_index():
    """
    Make sure that generateOEMolFromTopologyResidue adds the topology_index data
    """
    # Create a test set of molecules.
    molecules = [ createOEMolFromIUPAC(name) for name in IUPAC_molecule_names ]
    from openmoltools.forcefield_generators import generateTopologyFromOEMol, generateOEMolFromTopologyResidue
    topologies = [generateTopologyFromOEMol(molecule) for molecule in molecules]
    for topology in topologies:
        residue = list(topology.residues())[0] #there is only one residue
        regenerated_mol = generateOEMolFromTopologyResidue(residue)
        for i, top_atom in enumerate(topology.atoms()):
            oeatom = regenerated_mol.GetAtom(oechem.OEHasAtomIdx(top_atom.index))
            assert oeatom.GetData("topology_index")==top_atom.index

@pytest.mark.skipif(not HAVE_OE, reason="Cannot test without OpenEye tools")
def check_system_generator(ffxmls, forcefield_kwargs, system_name, **kwargs):
    """
    Check SystemGenerator on a specific topology.
    """
    import openmmtools
    try:
        constructor = getattr(openmmtools.testsystems, system_name)
        testsystem = constructor()
        topology = testsystem.topology
    except AttributeError:
        molecule = createOEMolFromIUPAC(system_name)
        topology = forcefield_generators.generateTopologyFromOEMol(molecule)
    system_generator = forcefield_generators.SystemGenerator(ffxmls,
                                                     forcefield_kwargs=forcefield_kwargs, **kwargs)
    system_generator.createSystem(topology)

def test_system_generator():
    """
    Test SystemGenerator.
    """
    from functools import partial
    # Vacuum tests.
    ffxmls = ['amber99sbildn.xml']
    forcefield_kwargs = {'nonbondedMethod': NoCutoff, 'implicitSolvent': None, 'constraints': None}
    for testsystem_name in ['AlanineDipeptideVacuum']:
        f = partial(check_system_generator, ffxmls, forcefield_kwargs, testsystem_name)
        f.description = 'Testing SystemGenerator on %s' % testsystem_name
        yield f

    # Implicit solvent tests.
    ffxmls = ['amber99sbildn.xml', 'amber99_obc.xml']
    forcefield_kwargs = {'nonbondedMethod': NoCutoff, 'implicitSolvent': OBC2, 'constraints': None}
    for testsystem_name in ['AlanineDipeptideImplicit']:
        f = partial(check_system_generator, ffxmls, forcefield_kwargs, testsystem_name)
        f.description = 'Testing SystemGenerator on %s' % testsystem_name
        yield f

    # Small molecule tests.
    gaff_xml_filename = utils.get_data_filename("parameters/gaff.xml")
    ffxmls = [gaff_xml_filename]
    forcefield_kwargs = {'nonbondedMethod': NoCutoff, 'implicitSolvent': None, 'constraints': None}
    for name in IUPAC_molecule_names:
        f = partial(check_system_generator, ffxmls, forcefield_kwargs, name, use_gaff=True)
        f.description = 'Testing SystemGenerator on %s' % name
        yield f

def imatinib_timing():
    print("Loading imatinib...")
    # Load the PDB file.
    from simtk.openmm.app import PDBFile
    pdb_filename = utils.get_data_filename("chemicals/imatinib/imatinib.pdb")
    pdb = PDBFile(pdb_filename)
    # Create a ForceField object.
    gaff_xml_filename = utils.get_data_filename("parameters/gaff.xml")
    forcefield = ForceField(gaff_xml_filename)
    # Add the residue template generator.
    from openmoltools.forcefield_generators import gaffTemplateGenerator
    forcefield.registerTemplateGenerator(gaffTemplateGenerator)
    # Parameterize system.
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 5.0 / unit.picoseconds, 1.0 * unit.femtoseconds)
    # Create Context
    context = openmm.Context(system, integrator)
    context.setPositions(pdb.positions)
    integrator.step(100)

    import time
    nsteps = 10000000
    initial_time = time.time()
    integrator.step(nsteps)
    state = context.getState().getPeriodicBoxVectors() # force dynamics
    final_time = time.time()
    elapsed_time = final_time / initial_time
    time_per_step = elapsed_time / float(nsteps)
    print('time per force evaluation is %.3f us' % (time_per_step*1e6))

class Timer(object):
    def __enter__(self):
        import time
        self.start = time.time()
        return self

    def __exit__(self, *args):
        import time
        self.end = time.time()
        self.interval = self.end - self.start

@pytest.mark.skipif(not HAVE_OE, reason="Cannot test without OpenEye tools")
class TestOEGAFFTemplateGenerator(unittest.TestCase):
    def setUp(self):
        from openeye import oechem
        nmolecules = 5 # number of molecules to read
        ifs = oechem.oemolistream(utils.get_data_filename("chemicals/minidrugbank/MiniDrugBank_tripos.mol2"))
        self.oemols = list() # list of molecules to use
        for index in range(nmolecules):
            oemol = oechem.OEMol()
            oechem.OEReadMolecule(ifs, oemol)
            self.oemols.append(oemol)
        ifs.close()

        # Suppress DEBUG logging from various packages
        import logging
        for name in ['parmed', 'matplotlib']:
            logging.getLogger(name).setLevel(logging.WARNING)

    def test_create(self):
        """Test creation of an OEGAFFTemplateGenerator"""
        from openmoltools.forcefield_generators import OEGAFFTemplateGenerator
        # Create an empty generator
        generator = OEGAFFTemplateGenerator()
        # Create a generator that knows about a few molecules
        generator = OEGAFFTemplateGenerator(oemols=self.oemols)
        # Create a generator that also has a database cache
        with tempfile.TemporaryDirectory() as tmpdirname:
            cache = os.path.join(tmpdirname, 'db.json')
            # Create a new database file
            generator = OEGAFFTemplateGenerator(oemols=self.oemols, cache=cache)
            del generator
            # Reopen it (with cache still empty)
            generator = OEGAFFTemplateGenerator(oemols=self.oemols, cache=cache)
            del generator

    def test_parameterize(self):
        """Test parameterizing molecules with OEGAFFTemplateGenerator"""
        from openmoltools.forcefield_generators import OEGAFFTemplateGenerator
        for gaff_version in ['gaff', 'gaff2']:
            # Create a generator that knows about a few molecules
            generator = OEGAFFTemplateGenerator(oemols=self.oemols)
            # Create a ForceField
            from simtk.openmm.app import ForceField
            gaff_xml_filename = utils.get_data_filename("parameters/{}.xml".format(gaff_version))
            forcefield = ForceField(gaff_xml_filename)
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize some molecules
            from simtk.openmm.app import NoCutoff
            from openmoltools.forcefield_generators import generateTopologyFromOEMol
            for oemol in self.oemols:
                topology = generateTopologyFromOEMol(oemol)
                with Timer() as t1:
                    system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)
                assert system.getNumParticles() == sum([1 for atom in oemol.GetAtoms()])
                # Molecule should now be cached
                with Timer() as t2:
                    system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)
                assert system.getNumParticles() == sum([1 for atom in oemol.GetAtoms()])
                assert (t2.interval < t1.interval)

    def test_add_oemols(self):
        """Test that OEMols can be added to OEGAFFTemplateGenerator later"""
        from openmoltools.forcefield_generators import OEGAFFTemplateGenerator
        gaff_version = 'gaff'
        # Create a generator that does not know about any molecules
        generator = OEGAFFTemplateGenerator()
        # Create a ForceField
        from simtk.openmm.app import ForceField
        gaff_xml_filename = utils.get_data_filename("parameters/{}.xml".format(gaff_version))
        forcefield = ForceField(gaff_xml_filename)
        # Register the template generator
        forcefield.registerTemplateGenerator(generator.generator)

        # Check that parameterizing a molecule fails
        oemol = self.oemols[0]
        from simtk.openmm.app import NoCutoff
        from openmoltools.forcefield_generators import generateTopologyFromOEMol
        try:
            # This should fail with an exception
            system = forcefield.createSystem(generateTopologyFromOEMol(oemol), nonbondedMethod=NoCutoff)
        except ValueError as e:
            # Exception 'No template found...' is expected
            assert str(e).startswith('No template found')

        # Now add the molecule to the generator and ensure parameterization passes
        generator.add_oemols(oemol)
        system = forcefield.createSystem(generateTopologyFromOEMol(oemol), nonbondedMethod=NoCutoff)
        assert system.getNumParticles() == sum([1 for atom in oemol.GetAtoms()])

        # Add multiple molecules, including repeats
        generator.add_oemols(self.oemols)

        # Ensure all molecules can be parameterized
        for oemol in self.oemols:
            system = forcefield.createSystem(generateTopologyFromOEMol(oemol), nonbondedMethod=NoCutoff)
            assert system.getNumParticles() == sum([1 for atom in oemol.GetAtoms()])

    def test_cache(self):
        """Test cache capability of OEGAFFTemplateGenerator"""
        from openmoltools.forcefield_generators import OEGAFFTemplateGenerator, generateTopologyFromOEMol
        from simtk.openmm.app import ForceField, NoCutoff
        gaff_version = 'gaff'
        gaff_xml_filename = utils.get_data_filename("parameters/{}.xml".format(gaff_version))
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Create a generator that also has a database cache
            cache = os.path.join(tmpdirname, 'db.json')
            generator = OEGAFFTemplateGenerator(oemols=self.oemols, cache=cache)
            # Create a ForceField
            forcefield = ForceField(gaff_xml_filename)
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize the molecules
            for oemol in self.oemols:
                forcefield.createSystem(generateTopologyFromOEMol(oemol), nonbondedMethod=NoCutoff)

            # Check database contents
            def check_database(generator):
                db_entries = generator._db.all()
                nentries = len(db_entries)
                nmolecules = len(self.oemols)
                assert (nmolecules == nentries), \
                    "Expected {} entries but database only has {}\n db contents: {}".format(nmolecules, nentries, db_entries)

            check_database(generator)
            # Clean up, forcing closure of database
            del forcefield, generator

            # Create a generator that also uses the database cache but has no molecules
            print('Creating new generator with just cache...')
            generator = OEGAFFTemplateGenerator(cache=cache)
            # Check database contents
            check_database(generator)
            # Create a ForceField
            forcefield = ForceField(gaff_xml_filename)
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize the molecules; this should succeed
            for oemol in self.oemols:
                from openeye import oechem
                forcefield.createSystem(generateTopologyFromOEMol(oemol), nonbondedMethod=NoCutoff)

@pytest.mark.skipif(not HAVE_OE, reason="Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_generate_ffxml_from_molecules():
    """
    Test generation of single ffxml file from a list of molecules
    """
    # Create a test set of molecules.
    molecules = [ createOEMolFromIUPAC(name) for name in IUPAC_molecule_names ]
    # Create an ffxml file.
    from openmoltools.forcefield_generators import generateForceFieldFromMolecules
    ffxml = generateForceFieldFromMolecules(molecules)
    # Create a ForceField.
    gaff_xml_filename = utils.get_data_filename("parameters/gaff.xml")
    forcefield = ForceField(gaff_xml_filename)
    try:
        forcefield.loadFile(StringIO(ffxml))
    except Exception as e:
        msg  = str(e)
        msg += "ffxml contents:\n"
        for (index, line) in enumerate(ffxml.split('\n')):
            msg += 'line %8d : %s\n' % (index, line)
        raise Exception(msg)

    # Parameterize the molecules.
    from openmoltools.forcefield_generators import generateTopologyFromOEMol
    for molecule in molecules:
        # Create topology from molecule.
        topology = generateTopologyFromOEMol(molecule)
        # Create system with forcefield.
        system = forcefield.createSystem(topology)
        # Check potential is finite.
        positions = extractPositionsFromOEMOL(molecule)
        check_potential_is_finite(system, positions)


if __name__ == '__main__':
    imatinib_timing()

    #test_PerceiveBondOrdersExplicitHydrogens(write_pdf=True)
    unittest.main()
