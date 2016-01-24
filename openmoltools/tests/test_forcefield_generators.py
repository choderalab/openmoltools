from nose.plugins.attrib import attr
import unittest
from unittest import skipIf
import tempfile
import os, sys
from openmoltools import utils, forcefield_generators
from simtk.openmm.app import ForceField, NoCutoff
if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO

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

IUPAC_molecule_names = ['ibuprofen', 'aspirin', 'imatinib', 'bosutinib']
def createOEMolFromIUPAC(iupac_name='ibuprofen'):
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
    omega.SetStrictStereo(False)
    omega(mol)

    return mol

@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def disable_PerceiveBondOrdersExplicitHydrogens(write_pdf=False):
    molecule = createOEMolFromIUPAC('ibuprofen')
    from openmoltools.forcefield_generators import PerceiveBondOrdersExplicitHydrogens
    mols = PerceiveBondOrdersExplicitHydrogens(molecule)

    if write_pdf:
        # DEBUG
        from openeye import oedepict
        options = oedepict.OEReportOptions()
        report = oedepict.OEReport(options)
        for mol in mols:
            oedepict.OEPrepareDepiction(mol)
        for mol in mols:
            cell = report.NewCell()
            opts = oedepict.OE2DMolDisplayOptions()
            disp = oedepict.OE2DMolDisplay(mol, opts)
            oedepict.OERenderMolecule(cell, disp)
        oedepict.OEWriteReport('output.pdf', report)

class TestForceFieldGenerators(unittest.TestCase):
    @skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
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
            # TODO: Check bonds.
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

@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_generateResidueTemplate():
    """
    Test GAFF residue template generation from OEMol molecules.
    """
    from openeye import oechem, oeiupac
    # TODO: Test more molecules.
    mol = createOEMolFromIUPAC('ibuprofen')
    # Generate an ffxml residue template.
    from openmoltools.forcefield_generators import generateResidueTemplate
    [template, ffxml] = generateResidueTemplate(mol)
    # Create a ForceField object.
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml', 'gaff.xml')
    # Add the additional parameters and template to the forcefield.
    forcefield.registerResidueTemplate(template)
    forcefield.load_file(StringIO(ffxml))
    # Create a Topology from the molecule.
    from openmoltools.forcefield_generators import generateOpenMMTopology
    topology = generateOpenMMTopology(molecule)
    # Parameterize system.
    system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)

def disable_gaffResidueTemplateGenerator():
    """
    Test the GAFF residue template generator.
    """

    #
    # Test where we generate parameters for only a ligand.
    #

    # Load the PDB file.
    from simtk.openmm.app import PDBFile
    pdb_filename = utils.get_data_filename("chemicals/proteins/T4-lysozyme-L99A-p-xylene-implicit.pdb")
    pdb = PDBFile(pdb_filename)
    # Create a ForceField object.
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml', 'gaff.xml')
    # Add the residue template generator.
    from openmoltools.forcefield_generators import gaffTemplateGenerator
    forcefield.registerTemplateGenerator(gaffTemplateGenerator)
    # Parameterize system.
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    # TODO: Test energies are finite?

#@unittest.skipIf(os.getenv('AMBERHOME') is None, 'Cannot test w/out Amber')
def testWriteXMLParametersGAFF():
    """ Test writing XML parameters loaded from Amber GAFF parameter files """

    # Generate ffxml file contents for parmchk-generated frcmod output.
    leaprc = StringIO("parm = loadamberparams gaff.dat")
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

if __name__ == '__main__':
    #test_PerceiveBondOrdersExplicitHydrogens(write_pdf=True)
    unittest.main()
