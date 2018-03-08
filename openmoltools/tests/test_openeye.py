from nose.plugins.attrib import attr
import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
import numpy as np
import re
from mdtraj.testing import eq
from unittest import skipIf
from openmoltools import utils, packmol
import os
import openmoltools.openeye
import pandas as pd
import mdtraj as md
from numpy.testing import assert_raises

smiles_fails_with_strictStereo = "CN1CCN(CC1)CCCOc2cc3c(cc2OC)C(=[NH+]c4cc(c(cc4Cl)Cl)OC)C(=C=[N-])C=[NH+]3"

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

try:
    import parmed
    HAVE_PARMED = True
except ImportError:
    HAVE_PARMED = False


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_butanol_keepconfs():
    m0 = openmoltools.openeye.iupac_to_oemol("butanol")
    m1 = openmoltools.openeye.get_charges(m0, keep_confs=1)
    eq(m0.NumAtoms(), m1.NumAtoms())
    assert m1.NumConfs() == 1, "This OEMol was created to have a single conformation."
    assert m1.NumAtoms() == 15, "Butanol should have 15 atoms"

@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.\n" + openeye_exception_message)
def test_butanol_unnormalized():
    m0 = openmoltools.openeye.iupac_to_oemol("butanol")
    m0.SetTitle("MyCustomTitle")
    m1 = openmoltools.openeye.get_charges(m0, normalize=False, keep_confs=1)
    eq(m0.NumAtoms(), m1.NumAtoms())
    assert m1.NumConfs() == 1, "This OEMol was created to have a single conformation."
    assert m1.NumAtoms() == 15, "Butanol should have 15 atoms"
    assert m0.GetTitle() == m1.GetTitle(), "The title of the molecule should not be changed by normalization."


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_output_mol2():
    molecule = openmoltools.openeye.iupac_to_oemol("cyclopentane")
    openmoltools.openeye.molecule_to_mol2(molecule, tripos_mol2_filename="testing mol2 output.tripos.mol2")


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_output_mol2_standardize():
    molecule = openmoltools.openeye.iupac_to_oemol("cyclopentane")
    list(molecule.GetAtoms())[0].SetName("MyNameIsAtom")
    openmoltools.openeye.molecule_to_mol2(molecule, tripos_mol2_filename="testing mol2 standardize output.tripos.mol2", standardize=True)
    with open("testing mol2 standardize output.tripos.mol2", "r") as outfile:
        text = outfile.read()
    # This should not find the text we added, to make sure the molecule is standardized.
    assert re.search("MyNameIsAtom", text) is None


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_output_mol2_no_standardize():
    molecule = openmoltools.openeye.iupac_to_oemol("cyclopentane")
    list(molecule.GetAtoms())[0].SetName("MyNameIsAtom")
    openmoltools.openeye.molecule_to_mol2(molecule, tripos_mol2_filename="testing mol2 nostandardize output.tripos.mol2", standardize=False)
    with open("testing mol2 nostandardize output.tripos.mol2", "r") as outfile:
        text = outfile.read()
    # This should find the text we added, to make sure the molecule is not standardized.
    assert re.search("MyNameIsAtom", text) is not None


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_output_mol2_multiple_confs():
    molecule = openmoltools.openeye.iupac_to_oemol("butanol")
    multiple_conformers = openmoltools.openeye.generate_conformers(molecule)
    openmoltools.openeye.molecule_to_mol2(multiple_conformers, tripos_mol2_filename="testing mol2 multiple conformers.tripos.mol2", conformer=None)
    with open("testing mol2 multiple conformers.tripos.mol2", "r") as outfile:
        text = outfile.read()
    # This should not find the text we added, to make sure the molecule is standardized.
    assert text.count("@<TRIPOS>MOLECULE") > 1


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_butanol():
    m0 = openmoltools.openeye.iupac_to_oemol("butanol")
    m1 = openmoltools.openeye.get_charges(m0, keep_confs=-1)
    eq(m0.NumAtoms(), m1.NumAtoms())
    assert m1.NumConfs() >= 2, "Butanol should have multiple conformers."
    assert m1.NumAtoms() == 15, "Butanol should have 15 atoms"

    all_data = {}
    for k, molecule in enumerate(m1.GetConfs()):
        names_to_charges, str_repr = openmoltools.openeye.get_names_to_charges(molecule)
        all_data[k] = names_to_charges
        eq(sum(names_to_charges.values()), 0.0, decimal=7)  # Net charge should be zero

    # Build a table of charges indexed by conformer number and atom name
    all_data = pd.DataFrame(all_data)

    # The standard deviation along the conformer axis should be zero if all conformers have same charges
    eq(all_data.std(1).values, np.zeros(m1.NumAtoms()), decimal=7)

    with utils.enter_temp_directory():
        # Try saving to disk as mol2
        openmoltools.openeye.molecule_to_mol2(m1, "out.mol2")
        # Make sure MDTraj can read the output
        t = md.load("out.mol2")
        # Make sure MDTraj can read the charges / topology info
        atoms, bonds = md.formats.mol2.mol2_to_dataframes("out.mol2")

        # Finally, make sure MDTraj and OpenEye report the same charges.
        names_to_charges, str_repr = openmoltools.openeye.get_names_to_charges(m1)
        q = atoms.set_index("name").charge
        q0 = pd.Series(names_to_charges)
        delta = q - q0  # An object containing the charges, with atom names as indices
        eq(delta.values, np.zeros_like(delta.values), decimal=4)


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_benzene():
    m0 = openmoltools.openeye.iupac_to_oemol("benzene")
    m1 = openmoltools.openeye.get_charges(m0)
    eq(m0.NumAtoms(), m1.NumAtoms())
    print(m1.NumConfs())
    assert m1.NumConfs() == 1, "Benezene should have 1 conformer"
    assert m1.NumAtoms() == 12, "Benezene should have 12 atoms"

    names_to_charges, str_repr = openmoltools.openeye.get_names_to_charges(m1)
    eq(sum(names_to_charges.values()), 0.0, decimal=7)  # Net charge should be zero

    with utils.enter_temp_directory():
        # Try saving to disk as mol2
        openmoltools.openeye.molecule_to_mol2(m1, "out.mol2")
        # Make sure MDTraj can read the output
        t = md.load("out.mol2")
        # Make sure MDTraj can read the charges / topology info
        atoms, bonds = md.formats.mol2.mol2_to_dataframes("out.mol2")

        # Finally, make sure MDTraj and OpenEye report the same charges.
        names_to_charges, str_repr = openmoltools.openeye.get_names_to_charges(m1)
        q = atoms.set_index("name").charge
        q0 = pd.Series(names_to_charges)
        delta = q - q0  # An object containing the charges, with atom names as indices
        eq(delta.values, np.zeros_like(delta.values), decimal=4)


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_link_in_utils():
    m0 = openmoltools.openeye.iupac_to_oemol("benzene")
    m1 = openmoltools.openeye.get_charges(m0)
    with utils.enter_temp_directory():
        # This function was moved from utils to openeye, so check that the old link still works.
        utils.molecule_to_mol2(m1, "out.mol2")


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_smiles():
    m0 = openmoltools.openeye.smiles_to_oemol("CCCCO")
    charged0 = openmoltools.openeye.get_charges(m0)

    m1 = openmoltools.openeye.iupac_to_oemol("butanol")
    charged1 = openmoltools.openeye.get_charges(m1)

    eq(charged0.NumAtoms(), charged1.NumAtoms())


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_ffxml():
    with utils.enter_temp_directory():
        m0 = openmoltools.openeye.smiles_to_oemol("CCCCO")
        charged0 = openmoltools.openeye.get_charges(m0)
        m1 = openmoltools.openeye.smiles_to_oemol("ClC(Cl)(Cl)Cl")
        charged1 = openmoltools.openeye.get_charges(m1)

        trajectories, ffxml = openmoltools.openeye.oemols_to_ffxml([charged0, charged1])


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_ffxml_simulation():
    """Test converting toluene and benzene smiles to oemol to ffxml to openmm simulation."""
    with utils.enter_temp_directory():
        m0 = openmoltools.openeye.smiles_to_oemol("Cc1ccccc1")
        charged0 = openmoltools.openeye.get_charges(m0)
        m1 = openmoltools.openeye.smiles_to_oemol("c1ccccc1")
        charged1 = openmoltools.openeye.get_charges(m1)
        ligands = [charged0, charged1]
        n_atoms = [15,12]

        trajectories, ffxml = openmoltools.openeye.oemols_to_ffxml(ligands)
        eq(len(trajectories),len(ligands))

        pdb_filename = utils.get_data_filename("chemicals/proteins/1vii.pdb")

        temperature = 300 * u.kelvin
        friction = 0.3 / u.picosecond
        timestep = 0.01 * u.femtosecond

        protein_traj = md.load(pdb_filename)
        protein_traj.center_coordinates()

        protein_top = protein_traj.top.to_openmm()
        protein_xyz = protein_traj.openmm_positions(0)

        for k, ligand in enumerate(ligands):
            ligand_traj = trajectories[k]
            ligand_traj.center_coordinates()

            eq(ligand_traj.n_atoms, n_atoms[k])
            eq(ligand_traj.n_frames, 1)

            #Move the pre-centered ligand sufficiently far away from the protein to avoid a clash.
            min_atom_pair_distance = ((ligand_traj.xyz[0] ** 2.).sum(1) ** 0.5).max() + ((protein_traj.xyz[0] ** 2.).sum(1) ** 0.5).max() + 0.3
            ligand_traj.xyz += np.array([1.0, 0.0, 0.0]) * min_atom_pair_distance

            ligand_xyz = ligand_traj.openmm_positions(0)
            ligand_top = ligand_traj.top.to_openmm()

            ffxml.seek(0)
            forcefield = app.ForceField("amber10.xml", ffxml, "tip3p.xml")

            model = app.modeller.Modeller(protein_top, protein_xyz)
            model.add(ligand_top, ligand_xyz)
            model.addSolvent(forcefield, padding=0.4 * u.nanometer)

            system = forcefield.createSystem(model.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=app.HAngles)

            integrator = mm.LangevinIntegrator(temperature, friction, timestep)

            simulation = app.Simulation(model.topology, system, integrator)
            simulation.context.setPositions(model.positions)
            print("running")
            simulation.step(1)


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_charge_fail1():
    with assert_raises(RuntimeError):
        with utils.enter_temp_directory():
            openmoltools.openeye.smiles_to_antechamber(smiles_fails_with_strictStereo, "test.mol2",  "test.frcmod", strictStereo=True)


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_charge_fail2():
    with assert_raises(RuntimeError):
        m = openmoltools.openeye.smiles_to_oemol(smiles_fails_with_strictStereo)
        m = openmoltools.openeye.get_charges(m, strictStereo=True, keep_confs=1)


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_charge_success1():
    with utils.enter_temp_directory():
        openmoltools.openeye.smiles_to_antechamber(smiles_fails_with_strictStereo, "test.mol2",  "test.frcmod", strictStereo=False)


@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_charge_success2():
    m = openmoltools.openeye.smiles_to_oemol(smiles_fails_with_strictStereo)
    m = openmoltools.openeye.get_charges(m, strictStereo=False)

@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_oeassigncharges_fail():
    with assert_raises(RuntimeError):
        # Fail test for OEToolkits (2017.2.1) new charging function
        m = openmoltools.openeye.smiles_to_oemol(smiles_fails_with_strictStereo)
        m = openmoltools.openeye.get_charges(m,  strictStereo=False, legacy=False)

@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
def test_oeassigncharges_success():
    # Success test for OEToolkits (2017.2.1) new charging function
    m = openmoltools.openeye.iupac_to_oemol("butanol")
    m = openmoltools.openeye.get_charges(m, legacy=False)

@skipIf(not HAVE_OE, "Cannot test openeye module without OpenEye tools.")
@skipIf(not HAVE_PARMED, "Cannot test without Parmed Chemistry.")
@skipIf(packmol.PACKMOL_PATH is None, "Skipping testing of packmol conversion because packmol not found.")
@attr("parmed")
def test_binary_mixture_rename():
    smiles_string0 = "CCCCCC"
    smiles_string1 = "CCCCCCCCC"

    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        mol2_filename0 = "./A.mol2"
        frcmod_filename0 = "./A.frcmod"

        mol2_filename1 = "./B.mol2"
        frcmod_filename1 = "./B.frcmod"

        gaff_mol2_filenames = [mol2_filename0, mol2_filename1]
        frcmod_filenames = [frcmod_filename0, frcmod_filename1]

        prmtop_filename = "./box.prmtop"
        inpcrd_filename = "./box.inpcrd"

        openmoltools.openeye.smiles_to_antechamber(smiles_string0, mol2_filename0, frcmod_filename0)
        openmoltools.openeye.smiles_to_antechamber(smiles_string1, mol2_filename1, frcmod_filename1)

        openmoltools.utils.randomize_mol2_residue_names(gaff_mol2_filenames)

        box_pdb_filename = "./box.pdb"

        gaff_mol2_filenames = [mol2_filename0, mol2_filename1]
        n_monomers = [10, 20]

        packed_trj = packmol.pack_box([md.load(mol2) for mol2 in gaff_mol2_filenames], n_monomers)
        packed_trj.save(box_pdb_filename)

        tleap_cmd = openmoltools.amber.build_mixture_prmtop(gaff_mol2_filenames, frcmod_filenames, box_pdb_filename, prmtop_filename, inpcrd_filename)

        prmtop = app.AmberPrmtopFile(prmtop_filename)
        inpcrd = app.AmberInpcrdFile(inpcrd_filename)

        system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*u.nanometers, constraints=app.HBonds)
