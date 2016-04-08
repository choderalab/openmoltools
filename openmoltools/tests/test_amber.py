import numpy as np
import mdtraj as md
from unittest import skipIf
import logging
from openmoltools import utils, amber,packmol
from distutils.spawn import find_executable
import shutil

logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")

@skipIf(packmol.PACKMOL_PATH is None, "Skipping testing of packmol conversion because packmol not found.")
def test_amber_box():
    etoh_filename = utils.get_data_filename("chemicals/etoh/etoh.mol2")
    trj_list = [md.load(etoh_filename)]
    
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        
        box_filename = "./box.pdb"
        box_trj = packmol.pack_box(trj_list, [50])
        box_trj.save(box_filename)
    
        gaff_mol2_filename1, frcmod_filename1 = amber.run_antechamber("etoh", etoh_filename, charge_method=None)
        
        mol2_filenames = [gaff_mol2_filename1]
        frcmod_filenames =  [frcmod_filename1]
        
        prmtop_filename = "./out.prmtop"
        inpcrd_filename = "./out.inpcrd"

        tleap_cmd = amber.build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename)
        print(tleap_cmd)


@skipIf(packmol.PACKMOL_PATH is None, "Skipping testing of packmol conversion because packmol not found.")
def test_amber_binary_mixture():
    sustiva_filename = utils.get_data_filename("chemicals/etoh/etoh.mol2")
    etoh_filename = utils.get_data_filename("chemicals/etoh/etoh_renamed.mol2")

    trj0, trj1 = md.load(sustiva_filename), md.load(etoh_filename)
    
    # Hack to assign unique residue names that are consistent with contents of mol2 files
    trj0.top.residue(0).name = "LIG"
    trj1.top.residue(0).name = "LI2"
    
    trj_list = [trj0, trj1]
    
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        
        box_filename = "./box.pdb"
        box_trj = packmol.pack_box(trj_list, [25, 50])
        box_trj.save(box_filename)
        
        gaff_mol2_filename0, frcmod_filename0 = amber.run_antechamber("sustiva", sustiva_filename, charge_method=None)
        gaff_mol2_filename1, frcmod_filename1 = amber.run_antechamber("etoh", etoh_filename, charge_method=None)
        
        mol2_filenames = [gaff_mol2_filename0, gaff_mol2_filename1]
        frcmod_filenames = [frcmod_filename0, frcmod_filename1]
        
        
        prmtop_filename = "./out.prmtop"
        inpcrd_filename = "./out.inpcrd"

        tleap_cmd = amber.build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename)
        print(tleap_cmd)

@skipIf(packmol.PACKMOL_PATH is None, "Skipping testing of packmol conversion because packmol not found.")
def test_amber_water_mixture():
    water_filename = utils.get_data_filename("chemicals/water/water.mol2")
    etoh_filename = utils.get_data_filename("chemicals/etoh/etoh.mol2")
    sustiva_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")

    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        shutil.copy( water_filename, 'c1.mol2' )
        shutil.copy( etoh_filename, 'c2.mol2' )
        shutil.copy( sustiva_filename, 'c3.mol2')
        water_filename = 'c1.mol2'
        etoh_filename = 'c2.mol2'
        sustiva_filename = 'c3.mol2'
        #Randomize residue names to avoid clashes
        utils.randomize_mol2_residue_names( [ water_filename, etoh_filename, sustiva_filename] )

        trj0, trj1, trj2 = md.load(water_filename), md.load(etoh_filename), md.load(sustiva_filename)

        trj_list = [trj0, trj1, trj2]
        
        box_filename = "./box.pdb"
        box_trj = packmol.pack_box(trj_list, [300, 25, 3])
        box_trj.save(box_filename)
        
        gaff_mol2_filename0, frcmod_filename0 = amber.run_antechamber("water", water_filename, charge_method=None)
        gaff_mol2_filename1, frcmod_filename1 = amber.run_antechamber("etoh", etoh_filename, charge_method=None)
        gaff_mol2_filename2, frcmod_filename2 = amber.run_antechamber("sustiva", sustiva_filename, charge_method=None)
        
        mol2_filenames = [gaff_mol2_filename0, gaff_mol2_filename1, gaff_mol2_filename2]
        frcmod_filenames = [frcmod_filename0, frcmod_filename1, frcmod_filename2]
        
        
        prmtop_filename = "./out.prmtop"
        inpcrd_filename = "./out.inpcrd"

        shutil.copy(box_filename, 'renamed.pdb')
        tleap_cmd = amber.build_mixture_prmtop(mol2_filenames, frcmod_filenames, 'renamed.pdb', prmtop_filename, inpcrd_filename)
        print(tleap_cmd)

        #Also do here for case of GAFF water
        tleap_cmd = amber.build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename, water_model = None)
        print(tleap_cmd)

        #Also do here for case of SPC
        tleap_cmd = amber.build_mixture_prmtop(mol2_filenames, frcmod_filenames, 'renamed.pdb', prmtop_filename, inpcrd_filename, water_model = 'SPC')
        print(tleap_cmd)


def test_run_antechamber():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        gaff_mol2_filename, frcmod_filename = amber.run_antechamber(molecule_name, input_filename, charge_method=None)

def test_run_antechamber_resname():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        gaff_mol2_filename, frcmod_filename = amber.run_antechamber(molecule_name, input_filename, charge_method=None, resname=True)
        with open(gaff_mol2_filename, 'r') as fin:
            fin.readline()
            assert fin.readline().strip() == molecule_name
                    
def test_run_tleap():
    molecule_name = "sustiva"
    input_filename = utils.get_data_filename("chemicals/sustiva/sustiva.mol2")
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        gaff_mol2_filename, frcmod_filename = amber.run_antechamber(molecule_name, input_filename, charge_method=None)
        prmtop, inpcrd = amber.run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename)

def test_run_antechamber_charges():
    molecule_name = "acetate"
    input_filename = utils.get_data_filename("chemicals/acetate/acetate.mol2")
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        gaff_mol2_filename, frcmod_filename = amber.run_antechamber(molecule_name, input_filename, charge_method=None, net_charge=-1)

def test_generate_amino_acids():
    """
    Try to generate capped versions of each amino acid
    """
    aminos = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    for aa in aminos:
        aa_list = ['ACE', aa, 'NME']
        topology, position = amber.build_peptide_tleap(aa_list)
        residues = list(topology.residues())
        for i, aa in enumerate(aa_list):
            assert aa == residues[i].name
