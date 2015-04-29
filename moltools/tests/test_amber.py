import numpy as np
import mdtraj as md
from unittest import skipIf
import logging
from gaff2xml import utils, amber,packmol
from distutils.spawn import find_executable


logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")

@skipIf(packmol.PACKMOL_PATH is None, "Skipping testing of packmol conversion because packmol not found.")
def test_amber_box():
    etoh_filename = utils.get_data_filename("chemicals/etoh/etoh.mol2")
    trj_list = [md.load(etoh_filename)]
    
    with utils.enter_temp_directory():  # Prevents creating tons of GAFF files everywhere.
        
        box_filename = "./box.pdb"
        box_trj = packmol.pack_box(trj_list, [50])
        box_trj.save(box_filename)
    
        gaff_mol2_filename1, frcmod_filename1 = utils.run_antechamber("etoh", etoh_filename, charge_method=None)
        
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
        
        gaff_mol2_filename0, frcmod_filename0 = utils.run_antechamber("sustiva", sustiva_filename, charge_method=None)
        gaff_mol2_filename1, frcmod_filename1 = utils.run_antechamber("etoh", etoh_filename, charge_method=None)
        
        mol2_filenames = [gaff_mol2_filename0, gaff_mol2_filename1]
        frcmod_filenames = [frcmod_filename0, frcmod_filename1]
        
        
        prmtop_filename = "./out.prmtop"
        inpcrd_filename = "./out.inpcrd"

        tleap_cmd = amber.build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename)
        print(tleap_cmd)
