from unittest import skipIf
from openmoltools import utils, amber, packmol, gromacs
from distutils.spawn import find_executable
import os

def test_gromacs_merge():
    etoh_filename = utils.get_data_filename("chemicals/etoh/etoh.mol2")
    benzene_filename = utils.get_data_filename("chemicals/benzene/benzene.mol2")

    with utils.enter_temp_directory(): #Prevents creating lots of tleap/antechamber files everywhere
        #Generate frcmod files, mol2 files
        gaff_mol2_filename1, frcmod_filename1 = amber.run_antechamber( "etoh", etoh_filename, charge_method = None)
        gaff_mol2_filename2, frcmod_filename2 = amber.run_antechamber( "benzene", benzene_filename, charge_method = None)

        #Set file names
        prmtop_filename1 = "./out1.prmtop"
        prmtop_filename2 = "./out2.prmtop"
        crd_filename1 = "./out1.inpcrd"
        crd_filename2 = "./out2.inpcrd"
        top_filename1 = "./out1.top"
        top_filename2 = "./out2.top"
        gro_filename1 = "./out1.gro"
        gro_filename2 = "./out2.gro"

        #Generate AMBER files
        amber.run_tleap( 'etoh', gaff_mol2_filename1, frcmod_filename1, prmtop_filename1, crd_filename1 )
        amber.run_tleap( 'benzene', gaff_mol2_filename2, frcmod_filename2, prmtop_filename2, crd_filename2 )

        #Convert to GROMACS
        amber.convert_via_acpype( "etoh", prmtop_filename1, crd_filename1, out_top = top_filename1, out_gro = gro_filename1 ) 
        amber.convert_via_acpype( "benzene", prmtop_filename2, crd_filename2, out_top = top_filename2, out_gro = gro_filename2 )

        #Merge topologies
        gromacs.merge_topologies( [ top_filename1, top_filename2], './combined.top', 'combined', molecule_numbers = [1, 5] )

        #Test editing of molecule numbers in topology file
        gromacs.change_molecules_section( './combined.top', './edited.top', ['etoh', 'benzene'], [10, 20] )

@skipIf(gromacs.GROMACS_PATH is None, "Skipping testing of GROMACS solvation because GROMACS not found.")
def test_gromacs_solvate():
    etoh_filename = utils.get_data_filename("chemicals/etoh/etoh.mol2")
    with utils.enter_temp_directory(): #Prevents creating lots of tleap/antechamber files everywhere
        #Generate frcmod files, mol2 files
        gaff_mol2_filename, frcmod_filename = amber.run_antechamber( "etoh", etoh_filename, charge_method = None)

        #Amber setup
        amber.run_tleap( 'etoh', gaff_mol2_filename, frcmod_filename, 'etoh.prmtop', 'etoh.crd' )
        #GROMACS conversion
        utils.convert_via_acpype( 'etoh', 'etoh.prmtop', 'etoh.crd', 'etoh.top', 'etoh.gro' )
        #Solvate
        gromacs.do_solvate( 'etoh.top', 'etoh.gro', 'etoh_solvated.top', 'etoh_solvated.gro', 1.2, 'dodecahedron', 'spc216', 'tip3p.itp' )
 
