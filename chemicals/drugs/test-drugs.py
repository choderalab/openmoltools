#!/usr/bin/env python

import os
import os.path
import tempfile
import sys

from openeye.oechem import *

import simtk.openmm
from simtk.openmm import app
import simtk.unit as units

from gaff2xml.amber_parser import AmberParser

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format="LOG: %(message)s")

def run_antechamber(molecule, charge_method=None):
    """
    Run AmberTools antechamber to create GAFF mol2 and frcmod files.

    Parameters
    ----------
    molecule : openeye.oechem.OEGraphMol 
        The molecule to be parameterized.
    charge_method : str, optional
        If not None, the charge method string will be passed to Antechamber.

    Returns
    -------
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    
    """
    # Get molecule name.
    molecule_name = molecule.GetTitle()
    logger.debug(molecule_name)

    # Write molecule as Tripos mol2.
    tripos_mol2_filename = molecule_name + '.tripos.mol2'
    ofs = oemolostream(tripos_mol2_filename)
    ofs.SetFormat(OEFormat_MOL2H)
    OEWriteMolecule(ofs, molecule)
    ofs.close()

    # Replace <0> substructure names with valid text.
    infile = open(tripos_mol2_filename, 'r')
    lines = infile.readlines()
    infile.close()
    newlines = [ line.replace('<0>', 'MOL') for line in lines ]
    outfile = open(tripos_mol2_filename, 'w')
    outfile.writelines(newlines)
    outfile.close()

    # Run Antechamber to generate parameters.
    gaff_mol2_filename = molecule_name + '.gaff.mol2'
    frcmod_filename = molecule_name + '.frcmod'
    cmd = "antechamber -i %s -fi mol2 -o %s -fo mol2 -s 2" % (tripos_mol2_filename, gaff_mol2_filename)
    if charge_method:
        cmd += ' -c %s' % charge_method
    logger.debug(cmd)
    output = os.system(cmd)
    logger.debug(output)
    cmd = "parmchk -i %s -f mol2 -o %s" % (gaff_mol2_filename, frcmod_filename)
    logger.debug(cmd)
    output = os.system(cmd)
    logger.debug(output)

    return (gaff_mol2_filename, frcmod_filename)

def run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename):
    prmtop_filename = molecule_name + '.prmtop' 
    inpcrd_filename = molecule_name + '.inpcrd'

    # Run tleap to generate prmtop/inpcrd.
    tleap_input = """
source leaprc.ff99SB
source leaprc.gaff
LIG = loadmol2 %(gaff_mol2_filename)s
check LIG
loadamberparams %(frcmod_filename)s
saveamberparm LIG %(prmtop_filename)s %(inpcrd_filename)s
quit

""" % vars()
    leap_input_filename = 'leap.in'
    outfile = open(leap_input_filename, 'w')
    outfile.writelines(tleap_input)
    outfile.close()
    cmd = "tleap -f %s " % leap_input_filename
    output = os.system(cmd)
    logger.debug(output)

    return (prmtop_filename, inpcrd_filename)

def create_ffxml_simulation(molecule, gaff_mol2_filename, frcmod_filename):
    molecule_name = molecule.GetTitle()

    # Generate ffxml file.
    amberhome_path = os.environ['AMBERHOME']
    amber_parser = AmberParser()
    gaff_dat_filename = os.path.join(amberhome_path, 'dat', 'leap', 'parm', 'gaff.dat')
    logger.debug(gaff_dat_filename)
    amber_parser.parse_filenames([gaff_dat_filename, gaff_mol2_filename, frcmod_filename])
    ffxml_stream = amber_parser.generate_xml()
    ffxml_filename = molecule_name + '.ffxml'
    outfile = open(ffxml_filename, 'w')
    outfile.write(ffxml_stream.read())
    outfile.close()

    # Read mol2 file.
    from gaff2xml.gafftools import Mol2Parser
    mol2 = Mol2Parser(gaff_mol2_filename)
    [topology, positions] = mol2.to_openmm()

    # Create System object.
    forcefield = app.ForceField(ffxml_filename)
    system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=None)
    
    # Create integrator.
    timestep = 1.0 * units.femtoseconds
    integrator = simtk.openmm.VerletIntegrator(timestep)
    
    # Create simulation.
    platform = simtk.openmm.Platform.getPlatformByName("Reference")
    simulation = app.Simulation(topology, system, integrator, platform=platform)
    simulation.context.setPositions(positions)

    return simulation

def create_leap_simulation(molecule, gaff_mol2_filename, frcmod_filename):
    molecule_name = molecule.GetTitle()

    # Parameterize system with LEaP.
    [prmtop_filename, inpcrd_filename] = run_tleap(molecule_name, gaff_mol2_filename, frcmod_filename)

    # Create System object.
    prmtop = app.AmberPrmtopFile(prmtop_filename)
    topology = prmtop.topology
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=None)

    # Read positions.
    inpcrd = app.AmberInpcrdFile(inpcrd_filename)
    positions = inpcrd.getPositions()

    # Create integrator.
    timestep = 1.0 * units.femtoseconds
    integrator = simtk.openmm.VerletIntegrator(timestep)
    
    platform = simtk.openmm.Platform.getPlatformByName("Reference")
    simulation = app.Simulation(topology, system, integrator, platform=platform)
    simulation.context.setPositions(positions)

    return simulation
    pass

def test_molecule(molecule, charge_method=None):
    
    # Create temporary directory.
    tmp_path = 'tmp' # DEBUG: Create an actual temporary directory.
    if not os.path.exists(tmp_path): os.makedirs(tmp_path)
    logger.debug('temporary directory created in %s' % tmp_path)
    cwd = os.getcwd()
    os.chdir(tmp_path)
    
    # Generate GAFF parameters.
    [gaff_mol2_filename, frcmod_filename] = run_antechamber(molecule)

    # Create simulations.
    simulation_ffxml = create_ffxml_simulation(molecule, gaff_mol2_filename, frcmod_filename)
    simulation_leap  = create_leap_simulation(molecule, gaff_mol2_filename, frcmod_filename)
    
    # Compare simulations.
    from gaff2xml.system_checker import SystemChecker
    syscheck = SystemChecker(simulation_ffxml, simulation_leap)
    syscheck.check_force_parameters()

    # TODO: Clean up temporary directory.

    # Restore current working directory.
    os.chdir(cwd)

    return

if __name__ == '__main__':

    # Test downloaded drug database.
    database_filename = 'Zdd.mol2.gz' # mol2 database source
    ifs = oemolistream(database_filename)
    for molecule in ifs.GetOEGraphMols():
        # Test molecule.
        test_molecule(molecule)

