
"""
antechamber -i sustiva.pdb -fi pdb -o sustiva.mol2 -fo mol2 -c bcc -s 2
parmchk -i sustiva.mol2 -f mol2 -o sustiva.frcmod

tleap -f leaprc.ff99SB
source leaprc.gaff
SUS = loadmol2 sustiva.mol2 
check SUS
loadamberparams sustiva.frcmod

"""

import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
import gafftools_new

temperature = 300 * u.kelvin
friction = 0.3 / u.picosecond
timestep = 0.1 * u.femtosecond

prmtop = app.AmberPrmtopFile("./sustiva.prmtop")
inpcrt = app.AmberInpcrdFile("./sustiva.inpcrd")

system_prm = prmtop.createSystem(nonbondedMethod=app.NoCutoff, nonbondedCutoff=1.0*u.nanometers, constraints=None)

mol2 = gafftools_new.Mol2Parser("./sample_files/sustiva.mol2")
top, xyz = mol2.to_openmm()

forcefield = app.ForceField("out.xml")

system_xml = forcefield.createSystem(top, nonbondedMethod=app.NoCutoff, nonbondedCutoff=1.0*u.nanometers, constraints=None)

forces_xml = system_xml.getForces()
forces_xml.sort(key=lambda x: type(x))

forces_prm = system_prm.getForces()
forces_prm.sort(key=lambda x: type(x))

def check_bonds(force0, force1):
    eq(force0.getNumBonds(), force1.getNumBonds())
