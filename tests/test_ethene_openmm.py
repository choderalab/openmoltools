import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
import gafftools
import system_checker

temperature = 300 * u.kelvin
friction = 0.3 / u.picosecond
timestep = 0.1 * u.femtosecond

prmtop = app.AmberPrmtopFile("./ethene/ethene.prmtop")
inpcrt = app.AmberInpcrdFile("./ethene/ethene.inpcrd")

system_prm = prmtop.createSystem(nonbondedMethod=app.NoCutoff, nonbondedCutoff=1.0*u.nanometers, constraints=None)

mol2 = gafftools.Mol2Parser("./ethene/ethene.mol2")
top, xyz = mol2.to_openmm()

forcefield = app.ForceField("ethene.xml")

system_xml = forcefield.createSystem(top, nonbondedMethod=app.NoCutoff, nonbondedCutoff=1.0*u.nanometers, constraints=None)

checker = system_checker.SystemChecker(system_prm, system_xml)
checker.check_forces()


integrator_xml = mm.LangevinIntegrator(temperature, friction, timestep)
simulation_xml = app.Simulation(top, system_xml, integrator_xml)
simulation_xml.context.setPositions(xyz)
state_xml = simulation_xml.context.getState(getEnergy=True)
state_xml.getPotentialEnergy()

integrator_prm = mm.LangevinIntegrator(temperature, friction, timestep)
simulation_prm = app.Simulation(prmtop.topology, system_prm, integrator_prm)
simulation_prm.context.setPositions(xyz)
state_prm = simulation_prm.context.getState(getEnergy=True)
state_prm.getPotentialEnergy()
