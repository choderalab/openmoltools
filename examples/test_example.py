import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
from gaff2xml import gafftools, system_checker

example = "imatinib"
temperature = 300 * u.kelvin
friction = 0.3 / u.picosecond
timestep = 0.1 * u.femtosecond

prmtop = app.AmberPrmtopFile("./chemicals/%s/%s.prmtop" % (example, example))
inpcrt = app.AmberInpcrdFile("./chemicals/%s/%s.inpcrd" % (example, example))

system_prm = prmtop.createSystem(nonbondedMethod=app.NoCutoff, nonbondedCutoff=1.0*u.nanometers, constraints=None)

mol2 = gafftools.Mol2Parser("./chemicals/%s/%s.mol2" % (example, example))
top, xyz = mol2.to_openmm()

forcefield = app.ForceField("./chemicals/%s/%s.xml" % (example, example))

system_xml = forcefield.createSystem(top, nonbondedMethod=app.NoCutoff, nonbondedCutoff=1.0*u.nanometers, constraints=None)

integrator_xml = mm.LangevinIntegrator(temperature, friction, timestep)
simulation_xml = app.Simulation(top, system_xml, integrator_xml)
simulation_xml.context.setPositions(xyz)

integrator_prm = mm.LangevinIntegrator(temperature, friction, timestep)
simulation_prm = app.Simulation(prmtop.topology, system_prm, integrator_prm)
simulation_prm.context.setPositions(xyz)

checker = system_checker.SystemChecker(simulation_xml, simulation_prm)
checker.check_force_parameters()
energy0, energy1 = checker.check_energies()

abs((energy0 - energy1) / u.kilojoules_per_mole)
