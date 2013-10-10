from simtk.openmm import app
import simtk.openmm as mm
import simtk.unit as u

temperature = 300 * u.kelvin
friction = 0.3 / u.picosecond
timestep = 1.0 * u.femtosecond

pdb = app.PDBFile('out.pdb')
forcefield = app.ForceField("./xml/gaff.xml", "out.xml")

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, nonbondedCutoff=1.0*u.nanometers, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, friction, timestep)

simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

simulation.reporters.append(app.PDBReporter("simulation.pdb", 50))
print("running")
simulation.step(5000)
