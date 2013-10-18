import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
import gafftools_new

temperature = 300 * u.kelvin
friction = 0.3 / u.picosecond
timestep = 0.1 * u.femtosecond

mol2 = gafftools_new.Mol2Parser("./sample_files/sustiva.mol2")
top, xyz = mol2.to_openmm()

forcefield = app.ForceField("out.xml")

system = forcefield.createSystem(top, nonbondedMethod=app.NoCutoff, nonbondedCutoff=1.0*u.nanometers, constraints=None)

integrator = mm.LangevinIntegrator(temperature, friction, timestep)

simulation = app.Simulation(top, system, integrator)
simulation.context.setPositions(xyz)

simulation.reporters.append(app.PDBReporter("simulation.pdb", 50))
print("running")
simulation.step(5000)
