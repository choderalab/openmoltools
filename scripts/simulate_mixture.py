import numpy as np
import mdtraj as md
import gaff2xml
import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm

m0 = gaff2xml.openeye.smiles_to_oemol("CCCCC")
charged0 = gaff2xml.openeye.get_charges(m0, keep_confs=1)

ligands = [charged0]

with gaff2xml.utils.enter_temp_directory():  # Supress antechamber output files
    trajectories, ffxml = gaff2xml.openeye.oemols_to_ffxml(ligands)

open("./ff.xml", 'w').write(ffxml.read())
traj = gaff2xml.packmol.pack_box(trajectories, [400])
traj.save("./box.pdb")

temperature = 300 * u.kelvin
friction = 1.0 / u.picosecond
timestep = 2.0 * u.femtosecond

positions = traj.openmm_positions(0)
topology = traj.top.to_openmm(traj)

ffxml.seek(0)
forcefield = app.ForceField(ffxml)

system = forcefield.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=app.HAngles)
system.addForce(mm.MonteCarloBarostat(1 * u.atmospheres, temperature, 25))

integrator = mm.LangevinIntegrator(temperature, friction, timestep / 5.)

simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

simulation.step(10000)
integrator.setStepSize(timestep)

simulation.reporters.append(app.DCDReporter('out.dcd', 500))
simulation.reporters.append(app.StateDataReporter("out.csv", 500, step=True, temperature=True, density=True, energy=True, totalSteps=1000, separator=","))

print("running")
simulation.step(10000000)
