import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
import gafftools
import mdtraj

ligand_name = "sustiva"

temperature = 300 * u.kelvin
friction = 0.3 / u.picosecond
timestep = 2.0 * u.femtosecond

protein_traj = mdtraj.load("./chemicals/1vii.pdb")
protein_traj.center_coordinates()

protein_top = protein_traj.top.to_openmm()
protein_xyz = protein_traj.openmm_positions(0)

mol2 = gafftools.Mol2Parser("./chemicals/%s/%s.mol2" % (ligand_name, ligand_name))
ligand_top = mol2.to_openmm()[0]

ligand_traj = mol2.to_mdtraj()
ligand_traj.center_coordinates()

#Move the pre-centered ligand sufficiently far away from the protein to avoid a clash.  
min_atom_pair_distance = ((ligand_traj.xyz[0] ** 2.).sum(1) ** 0.5).max() + ((protein_traj.xyz[0] ** 2.).sum(1) ** 0.5).max() + 0.3
ligand_traj.xyz += np.array([1.0, 0.0, 0.0]) * min_atom_pair_distance

ligand_xyz = ligand_traj.openmm_positions(0)

forcefield = app.ForceField("amber10.xml", "./chemicals/%s/%s.xml" % (ligand_name, ligand_name), "tip3p.xml")

model = app.modeller.Modeller(protein_top, protein_xyz)
model.add(ligand_top, ligand_xyz)
model.addSolvent(forcefield, padding=0.4 * u.nanometer)

system = forcefield.createSystem(model.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=app.HAngles)

integrator = mm.LangevinIntegrator(temperature, friction, timestep)

simulation = app.Simulation(model.topology, system, integrator)
simulation.context.setPositions(model.positions)

simulation.minimizeEnergy(tolerance=2.0)

simulation.reporters.append(app.PDBReporter("simulation.pdb", 250))
print("running")
simulation.step(2500)
