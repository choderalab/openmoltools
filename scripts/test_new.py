import simtk.openmm.app as app
import gafftools_new
import mdtraj

mol2 = gafftools_new.Mol2Parser("./sample_files/sustiva.mol2")
traj = mol2.to_mdtraj()

top, xyz = mol2.to_openmm()
