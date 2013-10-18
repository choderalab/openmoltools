import simtk.openmm.app as app
import gafftools
import os

mol2 = gafftools.Mol2Parser("./sustiva.mol2")
traj = mol2.to_mdtraj()
top, xyz = mol2.to_openmm()

xml = gafftools.generate_gaff_xml(mol2.atoms, mol2.bonds)
f = open("./sustiva.xml", 'w')
f.writelines(xml.readlines())

app.PDBFile.writeFile(top, traj.xyz[0], open("test.pdb", "w"))

os.system("python ../src/processAmberForceField.py sustiva.frcmod >> sustiva.frcmod.xml")
