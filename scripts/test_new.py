import gafftools_new

mol2 = gafftools_new.Mol2Parser("./sample_files/sustiva.mol2")
traj = mol2.to_mdtraj()
top, xyz = mol2.to_openmm()

xml = gafftools_new.generate_gaff_xml(mol2.atoms, mol2.bonds)
f = open("./out.xml", 'w')
f.writelines(xml.readlines())
