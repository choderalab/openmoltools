import gafftools

pybel_mol2, pandas_mol2 = gafftools.load_antechamber_mol2("./sustiva.mol2")
pybel_mol2.atoms[6].type
pybel_mol2.write("pdb", "out.pdb")

xml = gafftools.generate_gaff_xml(pandas_mol2, pybel_mol2)
f = open("./out.xml", 'w')
f.writelines(xml.readlines())
