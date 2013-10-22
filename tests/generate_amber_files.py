
"""
obabel -i sdf sustiva.sdf -o pdb > sustiva.pdb
antechamber -i sustiva.pdb -fi pdb -o sustiva.mol2 -fo mol2 -c bcc -s 2
parmchk -i sustiva.mol2 -f mol2 -o sustiva.frcmod

python ../../scripts/processAmberForceField.py ~/src/amber12/dat/leap/parm/gaff.dat ./sustiva.mol2 ./sustiva.frcmod > sustiva.xml

tleap -f leaprc.ff99SB
source leaprc.gaff
SUS = loadmol2 sustiva.mol2 
check SUS
loadamberparams sustiva.frcmod
saveoff SUS sus.lib 
saveamberparm SUS sustiva.prmtop sustiva.inpcrd

"""

"""
SUS = loadmol2 etoh.mol2 
check SUS
loadamberparams etoh.frcmod
saveoff SUS etoh.lib 
saveamberparm SUS etoh.prmtop etoh.inpcrd
"""


"""
obabel -i sdf benzene.sdf -o pdb > benzene.pdb
antechamber -i benzene.pdb -fi pdb -o benzene.mol2 -fo mol2 -c bcc -s 2
parmchk -i benzene.mol2 -f mol2 -o benzene.frcmod

python  ../../scripts/processAmberForceField.py ~/src/amber12/dat/leap/parm/gaff.dat ./benzene.mol2 ./benzene.frcmod > benzene.xml


tleap -f leaprc.ff99SB
source leaprc.gaff
SUS = loadmol2 benzene.mol2 
check SUS
loadamberparams benzene.frcmod
saveoff SUS benzene.lib 
saveamberparm SUS benzene.prmtop benzene.inpcrd
"""



"""
obabel -i sdf cyclopropane.sdf -o pdb > cyclopropane.pdb
antechamber -i cyclopropane.pdb -fi pdb -o cyclopropane.mol2 -fo mol2 -c bcc -s 2
parmchk -i cyclopropane.mol2 -f mol2 -o cyclopropane.frcmod

python ../../scripts/processAmberForceField.py ~/src/amber12/dat/leap/parm/gaff.dat ./cyclopropane.mol2 ./cyclopropane.frcmod > cyclopropane.xml


tleap -f leaprc.ff99SB
source leaprc.gaff
SUS = loadmol2 cyclopropane.mol2 
check SUS
loadamberparams cyclopropane.frcmod
saveoff SUS cyclopropane.lib 
saveamberparm SUS cyclopropane.prmtop cyclopropane.inpcrd
"""


"""
obabel -i sdf propene.sdf -o pdb > propene.pdb
antechamber -i propene.pdb -fi pdb -o propene.mol2 -fo mol2 -c bcc -s 2
parmchk -i propene.mol2 -f mol2 -o propene.frcmod

python ../../scripts/processAmberForceField.py ~/src/amber12/dat/leap/parm/gaff.dat ./propene.mol2 ./propene.frcmod > propene.xml


tleap -f leaprc.ff99SB
source leaprc.gaff
SUS = loadmol2 propene.mol2 
check SUS
loadamberparams propene.frcmod
saveoff SUS propene.lib 
saveamberparm SUS propene.prmtop propene.inpcrd
"""
