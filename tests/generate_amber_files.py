
"""
antechamber -i sustiva.pdb -fi pdb -o sustiva.mol2 -fo mol2 -c bcc -s 2
parmchk -i sustiva.mol2 -f mol2 -o sustiva.frcmod

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
