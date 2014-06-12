from unittest import skipIf
import tempfile
import os
from gaff2xml import utils

def test_drugs():
    path = tempfile.mkdtemp()
    database_filename = utils.get_data_filename("chemicals/drugs/Zdd.mol2.gz")
    cmd = "gunzip -c %s > %s/Zdd.mol2" % (database_filename, path)
    os.system(cmd)
    cmd = """awk '/MOLECULE/{close(x);x="%s/molecule_"i++".mol2"}{print > x}' %s/Zdd.mol2""" % (path, path)
    os.system(cmd)
    
    n_molecules = 3404
    if os.environ.get("TRAVIS", None) == 'true':
        n_molecules = 25  # If running on travis, only test the first 25 molecules due to speed.

    for k in range(n_molecules):
        molecule_name = "molecule_%d" % k
        mol2_filename = "%s/%s.mol2" % (path, molecule_name)
        cmd = """sed -i "s/<0>/LIG/" %s""" % mol2_filename
        os.system(cmd)  # Have to remove the <0> because it leads to invalid XML in the forcefield files.
        with utils.enter_temp_directory():        
             yield lambda : utils.test_molecule("LIG", mol2_filename)
