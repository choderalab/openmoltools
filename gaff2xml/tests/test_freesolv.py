from nose.plugins.attrib import attr
from unittest import skipIf
import tempfile
import os
from gaff2xml import utils
import numpy as np
import mdtraj as md
from distutils.spawn import find_executable
import tarfile
import pickle
import os
import numpy as np



@skipIf(find_executable('obabel') is None, 'You need obabel installed to run this test')
def _tester_load_freesolv_gaffmol2_vs_sybylmol2_vs_obabelpdb(charge_method="bcc"):
    with utils.enter_temp_directory():        
        
        tar_filename = utils.get_data_filename("chemicals/freesolv/freesolve_v0.3.tar.bz2")
        tar = tarfile.open(tar_filename, mode="r:bz2")
        tar.extractall()
        tar.close()

        database = pickle.load(open("./v0.3/database.pickle"))
        for key in database:
            for directory in ["mol2files_gaff", "mol2files_sybyl"]:
                gaff_filename = os.path.abspath("./v0.3/%s/%s.mol2" % (directory, key))
                
                cmd = """sed -i "s/<0>/LIG/" %s""" % gaff_filename
                os.system(cmd)  # Have to remove the <0> because it leads to invalid XML in the forcefield files.
                
                t_gaff = md.load(gaff_filename)

                with utils.enter_temp_directory():        
                     yield utils.tag_description(lambda : utils.test_molecule("LIG", gaff_filename, charge_method=charge_method), "Testing freesolv %s %s with charge model %s" % (directory, key, charge_method))

@attr("slow")
def test_load_freesolv_gaffmol2_vs_sybylmol2_vs_obabelpdb():
    _tester_load_freesolv_gaffmol2_vs_sybylmol2_vs_obabelpdb()

# Faster version because it skips AM1-BCC
def test_load_freesolv_gaffmol2_vs_sybylmol2_vs_obabelpdb_nobcc():
    _tester_load_freesolv_gaffmol2_vs_sybylmol2_vs_obabelpdb(charge_method=None)
