from nose.plugins.attrib import attr
from unittest import skipIf
import tempfile
import os
from gaff2xml import utils

def _drug_tester(n_molecules=3404, charge_method="bcc"):
    """Helper function for running various versions of the drug parameterization benchmark."""
    assert n_molecules <= 3404, "The maximum number of molecules is 3404."
    assert charge_method in ["bcc", None], "Charge method must be either None or 'bcc'"

    path = tempfile.mkdtemp()
    database_filename = utils.get_data_filename("chemicals/drugs/Zdd.mol2.gz")
    cmd = "gunzip -c %s > %s/Zdd.mol2" % (database_filename, path)
    os.system(cmd)
    cmd = """awk '/MOLECULE/{close(x);x="%s/molecule_"i++".mol2"}{print > x}' %s/Zdd.mol2""" % (path, path)
    os.system(cmd)

    for k in range(n_molecules):
        molecule_name = "molecule_%d" % k
        mol2_filename = "%s/%s.mol2" % (path, molecule_name)
        cmd = """sed -i "s/<0>/LIG/" %s""" % mol2_filename
        os.system(cmd)  # Have to remove the <0> because it leads to invalid XML in the forcefield files.
        with utils.enter_temp_directory():
            yield utils.tag_description(lambda : utils.test_molecule("LIG", mol2_filename, charge_method=CHARGE_METHOD), "Testing drugs %s with charge method %s" % (molecule_name, CHARGE_METHOD))

# This version is too slow to run during the conda-build post-test, on jenkins, or on travis.
# In fact, Travis is actually too slow to do even a single AM1-BCC calculation
@attr('slow')
def test_drugs_all():
    _drug_tester()

def test_drugs_fast():
    _drug_tester(n_molecules=20, charge_method=None)
