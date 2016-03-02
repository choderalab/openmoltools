#!/usr/bin/env python

"""Test functions in openmoltools.schrodinger."""

import unittest

from openmoltools.schrodinger import *


@unittest.skipIf(not is_schrodinger_suite_installed(), "This test requires Schrodinger's suite")
def test_structconvert():
    """Test run_structconvert() function."""
    benzene_path = utils.get_data_filename("chemicals/benzene/benzene.pdb")

    def collect_lines(pdb_path):
        """Collect all HETATM and CONNECT lines in the pdb."""
        all_lines = []
        with open(pdb_path, 'r') as f:
            for line in f:
                field = line[:6]
                n_atom = line[10:11]
                if n_atom == '4' or n_atom == '9':
                    continue  # Discard atoms 4 and 9 which have a -0.000
                if field == 'HETATM' or field == 'CONECT':
                    all_lines.append(line.strip())
        return all_lines

    with mdtraj.utils.enter_temp_directory():
        # Convert from pdb to mol2 and back
        run_structconvert(benzene_path, 'benzene.mol2')
        run_structconvert('benzene.mol2', 'benzene.pdb')
        new_lines = collect_lines('benzene.pdb')

    # The new pdb should be equal to the old one
    original_lines = collect_lines(benzene_path)
    assert original_lines == new_lines


@unittest.skipIf(not is_schrodinger_suite_installed(), "This test requires Schrodinger's suite")
def test_proplister():
    """Test run_proplister() function."""
    benzene_path = utils.get_data_filename("chemicals/benzene/benzene.sdf")
    properties = run_proplister(benzene_path)

    assert len(properties) == 1
    assert len(properties[0]) == 23

    # Test subset of properties
    expected = {'i_sd_PUBCHEM_COMPOUND_CID': '241',
                'r_sd_PUBCHEM_CONFORMER_RMSD': '0.4',
                'i_sd_PUBCHEM_CONFORMER_DIVERSEORDER': '1',
                's_sd_PUBCHEM_MMFF94_PARTIAL_CHARGES': ('12\n1 -0.15\n10 0.15\n11 0.15\n'
                                                        '12 0.15\n2 -0.15\n3 -0.15\n4 -0.15\n'
                                                        '5 -0.15\n6 -0.15\n7 0.15\n8 0.15\n'
                                                        '9 0.15')
                }
    assert set(expected.items()) < set(properties[0].items())


@unittest.skipIf(not is_schrodinger_suite_installed(), "This test requires Schrodinger's suite")
def test_epik_maesubset_autoconvert():
    """Test run_epik and run_maesubset functions and autoconvert_maestro decorator."""
    imatinib_path = utils.get_data_filename("chemicals/imatinib/imatinib.sdf")

    with mdtraj.utils.enter_temp_directory():
        run_structconvert(imatinib_path, 'imatinib.mae')
        run_epik('imatinib.mae', 'imatinib-epik.mae', ph=7.0)
        run_maesubset('imatinib-epik.mae', 'imatinib02.mae', range=[0, 2])
        run_structconvert('imatinib02.mae', 'imatinib02.sdf')

        # The 4 lines above should be equivalent to
        run_epik(imatinib_path, 'imatinib-auto02.sdf', ph=7.0, tautomerize=True,
                 extract_range=[0, 2])

        # Check that results contain indeed 2 molecules
        assert len(run_proplister('imatinib02.sdf')) == 2  # 2 molecules
        assert len(run_proplister('imatinib-auto02.sdf')) == 2

        # Check that results are identical
        with open('imatinib02.sdf', 'r') as f:
            lines = f.readlines()
        with open('imatinib-auto02.sdf', 'r') as f:
            assert f.readlines() == lines
