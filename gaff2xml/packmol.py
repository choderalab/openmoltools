import numpy as np
import shutil
import os
import mdtraj as md
from mdtraj.utils import enter_temp_directory
import tempfile
from distutils.spawn import find_executable

PACKMOL_PATH = find_executable("obabel")

HEADER_TEMPLATE = """
# Mixture 

tolerance %f
filetype pdb
output %s
add_amber_ter

"""

BOX_TEMPLATE = """
structure %s
  number %d 
  inside box 0. 0. 0. %f %f %f
end structure
"""

def pack_box(pdb_filenames, n_molecules_list, tolerance=2.0, box_size=40.):
    """Run packmol to generate a box containing a mixture of molecules.

    Parameters
    ----------
    pdb_filenames : list(str)
        List of pdb filenames for each component of mixture.
    n_molecules_list : list(int)
        The number of molecules of each mixture component.
    tolerance : float, optional, default=2.0
        The mininum spacing between molecules during packing.  In ANGSTROMS!
    box_size : float, optional, default=40.0
        The size of the box to generate.  In ANGSTROMS.

    Returns
    -------
    trj : MDTraj.Trajectory
        Single frame trajectory with mixture box.

    Notes
    -----
    Be aware that MDTraj uses nanometers internally, but packmol uses angstrom
    units.  The present function takes `tolerance` and `box_size` in 
    angstrom units, but the output trajectory will have data in nm.  
    Also note that OpenMM is pretty picky about the format of unit cell input, 
    so use the example in tests/test_packmol.py to ensure that you do the right thing.

    """    
    assert len(pdb_filenames) == len(n_molecules_list), "Must input same number of pdb filenames as num molecules"
    
    output_filename = tempfile.mktemp(suffix=".pdb")
    
    header = HEADER_TEMPLATE % (tolerance, output_filename)
    for k in range(len(pdb_filenames)):
        filename = pdb_filenames[k]
        n_molecules = n_molecules_list[k]
        header = header + BOX_TEMPLATE % (filename, n_molecules, box_size, box_size, box_size)
    
    pwd = os.getcwd()
    
    print(header)
    
    packmol_filename = "packmol_input.txt"
    packmol_filename = tempfile.mktemp(suffix=".txt")
    
    file_handle = open(packmol_filename, 'w')
    file_handle.write(header)
    file_handle.close()
    
    print(header)

    os.system("/home/kyleb/src/Software/packmol/packmol < %s" % packmol_filename)

    trj = md.load(output_filename)
    
    #Begin hack to introduce bonds for the MISSING CONECT ENTRIES THAT PACKMOL FAILS TO WRITE
    
    top, bonds = trj.top.to_dataframe()

    trj_i = [md.load(filename) for filename in pdb_filenames]
    bonds_i = [t.top.to_dataframe()[1] for t in trj_i]

    offset = 0
    bonds = []
    for i in range(len(pdb_filenames)):
        n_atoms = trj_i[i].n_atoms
        for j in range(n_molecules_list[i]):        
            bonds.extend(bonds_i[i] + offset)
            offset += n_atoms

    bonds = np.array(bonds)
    trj.top = md.Topology.from_dataframe(top, bonds)
    
    trj.unitcell_vectors = np.array([np.eye(3)]) * box_size / 10.
    
    return trj


def dipole_moment(traj, charges):
    xyz = md.compute_displacements(traj, np.array([[0, i] for i in range(traj.n_atoms)])).astype('float')
    return xyz.transpose(0, 2, 1).dot(charges)
