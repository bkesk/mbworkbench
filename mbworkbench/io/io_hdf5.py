import logging
import h5py as h5

"""
We read in data from binary hdf5 files.
Settings/parameters should be included in the yaml-based input file.

Examples of data to read from hdf5 files:
- orbitals
- matrix elements/intergrals
- transformaiton matrices
"""

def read_orbitals(orb_string, orb_name='/scf/mo_coeff'):
    '''
    read orbitals from an HDF5 file. Orbitals should be stored based
    on PySCF's conventions: '/scf/mo_coeff' by default.
    '''
    
    temp = orb_string.split(':')
    fname = temp[0]
    if len(temp) == 2:
        orb_name = temp[1]

    try:
        with h5.File(fname,'r') as f:
            orbitals = f[orb_name][...].copy()
        return orbitals
    except Exception as e:
        logging.error(f'Unable to read orbitals: {e}')
        return None



