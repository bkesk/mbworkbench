import logging

from yaml import safe_load
from pyscf import gto

"""
Module for parsing yaml-based input file.

Only parameters / settings should be included here and no large data arrays.

Examples of parameters to include in the input file:
- molecular geometry
- ecp/basis
- embedding settings ncore, nactive, etc.
- spin
- charge
- symmetry

"""

def read_input(fname):
    '''
    High-level function to read yaml input files
    '''
    raise NotImplementedError


def mol_from_input(fname):
    '''
    generate PySCF Mole instance based on input file.
    '''

    mol_params = _read_mol_params(fname)
    ecp, basis = _read_ecp_basis(fname)
    geom = _read_geom(fname)

    if geom is None:
        return None

    return gto.M(atom=geom, basis=basis, ecp=ecp, parse_arg=False, **mol_params)

def _read_ecp_basis(fname):
    '''
    read ecp and basis in NWChem format from yaml file

    Format details for a system with atomic species 'atom1', 'atom2', ... :
    
    ```yaml
    basis:
      [atom1]:
        pyscf_lib: True # basis included in pyscf's library? True/False
        data: ccpvdz #pyscf basis name!
      [atom2]:
        pyscf_lib: False
        data: |
          [atom2] S
              exp1     c1
              exp2     c2
              ... basis set data in NWChem format ...
      ... entry for each atom in system (each must have an entry here)
    ecp:  # optional section.
      [atom1]: |
        [atom1] nelec [n1]
        ... ecp data ...
        ... ecp data ...
      [atom2]: |
        [atom2] nelec [n2]
        ... ecp data ...
        ... ecp data ...
      ... entry for each atom with ecp
    ```

    '''

    def _process_basis(input_basis):
        '''
        process basis into format that PySCF can directly use.
        '''
        pyscf_basis = {}
        for key in input_basis:
            if input_basis[key]['pyscf_lib'] == False:
                logging.debug(f'parsing basis data for {key}')
                parsed_basis = gto.parse(input_basis[key]['data'])
                pyscf_basis[key] = parsed_basis
            else:
                logging.debug(f'use pyscf library basis {input_basis[key]["data"]} for {key}')
                pyscf_basis[key] = input_basis[key]['data']
        return pyscf_basis

    with open(fname,'r') as f:
        ecp_basis = safe_load(f) # limits load to simple python objects, just what we want!
    
    logging.debug(f' ecp_basis = {ecp_basis} ')
    
    try:
        ecp = ecp_basis['ecp']
    except KeyError:
        ecp = None
    
    basis = _process_basis(ecp_basis['basis'])

    logging.debug(f' ecp = {ecp} ')
    logging.debug(f' basis = {basis} ')

    return ecp, basis

def _read_geom(fname):
    '''
    read geometry in xyz format from yaml file

    Format details for a system with atomic species 'atom1', 'atom2', ... :
    
    ```yaml
    geom:
      comment: a description of the geometry
      atoms: |
        [atom1] [x coord.] [y coord.] [z coord.]
        [atom2] [x coord.] [y coord.] [z coord.]
        ... entry for each atom in system
    ```

    '''

    with open(fname,'r') as f:
        input_file = safe_load(f) # limits load to simple python objects, just what we want!
    
    logging.debug(f' input_file = {input_file} ')
    
    try:
        geom = input_file['geom']
        logging.debug(f' geom = {geom} ')

        if "comment" in geom.keys():
            print(f'\n [+] geometry comment: {geom["comment"]} ')
        
        
        assert "atoms" in geom.keys()

        atoms = geom["atoms"]
        return atoms

    except AssertionError:
        logging.error('Input file does not contain atomic positions under "geom:"')
        return None
    except KeyError:
        logging.error('could not read geometry from input file')
        return None

def _read_mol_params(fname):
    '''
    read molecule parameters from yaml file.
    All parameters are optional!

    Format details:

    ```yaml
    molecule:
      spin: [int Nup - Ndn]
      charge: [int total charge]
      symmetry: [str for desired symmetry]
    ```

    note: a value of None will be interpreted as 0 for spin, charge,
         bit
    '''

    with open(fname,'r') as f:
        input_file = safe_load(f) # limits load to simple python objects, just what we want!
    
    logging.debug(f' input_file = {input_file} ')
    
    try:
        # NOTE: this may not give us a dictionary, need to check!
        mol_params = input_file['molecule']
        logging.debug(f' mol_params = {mol_params} ')
        return mol_params
    except KeyError as e:
        print(f'\n [+] no molecule parameters found in input file : attempting to proceed with PySCF defaults')
        logging.debug(f'could not read molecule parameters: {e}')
        return None

def _read_embedding_params(fname):
    '''
    read embedding parameters from yaml file

    Format details:

    ```yaml
    embedding:
      afqmc:
        ncore: None
        nactive: None
        nelec: None
      scalar_shci:
        ncore: None
        nactive: None
        nelec: None
      soc_shci:
        ncore: None
        nactive: None
        nelec: None
    ```

    note: a value of None will be interpreted as 0 for 'ncore',
          or all orbitals active for 'nactive', and all electrons for 'nelec'
    '''

    with open(fname,'r') as f:
        input_file = safe_load(f) # limits load to simple python objects, just what we want!
    
    logging.debug(f' input_file = {input_file} ')
    
    try:
        emb = input_file['embedding']
        logging.debug(f' emb = {emb} ')
        return emb
    except KeyError as e:
        print(f'\n [+] no embedding parameters found in input file : proceding with no embedding ')
        logging.debug(f'could not read embedding paramers: {e}')
        return None


def read_workflow_def(fname):
    '''
    Read high-level workflow definition from input file.
    
    returns list of workflow block names
    '''
    with open(fname,'r') as f:
        input_file = safe_load(f)

    try:
        wkflw = input_file['workflow']
        logging.debug(f' workflow = {wkflw} ')
        return wkflw
    except KeyError as e:
        print(f'\n [+] no workflow definition found in input file')
        logging.debug(f'could not read workflow definition: {e}')
        return None
