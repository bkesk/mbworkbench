#!/bin/env python3
import argparse
import logging
import hashlib
import yaml

import numpy as np
import h5py as h5

from pyscf import gto, scf, mcscf
from pyscf.shciscf import shci

# Import specific functions / classes
from my_pyscf.scf.hf_custom import custom_ROHF, hsints_to_eri

from embedding.pyscf_embedding_lib import make_embedding_H_afqmclab
from embedding.cholesky.simple_cholesky import cholesky
from embedding.cholesky.integrals.gto import GTOIntegralGenerator

# TODO: make virtual env and use current python version(s)
# TODO: define spin and charge in input file (yaml)

def read_ecp_basis(fname):
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
        ecp_basis = yaml.safe_load(f) # limits load to simple python objects, just what we want!
    
    logging.debug(f' ecp_basis = {ecp_basis} ')
    
    try:
        ecp = ecp_basis['ecp']
    except KeyError:
        ecp = None
    
    basis = _process_basis(ecp_basis['basis'])

    logging.debug(f' ecp = {ecp} ')
    logging.debug(f' basis = {basis} ')

    return ecp, basis

def read_geom(fname):
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
        input_file = yaml.safe_load(f) # limits load to simple python objects, just what we want!
    
    logging.debug(f' input_file = {input_file} ')
    
    try:
        ecp = input_file['geom']
    except KeyError:
        logging.error('could not read geometry from input file')
        exit()

    geom = input_file['geom']
    logging.debug(f' geom = {geom} ')

    print(f'\n [+] reading geometry : comment {geom["comment"]} ')
    atoms = geom["atoms"]

    return atoms


def read_embedding_params(fname):
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
        input_file = yaml.safe_load(f) # limits load to simple python objects, just what we want!
    
    logging.debug(f' input_file = {input_file} ')
    
    try:
        ecp = input_file['embedding']
    except KeyError:
        print(f'\n [+] no embedding parameters found in input file : proceding with no embedding ')
        return None

    emb = input_file['embedding']
    logging.debug(f' emb = {emb} ')

    return emb


def get_args():
    '''
    gets and parses terminal arguments

    returns:
    args - an object containing the parsed arguments - see argparse documentation

    Arguments:
    - geometry (str) : name of file where geometry is located.
    - basis (str) : name of basis to use
    - spin (int) : Nup - Ndown for the desired total spin state
    - nroots (int) : number of roots for soc-shci
    - avg_nroots (int) : number of roots for scalar state-avg. SHCISCF.
    - run_scalar (bool) : whether to run a scalar SHCISCF step before SOC-SHCI
    - verbose (int) : verbosity level
    - test (bool) : run test
    '''
    
    parser = argparse.ArgumentParser(description='Compute truncated SOC-SHCI wavefunctions for AFQMC trial wavefunctions.')
    
    # required
    parser.add_argument('--input', '-i', metavar='input file', type=str,
                        action='store',
                        dest='input_file',
                        default='input.yaml',
                        help="name of input file (YAML format)")

    parser.add_argument(metavar='spin', type=int,
                        action='store',
                        dest='spin',
                        help='spin specified as Nup - Ndown')

    parser.add_argument('-nsoc', metavar='soc_nroots', type=int,
                        action='store',
                        dest='soc_nroots',
                        default=1,
                        help='number of roots to compute with SOC-SHCI')

    parser.add_argument('-navg', metavar='avg_nroots', type=int,
                        action='store',
                        dest='avg_nroots',
                        default=1,
                        help='number of roots to use in the state-avg. scalar SHCI-Step (if enabled)')

    # optional
    #parser.add_argument('--scalar-step', action=argparse.BooleanOptionalAction)
    parser.add_argument('--no-scalar-step','-nS',action='store_false',
                        dest='run_scalar',
                        help='if set, skips the scalar SHCISCF calculation before SOC-SHCI')

    parser.add_argument('--output','-o', metavar='output', type=str,
                        action='store',
                        dest='outname',
                        default='psi_soc_shci',
                        help='base name of output file for wavefunctions.')

    parser.add_argument('--verbose', '-v', action='count',
                        default=4,
                        dest='verb',
                        help='set verbosity level, use more \'v\' chars to make output more verbose (i.e. -vv ) default is -vvvv')

    parser.add_argument('--test','-t',action='store_true',
                        dest='run_test',
                        help='if set, runs test based on Br atom')

    args = parser.parse_args()
    return args

def SOC_SHCI(mf, norb, nelec, nroots=1, maxM=1000, tol=1.0e-8, ncore_afqmc=0, afqmc_basis=None, *args, **kwargs):
    
    mc = mcscf.CASCI(mf, norb, nelec, *args, **kwargs)
    if nroots > 1:
        mc.nroots = nroots

    if afqmc_basis is None:
        afqmc_basis = mc.mo_coeff
    
    shci.write_soc_ecp_integrals_mol(mc.mol, ncas=norb, ncore=ncore_afqmc+mc.ncore, mo_coeff=afqmc_basis)
    mc.fcisolver = shci.SHCI(mf.mol, maxM, tol=tol)
    mc.fcisolver.executable = '/home/bkeskr/sources/DICE-spec/Dice/ZDice2'
    mc.fcisolver.DoRDM = False

    if nroots > 1:
        mc.fcisolver.nroots = nroots

    mc.fcisolver.mpiprefix = "mpirun -np 1"
    mc.fcisolver.stochastic = True
    mc.fcisolver.nPTiter = 0  
    mc.fcisolver.sweep_iter = [0]
    mc.fcisolver.sweep_epsilon = [1.0e-3]

    mc.fcisolver.shci_extra_keyword.append(("writeBestDeterminants", "2000"))
    mc.fcisolver.shci_extra_keyword.append(("doSOC", ""))

    return mc

def main(geom,nroots=1,avg_nroots=None,run_scalar=True,emb_params=None,**mol_kwargs):
    '''
    
    Inputs:
    - geom (list, other): a list containing the atomic labels and positions (other: any geometry format that PySCF accepts) 
    - nroots (int) : number of roots to compute in the SOC-SHCI calculation
    - avg_nroots (int) : number of roots to average over in scalar SHCISCF (defaults to nroots)
    - run_scalar (bool) : if True, run scalar SHCISCF before SOC-SHCI

    TODO:
    - add option to force re-calculation of the Cholesky vectors in GTO basis
    '''

    if emb_params['afqmc']['ncore'] != 'None':
        ncore_afqmc = int(emb_params['afqmc']['ncore'])
    else:
        ncore_afqmc = 0

    

    mol = gto.M(atom=geom, **mol_kwargs)
    M_full = mol.nao_nr()

    if emb_params['scalar_shci']['nactive'] != 'None':
        norb_scalar_shci = int(emb_params['scalar_shci']['nactive'])
        logging.debug(f'norb_scalar_shci read from input : {norb_scalar_shci}')
    else:
        norb_scalar_shci = M_full - ncore_afqmc
        logging.info(f'norb_scalar_shci not set : using full basis minus afqmc core orbitals : {norb_scalar_shci}')


    if emb_params['scalar_shci']['nelec'] != 'None':
        nelec_scalar_shci = int(emb_params['scalar_shci']['nelec'])
        logging.debug(f'nelec_scalar_shci read from input : {nelec_scalar_shci}')
    else:
        norb_scalar_shci = sum(mol.nelec) - 2*ncore_afqmc
        logging.info(f'nelec_scalar_shci not set : using all electrons minus afqmc core electrons : {nelec_scalar_shci}')

    if emb_params['soc_shci']['nactive'] != 'None':
        norb_soc_shci = int(emb_params['soc_shci']['nactive'])
        logging.debug(f'norb_soc_shci read from input : {norb_soc_shci}')
    else:
        norb_soc_shci = M_full - ncore_afqmc
        logging.info(f'norb_soc_shci not set : using full basis minus afqmc core orbitals : {norb_soc_shci}')

    if emb_params['soc_shci']['nelec'] != 'None':
        nelec_soc_shci = int(emb_params['soc_shci']['nelec'])
        logging.debug(f'nelec_soc_shci read from input : {nelec_soc_shci}')
    else:
        nelec_soc_shci = sum(mol.nelec) - 2*ncore_afqmc
        logging.info(f'nelec_soc_shci not set : using all electrons minus afqmc core electrons : {nelec_soc_shci}')

    # 1. initial SCF
    mf = scf.ROHF(mol).newton()
    mf.chkfile = 'ROHF-direct.chk'
    mf.kernel()

    afqmc_basis = mf.mo_coeff

    # 2. make H^A - just based on energy to start with
    Sao = mol.intor('int1e_ovlp')
    k = mol.intor('int1e_kin')
    v = mol.intor('int1e_nuc')
    v += mol.intor('ECPscalar')

    try:
        print(f'\n [+] attempting to read Cholesky vectors from file ... ', end='')
        with h5.File('twoBodyAO','r') as f:
            choleskyAO = f['choleskyAO'][...]
        print(f'success!')
    except OSError:
        print(f'failed! computing directly')
        delta = 1.0E-6
        gto_gen = GTOIntegralGenerator(mol)
        _, choleskyAO = cholesky(gto_gen,tol=delta)
        with h5.File('twoBodyAO','w') as f:
            f.create_dataset('choleskyAO',data=choleskyAO)
            f.create_dataset('delta',data=delta)
            # TODO: check is the molecule is the same as the one used to compute the Cholesky vectors.
            #f.create_dataset('basis_hash',data=hashlib.sha256(str(mol.basis)).digest)


    twoBody,_,oneBody,S,E_const = make_embedding_H_afqmclab(nfc=ncore_afqmc, 
                                                            oneBody=k+v,
                                                            C = afqmc_basis,
                                                            twoBody=choleskyAO,
                                                            S=Sao,
                                                            tol=1.0E-5, 
                                                            transform_only=True)
    
    nelec = mol.nelec
    nelec_fc = (nelec[0] - ncore_afqmc, nelec[1] - ncore_afqmc)
    
    # 3. rohf for H^A
    mol.nelec = nelec_fc
    
    rohf_A = custom_ROHF(mol, S, E_const, oneBody, eri=hsints_to_eri(twoBody))

    if run_scalar:

        if avg_nroots is None:
            avg_nroots = nroots
        # 4. scalar SHCISCF in lmited active space
        if avg_nroots > 1:
            # TODO: need to average over total spin states (where relevant)
            print(f'\n [+] multiple roots requested : using State-Averaging with equal weights')
            weights = np.ones(avg_nroots) / avg_nroots
            mc = shci.SHCISCF(rohf_A, norb=norb_scalar_shci, nelec=nelec_scalar_shci).state_average_(weights)
        else:
            mc = shci.SHCISCF(rohf_A, norb=norb_scalar_shci, nelec=nelec_scalar_shci)

        mc.fcisolver.mpiprefix = "mpirun -np 4"
        mc.fcisolver.stochastic = True
        mc.fcisolver.sweep_iter = [0]
        mc.fcisolver.sweep_epsilon = [1.0e-3]
        mc.kernel()

        # TODO: sort orbitals by energy!

        # 5. update basis to scalar SHCISCF basis
        rohf_A.mo_coeff = mc.mo_coeff
        dm = rohf_A.make_rdm1()
        Evar = rohf_A.energy_elec(dm)
        print(f'\n  [+] E_var of leading determinant is {Evar[0] + E_const}')

    # 6. SOC-SHCI
    mc = SOC_SHCI(rohf_A, norb_soc_shci, nelec_soc_shci, ncore_afqmc=ncore_afqmc, nroots=nroots, afqmc_basis=afqmc_basis)
    mc.kernel()

    # 7. save WFN for AFQMC
    

    # 8. save H for AFQMC

def test():
    ecp = {'Br' : '''
    Br nelec 10
    Br ul
    2       1.0000000              0.0000000
    Br S
    2      70.0242570             49.9628340
    2      31.1784120            370.0142050
    2       7.1565930             10.2414390
    Br P
    2      46.7734710             99.1122440    -198.224488    
    2      46.1841200            198.2530460     198.253046
    2      21.7138580             28.2617400    -56.5234800 
    2      20.9417920             56.6233660     56.6233660
    Br D
    2      50.6988390            -18.6058530      18.605853 
    2      50.6447640            -27.9232800     -18.615520
    2      15.4475090             -0.3796930       0.379693
    2      15.5002590             -0.7805830      -0.520389
    2       2.8003910              0.0359680      -0.035968
    2       1.0774800              0.0943970       0.062931
    Br F
    2      14.4656060             -1.0912690       0.727513
    2      21.2340650             -2.8876910      -1.443846
    '''}

    basis='ccpvtzpp'

    geometry = [['Br',(0.000,0.000,0.000)]]

    nroots = 6 # number of states for state-averaging

    avg_nroots = 3

    emb_params = {'afqmc' : {'ncore' : 2,
                             'nactive' : None},
                  'scalar_shci' : {'nactive' : 9,
                                    'nelec' : 7},
                  'soc_shci' : {'nactive' : 32,
                                'nelec' : 21}}

    main(geometry,
         nroots=nroots,
         avg_nroots=avg_nroots,
         run_scalar=True,
         emb_params=emb_params,
         basis=basis,
         ecp=ecp,
         spin=1,
         verbose=4)


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s : [%(funcName)s - %(levelname)s] : %(message)s')

    cli = get_args()
    logging.debug(f'cli args = {cli}')
    if cli.run_test:
        test()
    else:
        ecp, basis = read_ecp_basis(cli.input_file)
        geom = read_geom(cli.input_file)
        emb_params = read_embedding_params(cli.input_file)

        main(geom=geom,
            nroots=cli.soc_nroots,
            avg_nroots=cli.avg_nroots,
            run_scalar=cli.run_scalar,
            emb_params=emb_params,
            spin=cli.spin,
            basis=basis,
            ecp=ecp,
            verbose=cli.verb)
