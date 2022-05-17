
import logging
import numpy as np

import mbworkbench.lib.cholesky.cholesky as ch

def make_embedding_H(nfc=0,nactive=None,Enuc=0.0,tol=1.0e-6,C=None,twoBody=None,oneBody=None,S=None,transform_only=False):
    '''
    high level function to produce the embedding / downfolding Hamiltonian
    saves the results to files

    gets CVs from V2b_source, which chould contain the GTO basis CVs
    similarly for V1b_source

    Inputs:
    
    Outputs:
    > saves the following files

    #TODO: add option to use a list of active orbital indices - internally, we should always do this. The problem then reduces to generating this list of indices.
    '''

    # 1. read in orbitals (go ahead and remove the last ntrim orbitals)
    try:
        assert C is not None
    except AssertionError:
        logging.error('orthonormal orbital basis is required')
        return None

    if nactive is None:
        nactive = S.shape[0] - nfc
    
    C = C[:,:nfc+nactive]
    
    # 2. read CVs from file, and transform to MO basis
    Alist = ch.ao2mo_cholesky(C,twoBody)
    Ncv = twoBody.shape[0]

    # in some cases, we only want to transform CVs to the MO basis
    if transform_only:
        print('Only transforming from GTO to orthonormal basis with no additional Cholesky decomposition', flush=True)
        NcvActive = Ncv
        twoBodyActive = Alist[:,nfc:,nfc:]
    else:        
        from mbworkbench.lib.cholesky.integrals import FactoredIntegralGenerator
        print(f'Performing Cholesky decomposition within the active space num. frozen occupied={nfc}, num. of active orbitals = {nactive}', flush=True)
        V = FactoredIntegralGenerator(Alist[:,nfc:,nfc:])
        twoBodyActive = ch.cholesky(integral_generator=V,tol=tol)
        NcvActive = twoBodyActive.shape[0]
        del(V)

    print('Computing one-body embedding terms', flush=True)    
    S_MO = ch.ao2mo_mat(C,S)
    oneBody_MO = ch.ao2mo_mat(C,oneBody)

    S_active = S_MO[nfc:,nfc:]
    oneBody_active = oneBody_MO[nfc:,nfc:]

    print(f'shape of oneBody_active is {oneBody_active.shape}')
    if nfc > 0:
        oneBody_active+=get_embedding_potential_CV(nfc, C, Alist, AdagList=None,is_complex=False)
    
    # 7. compute constant Energy
    #       - E_K = trace K over core orbitals
    #       - E_V = embedding constant from V2b
    print('Computing constant energy', flush=True)
    if nfc > 0:
        E_K = 2*np.einsum('ii->',oneBody_MO[:nfc,:nfc])
        E_V = get_embedding_constant_CV(C,Alist[:,:nfc,:nfc],AdagList=None,is_complex=False)
    else:
        E_K=0.0
        E_V=0.0
    E_const = Enuc + E_K + E_V
    print(f'E_0 = Enuc + E_K + E_V = {E_const} with:\n  - Enuc = {Enuc}\n  - E_K = {E_K}\n  - E_V = {E_V}')

    return twoBodyActive,NcvActive,oneBody_active,S_active,E_const




def get_embedding_potential_CV(nfc, C, Alist, AdagList, debug=False,is_complex=True):
    '''
    Computes the embedding potential from MO basis Cholesky vectors
    
    NOTES: make cuts before calling in C, Alist, AdagList
    '''
    
    if is_complex:
        G_core = np.eye(nfc,dtype='complex128')
        # compute the direct term as G_{I L} * V_{I j k L} -> Pyscf (Chemists') notation, want (IL|jk) mo integrals
        print('[+] computing Vd ...')
        Vd = np.einsum('il,gil,gjk->jk',G_core,Alist[:,:nfc,:nfc],AdagList[:, nfc:, nfc:],optimize='greedy')

        # compute the exchange term as G_{I L} * V_{i J k L} -> Pyscf (Chemists') notation, want (iL|Jk) mo integrals
        print('[+] computing Vx ...')
        Vx = np.einsum('jl,gil,gjk->ik',G_core,Alist[:,nfc:,:nfc],AdagList[:,:nfc,nfc:],optimize='greedy')
    else:
        G_core = np.eye(nfc)
        # compute the direct term as G_{I L} * V_{I j k L} -> Pyscf (Chemists') notation, want (IL|jk) mo integrals
        print('[+] computing Vd ...')
        Vd = np.einsum('il,gil,gjk->jk',G_core,Alist[:,:nfc,:nfc],Alist[:, nfc:, nfc:],optimize='greedy')

        # compute the exchange term as G_{I L} * V_{i J k L} -> Pyscf (Chemists') notation, want (iL|Jk) mo integrals
        print('[+] computing Vx ...')
        Vx = np.einsum('jl,gil,gjk->ik',G_core,Alist[:,nfc:,:nfc],Alist[:,:nfc,nfc:],optimize='greedy')
    
    return 2*Vd - Vx


def get_embedding_constant_CV(C, Alist, AdagList, debug=False, is_complex=True):
    '''
    Computes the embedding constant from MO basis Cholesky vectors
    
    NOTES:- make cuts before calling in C, Alist, AdagList
          - no need for C, its assumed that Alist / AdagList as in correct basis
    
    Inputs:
    C - array containing just the inactive orbitals
    Alist, AdagList - restricted to frozen orbitals
    '''
    if is_complex:
        Vd = np.einsum('gii,gjj->',Alist,AdagList,optimize='greedy')
        Vx = np.einsum('gij,gji->',Alist,AdagList,optimize='greedy')
    else:
        Vd = np.einsum('gii,gjj->',Alist,Alist,optimize='greedy')
        Vx = np.einsum('gij,gji->',Alist,Alist,optimize='greedy')
    
    return 2*Vd - Vx 