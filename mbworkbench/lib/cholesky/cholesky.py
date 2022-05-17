import numpy as np

from pyscf  import gto
from pyscf.gto import getints_by_shell

from .integrals import IntegralGenerator


def dampedPrescreenCond(diag, vmax, delta, s=None):
    '''
  
    Here, we evaluate the damped presreening condition as per:
    J. Comput. Chem., 31: 224-247. 2010. doi:10.1002/jcc.21318

    Also see J. Chem. Phys. 118, 9481 (2003)

    This function will return a prescreened version of the diagonal:

    a diagonal element, diag_(mu), will be set to zero if:

    s* sqrt( diag_(mu) * vmax ) <= delta

    where mu = (i,l) is a pair index, s is a numerical parameter, v is a Cholesky vector, 
    delta is the cholesky threshold and vmax is the maxium value on the diagonal

    s is the damping factor and typically:
    s = max([delta*1E9, 1.0])
    '''

    if s is None:
        s = max([delta*1E9, 1.0])

    negative = np.less(diag, 0.0) 
    diag[negative] = 0.0

    sDeltaSqr=(delta/s)*(delta/s)

    # the actual damped prescreening
    toScreen = np.less(diag*vmax, sDeltaSqr)
    diag[toScreen] = 0.0
    return diag, toScreen

def cholesky(integral_generator=None,tol=1.0E-8,prescreen=True,debug=False,max_cv=None):

    if not isinstance(integral_generator, IntegralGenerator): 
        raise TypeError('Invalid integral generator, must have base class IntegralGenerator')
    
    nbasis = integral_generator.nbasis

    delCol = np.zeros((nbasis,nbasis),dtype=bool)
    choleskyNum = 0
    
    if max_cv:
        choleskyNumGuess = max_cv
    else:
        choleskyNumGuess = 10*nbasis
    cholesky_vec = np.zeros((choleskyNumGuess,nbasis,nbasis))
    
    Vdiag = integral_generator.diagonal()

    if debug:
        print("Initial Vdiag: ", Vdiag)
        verb = 5
    else:
        verb = 4

    # np.argmax returns the 'flattened' max index -> need to 'unflatten'
    unflatten = lambda x : (x // nbasis, x % nbasis)
    
    # Note: screening should work still as long as we keep shape of Vdiag consistent!
    if prescreen: # zero small diagonal matrix elements - see J. Chem. Phys. 118, 9481 (2003)
        if debug:
            imax = np.argmax(Vdiag)
            print("imax , unflatten(imax) = ", imax, unflatten(imax))
        imax = np.argmax(Vdiag); vmax = Vdiag[unflatten(imax)]
        toScreen = np.less(Vdiag, tol*tol/vmax)
        Vdiag[toScreen] = 0.0
    
    while True:
        imax = np.argmax(Vdiag)
        vmax = Vdiag[unflatten(imax)]
        print( "Inside modified Cholesky {:<9} {:26.18e}".format(choleskyNum, vmax), flush=True )
        if(vmax<tol or choleskyNum==nbasis*nbasis or choleskyNum >= choleskyNumGuess):
            print( "Number of Cholesky fields is {:9}".format(choleskyNum) )
            print('\n')
            break
        else:
            if debug:
                print("\n*** getCholesky_OnTheFly: debugging info*** \n")
                print("imax = ", imax, " (i,l) ", (imax // nbasis, imax % nbasis))
                print("vmax = ", vmax)

            Vrow = integral_generator.row(imax, Alist=cholesky_vec[:choleskyNum,:,:])
            oneVec = Vrow/np.sqrt(vmax)
            if prescreen:
                oneVec[delCol]=0.0
            cholesky_vec[choleskyNum] = oneVec                
            if debug:
                print("Vrow", Vrow)
            choleskyNum+=1
            Vdiag -= oneVec**2
            if prescreen:
                Vdiag, removed = dampedPrescreenCond(Vdiag, vmax, tol)
                delCol[removed] = True 
            if debug:
                print("oneVec: ", oneVec)
                print("oneVec**2: ", oneVec**2)
                print("Vdiag: ", Vdiag)
                print("\n*** *** *** ***\n")

    return cholesky_vec[:choleskyNum,:,:]


#TODO: move these to their own submodule in mbworkbench/lib
def ao2mo_cholesky(C,choleskyVecAO,verb=False):
    '''
    Transforms the GTO basis Cholesky vectors to the MO basis
    
    Inputs:
       C - coefficient matrix which specifies the desired MOs in terms of the GTO basis funcs
             (index conventions: C_{mu i} mu - GTO index, i - MO index)
       choleskyVecAO - numpy array containing the Cholesky vectors represented in the GTO basis

    index convention for CVs: choleskyVecAO[gamma, mu, nu]
                with gamma - Cholesky vector index
                     mu,nu - GTO indices 
           * similar for MO basis mu,nu -> i,l

    Returns:
       chleskyVecMO - numpy array containing the Cholesky vectros represented in the MO basis
    '''
    ncv = choleskyVecAO.shape[0]
    MA = C.shape[1]
    Cdag = C.conj().T # for readability below!
    choleskyVecMO = np.zeros((ncv,MA,MA))
    for i in range(ncv):
        if verb:
            print(f'transforming vector {i}')
            if i % 100 == 0:
                print('',end='',flush=True)
        temp = np.matmul(Cdag,choleskyVecAO[i,:,:])
        choleskyVecMO[i,:,:] = np.matmul(temp,C)
    return choleskyVecMO

def ao2mo_mat(C, mat):
    '''
    Transforms the GTO basis matrix to the MO basis
    
    Inputs:
       C - coefficient matrix which specifies the desired MOs in terms of the GTO basis funcs
             (index conventions: C_{mu i} mu - GTO index, i - MO index)
       mat - numpy matrix (num GTO x num GTO) represented in the GTO basis

    Returns:
       matMO - matrix in MO basis
    '''
    #matMO = np.einsum('im,mn,nl->il', C.conj().T,mat,C)
    matMO = np.matmul(C.conj().T,mat)
    matMO = np.matmul(matMO,C)
    return matMO

    