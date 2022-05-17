'''
Linear algebra utilities and numerical checks
'''
import numpy as np

from mbworkbench.lib.colors import bcolors


def check_hermitian(mat,tol=1.0E-8,anti=False):
    '''
    Check if inpute matrix is Hermitian (or anti-Hermitian).

    Inputs:
    - mat (np.array) : matrix to test.
    - tol (float) : numerical relative tolerance for test
    - anti (bool) : if True, test for anti-hermiticity 

    Based on the test built into AFQMCLAB.
    '''

    rtol = tol

    M = mat.shape[0]
    
    error = 0.0
    norm = 0.0

    matDagger = mat.conj().T

    if anti:
        # check for anti-Hermiticity instead
        factor = -1
    else:
        factor = 1

    for j in range(M):
        for i in range(j, M):
            error+=np.abs(mat[i,j] - factor*matDagger[i,j])
            norm+=np.abs(mat[i,j])
    
    norm = norm/(M*M) 
    error = error/(M*M)

    if error/norm > rtol:
        print(f" Test Hermitian [{bcolors.FAIL}!!failed!!{bcolors.ENDC}] \n   -  summed error={error}\n   -  average matrix value={norm}\n   -  (summed err. / avg matrix value) ={error/norm}\n  - 0.5M^2 = {0.5*M*M}  ")
    else:
        print(f" Test Hermitian [{bcolors.OKGREEN}pass{bcolors.ENDC}] \n   -  summed error={error}\n   -  average matrix value={norm}\n   -  (summed err. / avg matrix value) ={error/norm}\n  - 0.5M^2 = {0.5*M*M}  ")


def check_no_lindep(S,tol=1.0E-6):
    '''
    Check for linear depedence based on the overlap matrix.

    Inputs:
    - S (np.array) : overlap matrix
    - tol (float) : numerical tolerance. In this case, it is the minimum singular value
                      which will be considered as "no linear depenence".
    '''

    _,svals,_ = np.linalg.svd(S, full_matrices=True, hermitian=True)

    min_s = svals[-1]
    print(f"  [+] check_noLinDep: minimum SVD eigenvalue of overlapt matrix is {min_s}")

    if min_s < tol:
        print(f"  [!!!] WARNING!! linear dependence detected {min_s} < tolerance={tol}")
        return False
    else:
        print(f"  [+] No linear dependence detected")
        return True
    

