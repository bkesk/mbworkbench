import random
import numpy as np
import h5py as h5

from mbworkbench.lib.linalg import check_hermitian

#TODO: This submodule should exist under the afqmc submodule

afqmclab_complex = np.dtype([("real","float64"),("imag","float64")])

def np_comp_to_afqmclab_comp(Mat):
    '''
    convert from default complex type to the one used in AFQMCLAb.
    This is important when writting to an hdf5 file.
    '''

    MatList = []
    dims = Mat.shape
    if len(dims) == 2:
        for j in range(dims[0]):
            row = []
            for k in range(dims[1]):
                row.append((Mat[j,k].real, Mat[j,k].imag))
            MatList.append(row)
    elif len(dims) ==1:
        for j in range(dims[0]):
            MatList.append((Mat[j].real, Mat[j].imag))
    return np.array(MatList,dtype=afqmclab_complex)


def writeModel(nElec,oneBody,twoBody,fname='model_param'):
    ''' Export Hamiltonian to the disk.'''
    
    oneBodyAdj = oneBody.copy()
    for i in range(twoBody.shape[0]):
        oneBodyAdj += (-0.5)*np.dot( twoBody[i], twoBody[i] )

    print("Testing if input oneBody is Hermitian ... ")
    check_hermitian(oneBody,tol=1.0E-8)
    print("Testing if adjusted oneBody is Hermitian ... ")
    check_hermitian(oneBodyAdj,tol=1.0E-8)

    choleskyN = twoBody.shape[0]
    basisN = twoBody.shape[1]

    f = h5.File(fname, "w")
    f.create_dataset("L",               (1,),                    data=[basisN],             dtype='int')
    f.create_dataset("Nup",             (1,),                    data=[nElec[0]],           dtype='int')
    f.create_dataset("Ndn",             (1,),                    data=[nElec[1]],           dtype='int')
    f.create_dataset("choleskyNumber",  (1,),                    data=[choleskyN],          dtype='int')
    f.create_dataset("t",               (basisN**2,),            data=oneBody.ravel(),      dtype='float64')
    f.create_dataset("K",               (basisN**2,),            data=oneBodyAdj.ravel(),   dtype='float64')
    f.create_dataset("choleskyVecs",    (choleskyN*basisN**2,),  data=twoBody.ravel(),      dtype='float64')
    f.create_dataset("choleskyBg",      (choleskyN,),            data=np.zeros(choleskyN),  dtype='float64')
    f.close()


def writeModelGHF(nElec,oneBodyIn,twoBody,is_complex=False,fname='model_param'):
    '''Write Hamiltonian to disk.'''
    
    # IMPORTANT! need to transpose the one body term in order to have correct 
    #   Hamiltonian in AFQMCLab
    oneBody = oneBodyIn.T 

    oneBodyAdj = oneBody.copy()

    for i in range(twoBody.shape[0]):
        oneBodyAdj += (-0.5)*np.dot( twoBody[i], twoBody[i] )

    choleskyN = twoBody.shape[0]
    basisN = twoBody.shape[1]

    print("Testing if input oneBody is Hermitian ... ")
    check_hermitian(oneBodyIn,tol=1.0E-8)
    print("Testing if adjusted oneBody is Hermitian ... ")
    check_hermitian(oneBodyAdj,tol=1.0E-8)

    f = h5.File(fname, "w")
    f.create_dataset("L",               (1,),                    data=[basisN],             dtype='int')
    f.create_dataset("N",             (1,),                    data=[nElec],           dtype='int')
    f.create_dataset("choleskyNumber",  (1,),                    data=[choleskyN],          dtype='int')
    if is_complex:
        f.create_dataset("t",               (basisN**2,),            data=np_comp_to_afqmclab_comp(oneBody).ravel(),      dtype=afqmclab_complex)
        f.create_dataset("K",               (basisN**2,),            data=np_comp_to_afqmclab_comp(oneBodyAdj).ravel(),   dtype=afqmclab_complex)
        f.create_dataset("choleskyBg",      (choleskyN,),            data=np_comp_to_afqmclab_comp(np.zeros(choleskyN)),  dtype=afqmclab_complex)
    else:
        f.create_dataset("t",               (basisN**2,),            data=oneBody.ravel(),      dtype='float64')
        f.create_dataset("K",               (basisN**2,),            data=oneBodyAdj.ravel(),   dtype='float64')
        f.create_dataset("choleskyBg",      (choleskyN,),            data=np.zeros(choleskyN),  dtype='float64')
    f.create_dataset("choleskyVecs",    (choleskyN*basisN**2,),  data=twoBody.ravel(),      dtype='float64')

    f.create_dataset("vecComplex",        (basisN**2,),         data=np_comp_to_afqmclab_comp(np.zeros(basisN**2)), dtype=afqmclab_complex)

    f.close()


def writeSD(wf, filename, N):
    ''' Write a generic determinant (GHF) to filename.'''

    f = open(filename, 'w')
    f.write('{:26.18e} {:26.18e} \n'.format(0.0,0.0)) # this is logw
    
    f.write('{:26d} \n'.format(2))
    f.write('{:26d} {:26d} \n'.format(wf.shape[0],N))
    
    if np.iscomplexobj(wf):
        print_form = lambda x : '{:26.18e} {:26.18e} \n'.format(x.real,x.imag)
    else:
        print_form = lambda x : '{:26.18e} {:26.18e} \n'.format(x,0.0)
    
    for i in range(N):
        for j in range(wf.shape[0]):
                f.write(print_form(wf[j,i]))
    f.close()


def writeUHF(wf, filename, nElec, noise=0.0):
    ''' Write a spin-unrestricted determinant to filename.'''
    if np.iscomplex(wf[0]).any() or np.iscomplex(wf[1]).any():
        print('Detected complex-valued orbitals: printing as complex')
        writeUHF_complex(wf, filename, nElec, noise)
    else:
        print('Detected real-valued orbitals: printing as real')
        writeUHF_real(wf, filename, nElec, noise)


def writeUHF_real(wf, filename, nElec, noise=0.0):
    ''' Write a spin-unrestricted determinant to filename.'''
    f = open(filename, 'w')
    f.write('{:26.18e} {:26.18e} \n'.format(0.0,0.0))

    for s in 0,1:
        f.write('{:26d} \n'.format(2))
        f.write('{:26d} {:26d} \n'.format(wf[s].shape[0],nElec[s]))
        for i in range(nElec[s]):
            for j in range(wf[s].shape[0]):
                f.write( '{:26.18e} {:26.18e} \n'.format( wf[s][j,i]+noise*random.random(),0.0 ) )
    f.close()


def writeUHF_complex(wf, filename, nElec, noise=0.0):
    ''' Write a spin-unrestricted determinant to filename.'''
    f = open(filename, 'w')
    f.write('{:26.18e} {:26.18e} \n'.format(0.0,0.0))

    for s in 0,1:
        f.write('{:26d} \n'.format(2))
        f.write('{:26d} {:26d} \n'.format(wf[s].shape[0],nElec[s]))
        for i in range(nElec[s]):
            for j in range(wf[s].shape[0]):
                f.write( '{:26.18e} {:26.18e} \n'.format( wf[s][j,i].real + noise*random.random(), wf[s][j,i].imag ) )
    f.close()


def readUHF(filename):
    ''' Read a spin-unrestricted determinant to filename.'''
    f = open(filename, 'r')

    logw = f.readline()
    
    wf_s = []
    nelec = []
    for s in 0,1:
        Ndim = f.readline()
        temp = f.readline().split()
        M = eval(temp[0])
        N = eval(temp[1])
        nelec.append(N)
        wf = np.zeros((M,N))
        for i in range(N):
            for j in range(M):
                temp = f.readline().split()
                real = eval(temp[0])
                imag = eval(temp[1]) # not using this for now.
                wf[j,i] = real
        wf_s.append(wf)
    f.close()

    wf = np.zeros((2,M,max(nelec)))
    for s in range(2):
        wf[s,:,:nelec[s]] = wf_s[s]

    return wf_s, nelec, M


def writeROHF(wf, filename, nElec, noise=0.0, is_complex=False):
    ''' Write a spin-unrestricted (but identical-spin) determinant to filename.'''
    f = open(filename, 'w')
    f.write('{:26.18e} {:26.18e} \n'.format(0.0,0.0))

    for s in 0,1:
        f.write('{:26d} \n'.format(2))
        f.write('{:26d} {:26d} \n'.format(wf.shape[0],nElec[s]))
        for i in range(nElec[s]):
            for j in range(wf.shape[0]):
                if is_complex:
                    f.write( '{:26.18e} {:26.18e} \n'.format( wf[j,i].real+noise*random.random(),wf[j,i].imag+noise*random.random() ) )
                else:
                    f.write( '{:26.18e} {:26.18e} \n'.format( wf[j,i]+noise*random.random(),0.0 ) )

    f.close()


