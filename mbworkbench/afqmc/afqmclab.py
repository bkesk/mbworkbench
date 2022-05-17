'''
Module to handle writting AFQMCLab input files.

need to write:

- model_param : Hamiltonian file
- afqmc_param : the input file for AFQMCLab
- phiT.dat : trial wavefunction
- phi.dat : initial walker
'''

import logging
import h5py as h5
from yaml import safe_load

from mbworkbench.lib import block as blk
from mbworkbench.io.afqmclab import writeModel, writeROHF

class Write_Afqmclab(blk.Block):

    def __init__(self, name):
        super().__init__(name)


    def step(self, data):
        '''
        Run the Write_Afqmclab block. Will build a settings dictionary,
        based on defaults and overwrite settings based on thie input file.
        '''

        settings = get_defaults()
        # update settings based on input file
        read_afqmc_input(data['input_file'], settings)

        try:
            assert self.status == blk.OKAY_2_RUN
        except AssertionError:
            logging.error("[dev] Attempted to run block without status OKAY_2_RUN ")
            self.status = blk.FAILED
            self.tear_down()
            return

        try:
            #TODO: add afqmc_basis to the data!
            logging.info('Setting AFQMC orbital basis')
            get_afqmc_basis(settings['__orbital_basis'],data)
            logging.info('building/write AFQMCLab model_param file')
            write_model_file(settings, data)
            logging.info('writing AFQMCLab afqmc_param file')
            write_afqmc_param(settings)
            self.status = blk.COMPLETE
            return
        except Exception as e:
            logging.error(f"Write_Afqmclab block failed to run. Exception: {e}")
            self.status = blk.FAILED
            self.tear_down()
            return

    def check_data(self, data):
        try:
            #TODO : add checks on AFQMC input data
            self.status = blk.OKAY_2_RUN
        except Exception as e:
            logging.error(f" {self.name} : Workflow data is not correct. Exception: {e}")
            self.status = blk.FAILED

    def write_to_chk(self, data, chkfile=None):
        '''
        Write data for Block to chkfile.

        Dev Note: specific Block types are responsible for serializing/deserializing
        any objects that they have added to 'data'.
        '''
        
        try:
            assert chkfile is not None
            pass
        except AssertionError:
            logging.error("no checkpoint file specified")

    def load_from_chk(self, data, chkfile=None):
        '''
        Read Block's data from chkfile, and add it to data

        Dev Note: specific Block types are responsible for serializing/deserializing
        any objects that they have added to 'data'.
        '''
        try:
            assert chkfile is not None
        except AssertionError:
            logging.error("no checkpoint file specified")

        try:
            with h5.File(chkfile, 'r') as f:
                generic = f[self.name + '/generic' ][...]
        except ValueError as e:
            logging(f"read from checkpoint failed. Excetion {e}")



def read_afqmc_input(fname, settings):
    '''
    read "write_afqmclab" parameters from yaml file, and add them to
    settings. This will overwrite existing settings.


    ```
    write_afqmclab:
      basis: [source]:[name]
      integrals: (molecule|file)
      write_wfns:
        - scf
        - ...
      run_params:
        - dt: 
        - .... 

      ...

    ```

    '''

    with open(fname,'r') as f:
        infile = safe_load(f) # limits load to simple python objects, just what we want!
    try:
        assert 'write_afqmclab' in infile.keys()

        afqmclab = infile['write_afqmclab']
        logging.info(f' found "write_afqmclab" block in input file ')
        logging.debug(f' input file: write_afqmclab  = {afqmclab} ')

        if 'basis' in afqmclab.keys():
            logging.debug(' found "basis" block in "write_afqmclab" input file block')
            settings['__orbital_basis'] = afqmclab['basis']

        if 'integrals' in afqmclab.keys():
            logging.debug(' found "integrals" block in "write_afqmclab" input file block')
            settings['__integrals'] = afqmclab['integrals']

        if 'write_wfns' in afqmclab.keys():
            logging.debug(' found "write_wfns" block in "write_afqmclab" input file block')
            settings['__write_wfns'] = afqmclab['write_wfns']

    except AssertionError:
        logging.info(f' no "write_afqmclab" block in input file ')
        return None, {}


def get_defaults():
    '''
    make default settings dictionary. Keys beginning with a double-underscore
    are settings related to the Block instances.
    Keys beginning with a single underscore will be computed internally
    based on other values.
    Keys that have no leading underscore are used directly.
    '''

    from random import randint

    default_settings = {'__model_name' : 'model_param',
                        '__orbital_basis' : 'scf',
                         '__integrals' : 'molecule',
                         '__write_wfns' : ['scf'],
                        'dt' : 0.01,
                        '_thermalization_size' : 1000, # units : steps
                        '_block_to_write' : 100, # units: _block_size 
                        '_block_size' : 20, # units _block_unit
                        '_block_unit' : 5, # units : steps
                        '_walkers_per_thread' : 10,
                        'force_type' : 'dynamicForce',
                        'force_cap' : 5,
                        'trial_wavefuntion' : 'readFromFile',
                        'initial_walkers' : 'sampleFromPhiT',
                        'mgs_interval' : 5, # units : steps
                        'pop_control_interval' : 10, # units : steps
                        'initial_pop_control_length' : 1000, # units : steps
                        'log_energy_cap' : 2.3,
                        'trial_energy' : 0.0, # should always be set!
                        'initial_hs_background' : 'EstimateFromPhiTWalker',
                        'adjust_trial_energy_interval' : 5, # units : steps
                        'adjust_trial_energy_legnth' : 800, # units : steps
                        'hs_background_growth_interval' : 5, # units : steps
                        'hs_background_growth_length' : 1000,
                        'rng_seed' : randint(000000000,999999999) }

    return default_settings

def write_afqmc_param(settings):
    '''
    write afqmc_param based on 'settings' dictionary
    '''
    
    with open('afqmc_param','w') as f:
        f.write('''
dt                                  {settings['dt']}
thermalSize                         {settings['_thermalization_size']}
writeNumber                         {settings['_block_to_write']}
measureNumberPerWrite               {settings['_block_size']}
measureSkipStep                     {settings['_block_unit']}

walkerSizePerThread                 {settings['_walkers_per_thread']}

forceType                           {settings['force_type']}
forceCap                            {settings['force_cap']}
initialPhiTFlag                     {settings['trial_wavefuntion']}
initialWalkerFlag                   {settings['initial_walkers']}

mgsStep                             {settings['mgs_interval']}
popControlStep                      {settings['pop_control_interval']}
initPopControlMaxSize               {settings['initial_pop_control_length']}

logEnergyCap                        {settings['log_energy_cap']}

ET                                  {settings['trial_energy']}
backGroundInit                      {settings['initial_hs_background']}
ETAdjustStep                        {settings['adjust_trial_energy_interval']}
ETAdjustMaxSize                     {settings['adjust_trial_energy_legnth']}
ETAndBackGroundGrowthEstimateStep   {settings['hs_background_growth_interval']}
ETAndBackGroundGrowthEstimateMaxSize {settings['hs_background_growth_length']}

seed                                {settings['rng_seed']} # ex: 905856018
''')


def write_model_file(settings, data):
    '''
    high-level function to write model_param file for AFQMCLab
    '''
    from mbworkbench.lib.embedding import make_embedding_H

    get_nelec(settings['__integrals'], data)
    get_one_body(settings['__integrals'], data)
    get_two_body(settings['__integrals'], data)
    
    
    # use the embedding lib convenience function to transform basis.
    twoBody,numCholeskyActive,oneBody,S,E_const = make_embedding_H(tol=1.0e-6,
                                                                   C=data['afqmc_basis'],
                                                                   S=data['s'],
                                                                   twoBody=data['cholesky_vecs'],
                                                                   oneBody=data['one_body'])

    writeModel(nElec=data['nelec'],
               oneBody=oneBody,
               twoBody=twoBody,
               fname=settings['__model_name'])


def get_afqmc_basis(source, data):
    '''
    '''
    if source == 'scf':
        # get from scf block orbitals
        data['afqmc_basis'] = data['scf/mo_coeff']
    else:
        raise ValueError("Unkown source for nelec")

def get_nelec(source, data):
    '''
    get nelec based on the requested source using workflow data
    '''
    if source == 'molecule':
        nelec_from_mol(data)
    elif source == 'file':
        raise NotImplementedError
    else:
        raise ValueError("Unkown source for nelec")


def nelec_from_mol(data):
    '''
    get nelec from molecule object in workflow data
    '''
    try:
        assert 'mol' in data.keys()
        data['nelec'] = data['mol'].nelec
    except AssertionError:
        raise ValueError("No molecule object in workflow data.")


def get_one_body(source, data):
    '''
    get one_body integrals based on the requested source using workflow data
    '''
    if source == 'molecule':
        get_one_body_from_mol(data)
    elif source == 'file':
        raise NotImplementedError
    else:
        raise ValueError("Unkown source for one body integrals")


def get_one_body_from_mol(data):
    '''
    '''
    mol = data['mol']

    s = mol.intor('int1e_ovlp')
    data['s'] = s

    h1 = mol.intor('int1e_kin')
    h1 += mol.intor('int1e_nuc')

    if mol.ecp:
        h1 += mol.intor('ECPscalar')

    data['one_body'] = h1


def get_two_body(source, data):
    '''
    get two_body integrals based on the requested source using workflow data
    
    TODO: add option for density fitting
    '''
    if source == 'molecule':
        get_two_body_from_mol(data)
    elif source == 'file':
        raise NotImplementedError
    else:
        raise ValueError("Unkown source for one body integrals")


def get_two_body_from_mol(data):
    '''    
    '''
    from mbworkbench.lib.cholesky.cholesky import cholesky
    from mbworkbench.lib.cholesky.integrals import GTOIntegralGenerator

    mol = data['mol']
    verb = mol.verbose
    mol.verbose = 4
    gto_gen = GTOIntegralGenerator(mol)
    choleskyAO = cholesky(gto_gen,tol=1.0E-6)
    data['cholesky_vecs'] = choleskyAO
    mol.verbose = verb

