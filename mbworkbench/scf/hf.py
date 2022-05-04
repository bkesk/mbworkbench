import logging
import argparse
import h5py as h5

from pyscf.gto import Mole

from mbworkbench.io.input_file import mol_from_input
from mbworkbench.scf.io import read_scf_input

from mbworkbench.lib import block as blk

'''
run Basic SCF calculations / save results for latter use in mbworkbench

Can accept input from the command line interface (cli) - or - from the YAML input file.
The cli options will overwrite the YAML options.
'''

def get_cli_args():
    '''
    Read CLI arguments for running SCF.

    TODO: add a DFT XC field to SCF section of input file.
    TODO: 'with SOC' field for input file - default is False
    '''

    parser = argparse.ArgumentParser(description='Run an SCF calculation (PySCF is the back end)')
    
    parser.add_argument('--type','-t', 
                        choices=['infile', 'rohf','uhf','ghf','socghf'],
                        default='infile',
                        dest='type',
                        help='type of SCF to run. "ghf" does not include SOC -> use "socghf"!')

    parser.add_argument('--input', '-i', metavar='input file', type=str,
                        action='store',
                        dest='input_file',
                        default='input.yaml',
                        help="name of input file (YAML format)")

    parser.add_argument('--verbose', '-v', action='count',
                        default=4,
                        dest='verb',
                        help='set verbosity level, use more \'v\' chars to make output more verbose (i.e. -vv ) default is -vvvv')

    args = parser.parse_args()
    return args

def scf_core(data, scf_type, scf_params):
    '''
    run SCF within a Workflow context.
    '''
    from pyscf import scf

    try:
        mol = data['mol']
    except KeyError:
        logging.error('[dev] mol not in Workflow data.')
        raise ValueError("Workflow data incomplete")
    
    if scf_type == 'rohf':
        mf = scf.ROHF(mol)
    elif scf_type == 'uhf':
        mf = scf.UHF(mol)
    elif scf_type == 'ghf':
        mf = scf.GHF(mol)
    elif scf_type == 'socghf':
        mf = scf.GHF(mol)
        mf.with_soc = True
    else:
        logging.error(f'[dev] unkown SCF type = {scf_type}')
        raise ValueError
    
    print(f'      [+] PySCF SCF params:')
    for key in scf_params.keys():
        print(f'       {key} : {scf_params[key]}')
        setattr(mf, key, scf_params[key])

    print(f'      [+] Running SCF ... ')
    mf.kernel()



class Scf_Block(blk.Block):
    # TODO: implement SCF block class
    def __init__(self, name):
        super().__init__(name)
        
    def step(self, data):

        scf_type, scf_params = read_scf_input(data['input_file'])

        # make the memomry intensive steps here!
        try:
            assert self.status == blk.OKAY_2_RUN
        except AssertionError:
            logging.error("[dev] Attempted to run block without status OKAY_2_RUN ")
            self.status = blk.FAILED
            self.tear_down()
            return

        try:
            scf_core(data, scf_type, scf_params)
            self.status = blk.COMPLETE
            return
        except Exception as e:
            logging.error(f"SCF block failed to run. Exception: {e}")
            self.status = blk.FAILED
            self.tear_down()
            return

    def check_data(self, data):
        # Q: Are there any general validation checks? or all Block-specific?
        #super().check_data(data)
        try:
            assert 'mol' in data.keys()
            assert isinstance(data['mol'],Mole)
            self.status = blk.OKAY_2_RUN
        except AssertionError as e:
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
            with h5.File(chkfile, 'a') as f:
                generic = 'test'
                if self.name + '/generic' in chkfile:
                    temp = f[self.name + '/generic' ]
                    temp[...] = generic
                else:
                    f.create_dataset(self.name + '/generic', data=generic)
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


def main():
    '''
    Run an standalone SCF calculation based on settings from YAML / cli
    This is meant to be used outside of 'Workflow'. This can be usefule if
    there is difficulty converging an SCF result, for example.
    '''
    from pyscf import scf

    cli = get_cli_args()
    mol  = mol_from_input(cli.input_file)
    if mol is None:
        logging.error('No molecule was recieved, exiting!')
        return

    scf_type, scf_params = read_scf_input(cli.input_file)

    # cli has priority over input file
    if cli.type is not 'infile':
        scf_type = cli.type

    # TODO use input file for verbosity!
    mol.verbose = int(cli.verb)

    # TODO echo settings!
    print("     ====== Run SCF calculation (PySCF) ====== ")
    print(f'      [+] SCF Type : {scf_type}')

    if scf_type == 'rohf':
        mf = scf.ROHF(mol)
    elif scf_type == 'uhf':
        mf = scf.UHF(mol)
    elif scf_type == 'ghf':
        mf = scf.GHF(mol)
    elif scf_type == 'socghf':
        mf = scf.GHF(mol)
        mf.with_soc = True
    else:
        logging.error(f'unkown SCF type = {cli.type}')
    
    print(f'      [+] PySCF SCF params:')
    for key in scf_params.keys():
        print(f'       {key} : {scf_params[key]}')
        setattr(mf, key, scf_params[key])

    print(f'      [+] Running SCF ... ')
    mf.kernel()

if __name__ == '__main__':
    main()
