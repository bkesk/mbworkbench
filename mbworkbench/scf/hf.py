import logging
import argparse

from mbworkbench.io.input_file import mol_from_input
from mbworkbench.scf.io import read_scf_input

'''
run Basic SCF calculations / save results for latter use in mbworkbench

Can accept input from the command line interface (cli) - or - from the YAML input file.
The cli options will overwrite the YAML options.
'''

def get_cli_args():
    '''
    Read CLI arguments for running SCF.

    TODO: consider including an 'SCF' section to the input file.
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

def main():
    '''
    Run an SCF calculation based on settings from YAML / cli
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
