import logging
import argparse

from mbworkbench.io.input_file import mol_from_input

'''
run Basic SCF calculations / save results for latter use in mbworkbench
'''


def get_cli_args():
    '''
    Read CLI arguments for running SCF.

    TODO: consider including an 'SCF' section to the input file.
    TODO: add a DFT XC field to SCF section of input file.
    TODO: 'with SOC' field for input file - default is False
    '''

    parser = argparse.ArgumentParser(description='Run an SCF calculation (PySCF is the back end)')
    
    parser.add_argument('type', 
                        choices=['rohf','uhf','ghf','socghf'],
                        default='rohf',
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
        logging.error('No molecule was generated, exiting')
        return 10

    mol.verbose = int(cli.verb)

    # TODO echo settings!

    # TODO: try to read SCF type from input file - keep cli. CLI should be prioritized.

    if cli.type == 'rohf':
        mf = scf.ROHF(mol)
    elif cli.type == 'uhf':
        mf = scf.UHF(mol)
    elif cli.type == 'ghf':
        mf = scf.GHF(mol)
    elif cli.type == 'socghf':
        mf = scf.GHF(mol)
        mf.with_soc = True
    else:
        logging.error(f'unkown SCF type = {cli.type}')
        
    mf.kernel()

if __name__ == '__main__':
    main()
