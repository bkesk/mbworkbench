'''
mbworkbench command line interface.

Whenever possible, getting settings from the input file should be 
prefered to CLI arguments in order to maintain transparency/reproducibility.
Acceptable CLI args include:

1. Input file name
2. possibly some output file names; however, only output files that do not effect the 
    workflow should be editable via the CLI. Otherwise, workflows can easily break.
3. logging configuration

'''

import argparse

def get_cli_args():
    '''
    Get high-level cli

    #TODO: add subcommands for specific tasks
    '''

    parser = argparse.ArgumentParser(description='Run a many-body workflow')
    
    parser.add_argument('--input', '-i', metavar='input file', type=str,
                        action='store',
                        dest='input_file',
                        default='input.yaml',
                        help="name of input file (YAML format)")

    parser.add_argument('--log-level','-ll', 
                        choices=['error', 'warning','info','debug'],
                        default='warning',
                        dest='log_level',
                        help='logging level for Python. order of increasing verbosity : error -> warning -> info -> debug')

    args = parser.parse_args()
    return args
