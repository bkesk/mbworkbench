import logging

from mbworkbench.workflow.workflow import Workflow
from mbworkbench.io.cli import get_cli_args


LOG_LEVELS = {'error' : logging.ERROR,
              'warning' : logging.WARNING,
              'info' : logging.INFO,
              'debug' : logging.DEBUG }

LOG_FORMATS = {'error' : '[%(levelname)s] %(message)s',
              'warning' : '[%(levelname)s] %(message)s',
              'info' : '[%(levelname)s] %(message)s',
              'debug' : '[%(levelname)s] %(created)f : %(funcName)s says : %(message)s' }

def main():
    '''
    mbworkbench main program.
    '''

    # parse cli args
    cli = get_cli_args()

    

    logging.basicConfig(level=LOG_LEVELS[cli.log_level], 
                        format=LOG_FORMATS[cli.log_level])

    logging.debug('Building workflow')
    wkflow = Workflow(cli.input_file)
    wkflow.build()

    # TODO: parse subcommands: run, etc.
    wkflow.run()

if __name__ == '__main__':
    main()
