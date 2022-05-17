import logging

LOG_LEVELS = {'error' : logging.ERROR,
              'warning' : logging.WARNING,
              'info' : logging.INFO,
              'debug' : logging.DEBUG }

LOG_FORMATS = {'error' : '[%(levelname)s] %(message)s',
              'warning' : '[%(levelname)s] %(message)s',
              'info' : '[%(levelname)s] %(message)s',
              'debug' : '[%(levelname)s] %(created)f : %(funcName)s says : %(message)s' }

#logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(created)f : %(funcName)s says : %(message)s')

# TODO: fix packaging so that we can use 'from mbworkbench.workflow import Workflow', etc.
from mbworkbench.workflow.workflow import Workflow
from mbworkbench.io.cli import get_cli_args


def main():
    '''
    mbworkbench main program.
    '''

    # parse cli args
    cli = get_cli_args()


    logging.debug('Building workflow')
    wkflow = Workflow(cli.input_file)
    wkflow.build()

    # TODO: parse subcommands: run, etc.
    wkflow.run()
    print('Workflow complete.')

if __name__ == '__main__':
    main()
