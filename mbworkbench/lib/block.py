import logging
import h5py as h5

# constant status codes
NOT_RUN = 0
OKAY_2_RUN = 1
RUNNING = 2
COMPLETE = 3
FAILED = -1

class Block:
    '''
    base class for Workflow "Blocks".

    Open Q's / problems:
    1. should each block have an attribute that tells it what
    data it needs as input, and what data it will output to the
    workflow data?

    '''

    def __init__(self, name):
        self.status = NOT_RUN
        # TODO: generate random unique name if name is None
        self.name = name

    def step(self, data, **kwargs):
        self.status = RUNNING
        try:
            print(f"Stepping... {self.name}")
            self.status = COMPLETE
        except:
            self.status = FAILED
            self.tear_down()

    def check_data(self, data):
        # TODO: add general checks on 'data'
        try:
            print("Data Checked...")
            self.status = OKAY_2_RUN
        except:
            self.status = FAILED

    def tear_down(self):
        '''
        defines what to do if a probolem is encountered
        '''
        print("Tore Down...")

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
            logging.error(f"read from checkpoint failed. Excetion {e}")
