# constant status codes
NOT_RUN = 0
OKAY_2_RUN = 1
RUNNING = 2
COMPLETE = 3
FAILED = -1

class Block:
    '''
    base class for Workflow "Blocks".
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

