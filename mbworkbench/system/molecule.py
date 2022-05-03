import logging

from mbworkbench.lib import block as blk
from mbworkbench.io.input_file import mol_from_input

class Molecule_GTO(blk.Block):
    '''
    Block subclass for building a molecular system in a GTO basis.
    '''

    def __init__(self, name):
        super().__init__(name)

    def step(self, data, **kwargs):
        try:
            data['mol'] = mol_from_input(fname=data['input_file'])
            self.status = blk.COMPLETE
        except:
            logging.error('unable to build molecule')
            self.status = blk.FAILED
