import logging
import h5py as h5

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

    def write_to_chk(self, data, chkfile=None):
        '''
        Write data for Molecule_GTO Block to chkfile.
        '''
        from pyscf.gto import dumps

        try:
            assert chkfile is not None
            with h5.File(chkfile, 'a') as f:
                if '/mol' in f:
                    temp = f['mol']
                    temp[...] = dumps(data['mol'])
                else:
                    f.create_dataset('mol', data=dumps(data['mol']))
        except AssertionError:
            logging.error("no checkpoint file specified")

    def load_from_chk(self, data, chkfile=None):
            '''
            Read Block's data from chkfile, and add it to data

            Dev Note: specific Block types are responsible for serializing/deserializing
            any objects that they have added to 'data'.
            '''
            from pyscf.gto import loads

            try:
                assert chkfile is not None
            except AssertionError:
                logging.error("no checkpoint file specified")

            try:
                with h5.File(chkfile, 'r') as f:
                    if 'mol' in f:
                        data['mol'] = loads(f['mol'][...])
                    else:
                        raise ValueError("No mol in checkpoint file")
            except ValueError as e:
                logging.error(f"read from checkpoint failed. Excetion {e}")
