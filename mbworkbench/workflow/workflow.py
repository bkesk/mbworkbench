import logging
import h5py as h5

from mbworkbench.lib import block as blk
from mbworkbench.io.input_file import read_workflow_def

class Workflow:
    '''
    Workflow class manages a 'workflow' which consists of
    several 'Blocks' which each perform some step in the
    workflow.
    '''

    def __init__(self, infile=None):

        try:
            assert infile is not None
        except AssertionError:
            logging.error("Workflow was not given an input file")

        self.data = { 'input_file' : infile }
        self.blocks = [] # defines each block
        self.current_block = 0  # points to current block


    def run(self):
        for block in self.blocks[self.current_block:]:

            logging.info("Checking if block can run")
            block.check_data(self.data)

            if block.status == blk.OKAY_2_RUN:
                logging.info("running block")
                block.step(self.data)
                self.current_block += 1

            if block.status == blk.FAILED:
                logging.error("Block failed, pausing workflow")
                self.pause()
                break

        self.write_chk()

    def build(self, definition=None):

        if definition is None:
            definition = read_workflow_def(self.data['input_file'])

        self.blocks = generate_blocks(definition)


    def restart(self):
        self.load_chk()
        self.run()


    def pause(self):
        self.write_chk()


    def write_chk(self):
        '''
        Write's workflow info to checkpoint, and calls the current block's 'write_to_chk' method.

        Each Block class is responsible for reading/writing it's own data, and for serializing/deserializing data.
        '''
        logging.info('Writing workflow.chk')
        def _try_write(f, key, data):
            try:
                f.create_dataset(key, data=data)
            except ValueError:
                current_data = f[key]
                current_data[...] = data
            except TypeError as e:
                logging.warning(f"Couldn't write {key} to checkpoint. Exception: {e} ")

        with h5.File('workflow.chk','a') as f:            
            #for key in self.data:
            #    _try_write(f, 'Data/' + key, self.data[key])
            _try_write(f, 'current_block', self.current_block)
            
        for block in self.blocks[:self.current_block]:
            block.write_to_chk(self.data, 'workflow.chk')

    def load_chk(self):
        '''
        Load's Workflow data from checkpoint, and all block data
        *prior* to the current block.
        '''
        logging.info('Loading workflow.chk')
        try:
            with h5.File('workflow.chk','r') as f:
                self.current_block = f['current_block'][...]
                
            for block in self.blocks[:self.current_block]:
                block.load_from_chk(self.data, 'workflow.chk')

        except Exception as e:
            logging.warning(f"read from checkpoint workflow.chk failed: {e}")

############################################
#                                          #
#       Define known Block classes         #  
#                                          #
############################################

from mbworkbench.scf.hf import Scf_Block
from mbworkbench.system.molecule import Molecule_GTO

KNOWN_BLOCKS = {'base' : blk.Block,
                'molecule' : Molecule_GTO,
                'scf' : Scf_Block}

def generate_blocks(definition):
    '''
    generate blocks based on the givien definition.
    
    notes: we may want to try to look these up as we are trying to setup/run them,
    not all at once! Look at this on the dry erase board first!
    '''
    
    blocks = []

    for i,block in enumerate(definition):
        logging.info(f'Adding Block {i} : {block}')
        try:
            blocks.append(KNOWN_BLOCKS[block](name=block))
        except KeyError:
            logging.error("Invalid block name in workflow")
        except Exception as e:
            logging.error(f"Incountered Error {e}")

    return blocks

