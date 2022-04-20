import logging

from yaml import safe_load

'''
scf-specific input/output operations
'''

def read_scf_input(fname):
    '''
    read scf parameters from yaml file

    scf parameters should be defined in an scf yaml block:

    scf:
      type: [ rohf | uhf | ghf | socghf ]
      [params:]
        [scf keyword 1 : value 1 ]
        [scf keyword 2 : value 2 ]
        [scf keyword 3 : value 3 ]
      ...

    params is an optional block where PySCF keyword arguments may be supplied.
    where '[scf keyword x : value x ]' indicates that PySCF keywords can optionally
    be included.
    The 'type' field is a required argument, unless the cli is used to set it later.
    '''

    with open(fname,'r') as f:
        infile = safe_load(f) # limits load to simple python objects, just what we want!
    try:
        assert 'scf' in infile.keys()

        scf = infile['scf']
        logging.info(f' found "scf" block in input file ')
        logging.debug(f' input file: scf  = {scf} ')

        if 'type' in scf.keys():
            logging.debug(' found "type" block in "scf" input file block')
            scf_type = scf['type']
        else:
            scf_type = None

        if 'params' in scf.keys():
            logging.debug(' found "params" block in "scf" input file block')
            scf_params = scf['params']
        else:
            scf_params = {}

        return scf_type, scf_params

    except AssertionError:
        logging.info(f' no "scf" block in input file ')
        return None, {}
    
