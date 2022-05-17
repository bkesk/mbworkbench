'''
AFQMC interfaces

here, we will import only the interfaces to AFQMC codes that have
a Python interface currently installed.
'''
import logging
from os import environ

AFQMCLAB_DIR = environ.get('AFQMCLAB_DIR')

if AFQMCLAB_DIR is not None:
    logging.debug("found AFQMCLAB_LAB in environment")
    from afqmclab import *
else:
    logging.debug("did not find AFQMCLAB_DIR in environment")

from .afqmclab import Write_Afqmclab
