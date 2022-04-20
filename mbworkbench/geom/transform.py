
def shift_geom(atoms, shift):
    '''
    translate geometry (in same format as input) by 'shift' vector
    '''
    trans_atoms = []

    raise NotImplementedError

def _rotate_geom_zprime(atoms, zprime):
    '''
    rotate geometry such that zprime is the new z axis.

    NOTE: the next x-axis, $x' = z \times z'$, and y' = z' \times x'$
    '''

    raise NotImplementedError


