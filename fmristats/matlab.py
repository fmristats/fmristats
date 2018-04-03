"""
Will parse a the output of ... to a irritation.
"""

########################################################################
#
# Tool
#
########################################################################

from .name import Identifier

from .irritation import Block

def mat2irritation(mat, name):
    """
    Will parse a the output of ... to a block irritation.

    Parameters
    ----------
    mat : dict
        dict returned by scipy.io.loadmat(args.input)
    name : Identifier
        identifier of the sesson.
    """
    assert type(name) is Identifier, 'name must be of type Identifier'
    block_number = mat['names'].shape[-1]
    names = [mat['names'].ravel()[i][0] for i in range(block_number)]
    onsets = { names[i] : mat['onsets'].ravel()[i].ravel() for i in range(block_number)}
    durations = { names[i] : mat['durations'].ravel()[i].ravel()[0] for i in range(block_number)}
    return Block(name=name, names=names, onsets=onsets, durations=durations)
