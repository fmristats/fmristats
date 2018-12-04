# Copyright 2016-2018 Thomas W. D. Möbius
#
# This file is part of fmristats.
#
# fmristats is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# fmristats is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# It is not allowed to remove this copy right statement.

"""

Will parse a MATLAB coded stimulus design and create a block design

From: 8.3 Factorial Design - `Multiple Conditions`_ (Page 65):

    »If you have multiple conditions then entering the details a
    condition at a time is very inefficient. This option can be used to
    load all the required information in one go.

    You will need to create a *.mat file containing the relevant
    information. This *.mat file must include the following cell arrays:
    names, onsets and durations eg. names{2}=’SSent-DSpeak’,
    onsets{2}=[3 5 19 222], durations{2}=[0 0 0 0] contain the required
    details of the second condition. These cell arrays may be made
    available by your stimulus delivery program eg. COGENT. The duration
    vectors can contain a single entry if the durations are identical
    for all events.«

.. _`Multiple Conditions`: http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf

"""

########################################################################
#
# Tool
#
########################################################################

from .name import Identifier

from .stimulus import Block

def mat2block(mat, name):
    """
    Will parse mat and create a block stimulus.

    Parameters
    ----------
    mat : dict
        dict returned by scipy.io.loadmat(args.input)
    name : Identifier
        Unique identifier of the subject.
    """
    assert type(name) is Identifier, 'name must be of type Identifier'
    block_number = mat['names'].shape[-1]
    names = [mat['names'].ravel()[i][0] for i in range(block_number)]
    onsets = { names[i] : mat['onsets'].ravel()[i].ravel() for i in range(block_number)}
    durations = { names[i] : mat['durations'].ravel()[i].ravel()[0] for i in range(block_number)}
    return Block(name=name, names=names, onsets=onsets, durations=durations)
