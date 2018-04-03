# Copyright 2016-2017 Thomas W. D. Möbius
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

Handling observations

"""

import numpy as np

def create_slice_time_matrix(shape, epi_code, temporal_resolution,
        interleaved=None):
    """
    Creates a time vector from parameters

    Parameters
    ----------
    shape : tuple, type: int
    epi_code : int or None
    temporal_resolution : float
        The temporal resolution.
    interleaved : bool

    Returns
    -------
    slice_time : ndarray of float, shape (M,) or (M,N,O,P)
        The time vector that holds the slices time information. The
        vector holds the start point of the time at which a particular
        slice has been measured.  If ep is None, slice_time has
        shape (M,) with M equal to the number of full scan cycles.  If
        ep is int, then slice_time has
        the same shape as the data.

    """
    # TODO: implement interleaved
    assert not interleaved, 'sorry, interleaved not implemented yet'
    n = shape[0]
    m = shape[abs(epi_code)]
    slice_time = np.arange(
            start = 0,
            stop  = n*temporal_resolution,
            step  = temporal_resolution/m).reshape((n,m))
    return slice_time
