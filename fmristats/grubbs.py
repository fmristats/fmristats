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

Methods for outlier detection.

"""

import numpy as np

from scipy.stats.distributions import t

def grubbs_test(xvec, sgnf, null='equal'):
    """
    The Grubbs' test for outlyer detection

    Parameters
    ----------
    xvec : ndarray, shape (n,)
    sgnf : float
    null : str

    Notes
    -----
    Specification of the null hypothesis can either be 'equal'
    (default), 'greater', or 'less'. Only 'equal' is implemented.

    Missigs (nan values) are ignored in the test.
    """
    assert np.not_equal(np.nanvar(xvec), 0), 'variance of observations (xvec) is equal to 0'

    n = sum(~np.isnan(xvec))
    m = np.nanmean(xvec)
    s = np.sqrt(np.nanvar(xvec))
    candidates = abs(xvec - m)
    gstat = np.nanmax(candidates) / s
    tstat = t.ppf(q=1-(sgnf / (2*n)), df=n-2)**2
    crval = ((n-1)/np.sqrt(n)) * np.sqrt(tstat / (n-2+tstat))
    if (gstat > crval):
        return np.nanargmax(candidates)
    else:
        return None

def grubbs(xvec, sgnf, inplace=False):
    """
    Iterate Grubbs' test

    Parameters
    ----------
    xvec : ndarray, shape (n,)
    sgnf : float
    inplace : bool

    Notes
    -----
    Missigs (nan values) are set to True.
    """
    if not inplace:
        xvec = xvec.copy()
    args = np.repeat(False, len(xvec))
    args [ np.isnan(xvec) ] = True
    while (sum(~np.isnan(xvec)) > 6 & np.not_equal(np.nanvar(xvec), 0)):
        arg = grubbs_test(xvec, sgnf)
        if arg == None:
            break
        xvec [arg] = np.nan
        args [arg] = True
    return xvec, args
