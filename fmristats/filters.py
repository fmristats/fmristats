"""

This module contains methods for foreground detection.

"""

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

import numpy as np

from skimage import filters

from pandas import DataFrame

def fit_foreground(data, ep):
    """
    Fit foreground

    Parameters
    ----------
    data : ndarray, shape (n,a,b,c)
    ep : int

    Returns
    -------
    DataFrame :
        Summary of thresholds, mean signals, and coefficients of
        variation.

    Notes
    -----
    This will overwrite your data and set all estimated background
    pixels to nan.
    """
    if np.isnan(data).any():
        print('  …data contain {} nan values, set to 0'.format(np.isnan(data).sum()))
        data[np.isnan(data)] = 0

    dat = np.moveaxis(data, ep+1, 1)
    cn = dat.shape[0]
    pn = dat.shape[1]

    ct = np.zeros(cn)
    pt = np.zeros((pn,cn))

    for c in range(cn):
        ct[c] = filters.threshold_otsu(dat[c])
        for p in range(pn):
            pt[p,c] = filters.threshold_otsu(dat[c,p])

    dat = np.moveaxis(dat, (0,1), (-1,-2))
    msk = (dat < ct) | (dat < pt)
    dat[msk] = np.nan
    return ct, pt
