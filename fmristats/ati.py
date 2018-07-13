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

Fit signal model to a 3D image

"""

import statsmodels.api as sm

import statsmodels.formula.api as smf

import numpy as np

########################################################################

def extract_field(field, param, value, parameter_dict, value_dict):
    return field[..., value_dict[value], parameter_dict[param]]

def fit_field(coordinates, mask, endog, exog, agc,
        scale:float, radius:float, verbose=True, name=None):

    assert agc.shape[0] == endog.shape[0], 'shapes do not match'
    assert agc.shape[0] == exog.shape[0], 'shapes do not match'

    double_squared_scale = -2*scale**2

    p = exog.shape[1]
    result = np.zeros(coordinates.shape[:-1] + (4,max(p,3),), dtype=float)
    result[...,0,:3] = coordinates.copy()

    if mask is None:
        points_to_estimate = result.reshape((-1,) + result.shape[-2:])
        partial_fit = False

        if verbose:
            print("""{}:
            …coordinates in the domain to estimate: {:>10,d}""".format(
                name, points_to_estimate.shape[0]))
    else:
        assert type(mask) is np.ndarray, 'mask must be an ndarray'
        assert mask.dtype == bool, \
                'mask is not None it must be of dtype bool'
        assert mask.shape == coordinates.shape[:-1], \
                'image dimensions of coordinates and mask must match'
        result[~mask] = np.nan
        points_to_estimate = result[mask]
        partial_fit = True

        if verbose:
            print("""{}:
            …coordinates in the domain to estimate: {:>10,d}
            …coordinates in the domain not to:      {:>10,d}""".format(
                name, points_to_estimate.shape[0], (~mask).sum()))

    for x in iter(points_to_estimate):
        distances = ((agc - x[...,0,:3])**2).sum(axis=1)
        valid = np.less(distances, radius**2)
        if valid.sum() > 120:
            try:
                weights = np.exp(distances[valid] / double_squared_scale)
                fit = sm.WLS(
                    endog    = endog[valid].copy(),
                    exog     = exog[valid].copy(),
                    weights  = weights,
                    hasconst = True).fit()

                x[0]   = fit.params
                x[1]   = fit.bse
                x[2]   = fit.tvalues
                x[3,0] = fit.mse_resid
                x[3,1] = fit.df_resid
            except Exception as e:
                print("""{}: Exception at {}: {}""".format(
                    name, x[...,0,:3], e))
                x[...] = np.nan
        else:
            x[...] = np.nan

    if partial_fit:
        result[mask] = points_to_estimate

    # first position
    value_dict = {'point':0, 'stderr':1, 'tstatistic':2,
            'mse':3,
            'df':3, 'degrees_of_freedom':3}

    # second position
    parameter_dict = {
            'mse':0,
            'df':1, 'degrees_of_freedom':1}

    return result, parameter_dict, value_dict
