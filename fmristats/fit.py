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

Fit signal model to signal data

"""

import statsmodels.api as sm

import statsmodels.formula.api as smf

from statsmodels.stats.stattools import durbin_watson

import numpy as np

########################################################################

def extract_field(field, param, value, parameter_dict, value_dict):
    return field[..., value_dict[value], parameter_dict[param]]

def fit_field(coordinates, mask, data, exog,
        ep:int, scale:float, radius:float, verbose=True,
        name:str='undefined'):
    """
    Parameters
    ----------
    coordinates : ndarray, shape (…,3)
        The coordinates at which the models shall be fitted. Must be a
        numpy array of which the last dimension must be 3.
    mask : ndarray, shape (…)
        A boolean array of the same »layout« as coordinates. The first
        dimensions must match the dimensions of coordinates.
    data : ndarray, shape (…,9), dtype: float
        The array of observations:
            - [...,:3] = coordinates of observation
            - [..., 3] = MR signal response
            - [..., 4] = time of observation
            - [..., 5] = task during time of observation
            - [..., 6] = block number during time of observation
            - [..., 7] = scan cycle
            - [..., 8] = slice number
    ep : int
    scale : float
    radius : float
    verbose : boo
    name : str
    """
    assert ep in [0,1,2], 'ep must be one of 0, 1, or 2'

    r = radius**2
    s = -2*scale**2

    result = np.zeros(coordinates.shape[:-1] + \
                (4,max(exog.shape[-1],3),), dtype=float)
    result[...] = np.nan

    %%timeit
    ni, nj, nk, _ = coordinates.shape
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                if mask[i,j,k]:
                    squared_distances = ((data[...,:3] - coordinates[i,j,k])**2).sum(axis=1)
                    valid = np.where(squared_distances < r)
                    if valid[0].size > 120:
                        weights = np.exp(squared_distances[valid] / s)
                        fit = sm.WLS(
                            endog    = data[valid][...,3],
                            exog     = exog[valid],
                            weights  = weights,
                            hasconst = True).fit()
                        result[i,j,k,0]   = fit.params
                        result[i,j,k,1]   = fit.bse
                        result[i,j,k,2]   = fit.tvalues
                        result[i,j,k,3,0] = fit.mse_resid
                        result[i,j,k,3,1] = fit.df_resid








    assert agc.shape[0] == endog.shape[0], 'shapes do not match'
    assert agc.shape[0] == exog.shape[0], 'shapes do not match'

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

    if ep == 2:
        sortvar = ['time', 'k', 'i', 'j']
    elif ep == 1:
        sortvar = ['time', 'j', 'k', 'i']
    elif ep == 0:
        sortvar = ['time', 'i', 'j', 'k']

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

                df = dataframe[valid].copy()
                df.reweighted_residual = weights * fit.resid
                df.sort_values(by=sortvar, inplace=True)

                x[0]   = fit.params
                x[1]   = fit.bse
                x[2]   = fit.tvalues
                x[3,0] = fit.mse_resid
                x[3,1] = fit.df_resid
                x[3,2] = durbin_watson(df.reweighted_residual)
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
            'df':3, 'degrees_of_freedom':3,
            'dw':3, 'durbin_watson':3}

    # second position
    parameter_dict = {
            'mse':0,
            'df':1, 'degrees_of_freedom':1,
            'dw':2, 'durbin_watson':2}

    return result, parameter_dict, value_dict






def fit_field(coordinates, mask, endog, exog, agc, dataframe, ep,
        scale:float, radius:float, verbose=True, name=None):
    """
    Parameters
    ----------
    coordinates : ndarray, shape (…,3)
        The coordinates at which the models shall be fitted. Must be a
        numpy array of which the last dimension must be 3.
    mask : ndarray, shape (…)
        A boolean array of the same »layout« as coordinates. The first
        dimensions must match the dimensions of coordinates.
    endog : ndarray, shape (n,)
        The vector of observations
    """





    assert ep in [0,1,2], 'ep must be one of 0, 1, or 2'
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

    if ep == 2:
        sortvar = ['time', 'k', 'i', 'j']
    elif ep == 1:
        sortvar = ['time', 'j', 'k', 'i']
    elif ep == 0:
        sortvar = ['time', 'i', 'j', 'k']

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

                df = dataframe[valid].copy()
                df.reweighted_residual = weights * fit.resid
                df.sort_values(by=sortvar, inplace=True)

                x[0]   = fit.params
                x[1]   = fit.bse
                x[2]   = fit.tvalues
                x[3,0] = fit.mse_resid
                x[3,1] = fit.df_resid
                x[3,2] = durbin_watson(df.reweighted_residual)
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
            'df':3, 'degrees_of_freedom':3,
            'dw':3, 'durbin_watson':3}

    # second position
    parameter_dict = {
            'mse':0,
            'df':1, 'degrees_of_freedom':1,
            'dw':2, 'durbin_watson':2}

    return result, parameter_dict, value_dict

def data_at(coordinate, data, scale : float, radius : float):
    agc = data[['i','j','k']].values
    distances = ((agc - coordinate)**2).sum(axis=1)
    valid = np.less(distances, radius**2)

    double_squared_scale = -2*(scale**2)
    weights  = np.exp( distances[valid] / double_squared_scale )

    data = data[valid].copy()
    data['weight'] = weights
    return data

def model_at(formula, timevec=False, **kwargs):
    data = data_at(**kwargs)
    m = smf.wls(formula, weights=data.weight, data=data)

    if timevec:
        m.timevec = data['time']

    return m
