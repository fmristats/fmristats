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

from numba import jit, njit

import numpy as np

from numpy.linalg import solve, inv

from pandas import DataFrame

import statsmodels.api as sm

import statsmodels.formula.api as smf

from statsmodels.stats.stattools import durbin_watson

########################################################################

def data_at(coordinate, data, epi_code:int, scale:float, radius:float):
    r = radius**2
    s = -2*scale**2
    squared_distances = ((data[...,:3] - coordinate)**2).sum(axis=1)
    valid = np.where(squared_distances < r)
    weights = np.exp(squared_distances[valid] / s)
    data = data[valid]

    df = DataFrame({
        'x'      : data[...,0],
        'y'      : data[...,1],
        'z'      : data[...,2],
        'signal' : data[...,3],
        'time'   : data[...,4],
        'task'   : data[...,5],
        'block'  : data[...,6],
        'cycle'  : data[...,7],
        'slice'  : data[...,8],
        'weight' : weights})

    if abs(epi_code) == 3:
        sortvar = ['time', 'z', 'x', 'y']
    elif abs(epi_code) == 2:
        sortvar = ['time', 'y', 'z', 'x']
    elif abs(epi_code) == 1:
        sortvar = ['time', 'x', 'y', 'z']

    df.sort_values(by=sortvar, inplace=True)

    return df

def model_at(formula, **kwargs):
    data = data_at(**kwargs)
    model = smf.wls(formula, weights=data.weight, data=data)
    return model, data

def fit_at(formula, **kwargs):
    model, data = model_at(formula=formula, **kwargs)
    fit = model.fit()

    data['expected_signal'] = fit.predict()
    data['residual']          = fit.resid
    data['weighted_residual'] = data.weight * fit.resid
    fit.durbin_watson = durbin_watson(data.weighted_residual)

    return fit, model, data

def extract_field(field, param, value, parameter_dict, value_dict):
    return field[..., value_dict[value], parameter_dict[param]]

def fit_field(coordinates, mask, data, design, epi_code:int,
        scale:float, radius:float, verbose=True, backend='numba'):
    """
    Parameters
    ----------
    coordinates : ndarray, shape (…,3)
        The coordinates at which the models shall be fitted. Must be a
        numpy array of which the last dimension must be 3.
    mask : None or ndarray, shape (…)
        A boolean array of the same »layout« as coordinates. The first
        dimensions must match the dimensions of coordinates. If None,
        the model is fitted at all coordinates.
    data : ndarray, shape (…,9), dtype: float
        The array of observations:
            - [...,:3] = coordinates of observation
            - [..., 3] = MR signal response
            - [..., 4] = time of observation
            - [..., 5] = task during time of observation
            - [..., 6] = block number during time of observation
            - [..., 7] = scan cycle
            - [..., 8] = slice number
    design : ndarray, shape (…,p), dtype: float
        The design matrix.
    epi_code : int
    scale : float
    radius : float
    verbose : boo
    """

    ###################################################################
    # Asserts
    ###################################################################

    assert coordinates.shape[:-1] == mask.shape, \
            'shapes of coordinates and mask do not match'
    assert data.shape[:-1] == design.shape[:-1], \
            'shapes of data and design do not match'
    assert epi_code in [-3,-2,-1,1,2,3], 'epi_code must be within [-3,3] but not 0'

    ###################################################################
    # In case you need the Durbin-Watson statistics
    ###################################################################

    if backend == 'statsmodels':
        if abs(epi_code) == 3:
            sortvar = ['time', 'z', 'x', 'y']
        elif abs(epi_code) == 2:
            sortvar = ['time', 'y', 'z', 'x']
        elif abs(epi_code) == 1:
            sortvar = ['time', 'x', 'y', 'z']

    ###################################################################
    # Position of statistics in the statistics field
    ###################################################################

    value_dict = {'point':0, 'stderr':1, 'tstatistic':2,
            'mse':3,
            'df':3, 'degrees_of_freedom':3,
            'dw':3, 'durbin_watson':3}

    if backend == 'statsmodels':
        parameter_dict = {
                'mse':0,
                'df':1, 'degrees_of_freedom':1,
                'dw':2, 'durbin_watson':2}
    else:
        parameter_dict = {
                'mse':0,
                'df':1, 'degrees_of_freedom':1}

    ###################################################################
    # Hyperparameters
    ###################################################################

    r = radius**2
    s = -2*scale**2

    ###################################################################
    # Statistics field
    ###################################################################

    result = np.zeros(coordinates.shape[:-1] + \
                (4,max(design.shape[-1],3),), dtype=float)
    result[...] = np.nan

    ###################################################################
    # Reshape result, coordinates, and mask
    ###################################################################

    rcoordinates = coordinates.reshape((-1,3))
    rresult = result.reshape((-1,4,max(design.shape[-1],3)))

    if mask is None:
        to_fit = np.ones(rcoordinates.shape[0]).astype(bool)
    else:
        to_fit = mask.reshape((-1,))

    # if mask is None:
    #     to_fit = range(rcoordinates.shape[0])
    # else:
    #     to_fit, = np.where(mask.reshape((-1,)))
    #     to_fit  = to_fit.tolist()

    ###################################################################
    # Fit the model
    ###################################################################

    if backend == 'statsmodels':
        fit_sm(rresult, rcoordinates, to_fit, data, design, r, s, sortvar)

    else:
        fit_nb(rresult, rcoordinates, to_fit, data, design, r, s)

    return result, parameter_dict, value_dict

###################################################################
# Backends
###################################################################

def fit_sm(result, coordinates, to_fit, data, design, r, s, sortvar):
    for i in range(coordinates.shape[0]):
        if to_fit[i]:
            squared_distances = ((data[...,:3] - coordinates[i])**2).sum(axis=1)
            valid = np.where(squared_distances < r)
            if valid[0].size > 120:
                weights = np.exp(squared_distances[valid] / s)

                fit = sm.WLS(
                    endog    = data[valid][...,3],
                    exog     = design[valid],
                    weights  = weights,
                    hasconst = True).fit()

                result[i,0]   = fit.params
                result[i,1]   = fit.bse
                result[i,2]   = fit.tvalues
                result[i,3,0] = fit.mse_resid
                result[i,3,1] = fit.df_resid

                df = DataFrame({
                    'x'      : data[valid][...,0],
                    'y'      : data[valid][...,1],
                    'z'      : data[valid][...,2],
                    'time'   : data[valid][...,4]})
                df['weighted_residual'] = weights * fit.resid
                df.sort_values(by=sortvar, inplace=True)
                result[i,3,2] = durbin_watson(df.weighted_residual)


@jit(nopython=True)
def penrose_fit(endog, exog, weights):
    # set up
    w_half        = np.sqrt(weights)
    wendog        = w_half * endog
    wexog         = w_half.reshape((-1,1)) * exog
    # calculation
    pinv_wexog    = np.linalg.pinv(wexog)
    params        = pinv_wexog.dot(wendog)
    fitted_values = exog.dot(params)
    wresid        = wendog - wexog.dot(params)
    df_resid      = wexog.shape[0] - wexog.shape[1]
    mse           = np.dot(wresid, wresid) / df_resid
    # variances, covariances, and standard errors
    cov_params    = mse * np.dot(pinv_wexog, np.transpose(pinv_wexog))
    bse           = np.sqrt(np.diag(cov_params))
    return params, mse, df_resid, bse

@jit(nopython=True, fastmath=True) #, parallel=True)
def fit_nb(result, coordinates, to_fit, data, design, r, s):
    for i in range(coordinates.shape[0]):
        if to_fit[i]:
            squared_distances = ((data[...,:3] - coordinates[i])**2).sum(axis=1)
            valid = np.where(squared_distances < r)
            weights = np.exp(squared_distances[valid] / s)
            endog   = data[valid][...,3]
            exog    = design[valid]
            n = exog.shape[0]
            p = exog.shape[1]
            if (p < n-1) and (np.linalg.matrix_rank(exog) == p):
                params, mse, degrees_of_freedom, bse = penrose_fit(endog, exog, weights)
                result[i,0]   = params
                result[i,1]   = bse
                result[i,2]   = params / bse
                result[i,3,0] = mse
                result[i,3,1] = degrees_of_freedom
