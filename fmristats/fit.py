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

import statsmodels.api as sm

import statsmodels.formula.api as smf

from statsmodels.stats.stattools import durbin_watson

########################################################################

# TODO:
def data_at(coordinate, data, scale:float, radius:float):
    agc = data[['i','j','k']].values
    distances = ((agc - coordinate)**2).sum(axis=1)
    valid = np.less(distances, radius**2)

    double_squared_scale = -2*(scale**2)
    weights  = np.exp( distances[valid] / double_squared_scale )

    data = data[valid].copy()
    data['weight'] = weights
    return data

# TODO:
def model_at(formula, timevec=False, **kwargs):
    data = data_at(**kwargs)
    m = smf.wls(formula, weights=data.weight, data=data)

    if timevec:
        m.timevec = data['time']

    return m

def extract_field(field, param, value, parameter_dict, value_dict):
    return field[..., value_dict[value], parameter_dict[param]]

def fit_field(coordinates, mask, data, design, ep:int, scale:float,
        radius:float, verbose=True,
        durbin_watson=False, backend='jit'):
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
    ep : int
    scale : float
    radius : float
    verbose : boo
    """

    ###################################################################
    # Asssets
    ###################################################################

    assert coordinates.shape[:-1] == mask.shape, \
            'shapes of coordinates and mask do not match'
    assert data.shape[:-1] == design.shape[:-1], \
            'shapes of data and design do not match'
    assert ep in [0,1,2], 'ep must be one of 0, 1, or 2'

    ###################################################################
    # In case you need the Durbin-Watson statistics
    ###################################################################

    if durbin_watson:
        if ep == 2:
            sortvar = ['time', 'k', 'i', 'j']
        elif ep == 1:
            sortvar = ['time', 'j', 'k', 'i']
        elif ep == 0:
            sortvar = ['time', 'i', 'j', 'k']

    ###################################################################
    # Position of statistics in the statistics field
    ###################################################################

    value_dict = {'point':0, 'stderr':1, 'tstatistic':2,
            'mse':3,
            'df':3, 'degrees_of_freedom':3,
            'dw':3, 'durbin_watson':3}

    if durbin_watson:
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
    # Fit the model
    ###################################################################

    if backend == 'statsmodels':
        if durbin_watson:
            if mask is None:
                fit_sm_nm_dw(result, coordinates, data, design, r, s, sortvar)
            else:
                fit_sm_wm_dw(result, coordinates, mask, data, design, r, s, sortvar)
        else:
            if mask is None:
                fit_sm_nm(result, coordinates, data, design, r, s)
            else:
                fit_sm_wm(result, coordinates, mask, data, design, r, s)

    else:
        if mask is None:
            print('JIT with no mask')
            fit_nm(result, coordinates, data, design, r, s)
        else:
            print('JIT with mask')
            fit_wm(result, coordinates, mask, data, design, r, s)

    return result, parameter_dict, value_dict

###################################################################
# Backends
###################################################################

def fit_sm_nm_dw(result, coordinates, data, design, r, s, sortvar):
    ni, nj, nk, _ = coordinates.shape
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                squared_distances = ((data[...,:3] - coordinates[i,j,k])**2).sum(axis=1)
                valid = np.where(squared_distances < r)
                if valid[0].size > 120:
                    weights = np.exp(squared_distances[valid] / s)
                    fit = sm.WLS(
                        endog    = data[valid][...,3],
                        exog     = design[valid],
                        weights  = weights,
                        hasconst = True).fit()
                    result[i,j,k,0]   = fit.params
                    result[i,j,k,1]   = fit.bse
                    result[i,j,k,2]   = fit.tvalues
                    result[i,j,k,3,0] = fit.mse_resid
                    result[i,j,k,3,1] = fit.df_resid
                    #df = dataframe[valid].copy()
                    #df.reweighted_residual = weights * fit.resid
                    #df.sort_values(by=sortvar, inplace=True)
                    #x[3,2] = durbin_watson(df.reweighted_residual)

def fit_sm_wm_dw(result, coordinates, mask, data, design, r, s, sortvar):
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
                            exog     = design[valid],
                            weights  = weights,
                            hasconst = True).fit()
                        result[i,j,k,0]   = fit.params
                        result[i,j,k,1]   = fit.bse
                        result[i,j,k,2]   = fit.tvalues
                        result[i,j,k,3,0] = fit.mse_resid
                        result[i,j,k,3,1] = fit.df_resid
                        #df = dataframe[valid].copy()
                        #df.reweighted_residual = weights * fit.resid
                        #df.sort_values(by=sortvar, inplace=True)
                        #x[3,2] = durbin_watson(df.reweighted_residual)

def fit_sm_nm(result, coordinates, data, design, r, s):
    ni, nj, nk, _ = coordinates.shape
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                squared_distances = ((data[...,:3] - coordinates[i,j,k])**2).sum(axis=1)
                valid = np.where(squared_distances < r)
                if valid[0].size > 120:
                    weights = np.exp(squared_distances[valid] / s)
                    fit = sm.WLS(
                        endog    = data[valid][...,3],
                        exog     = design[valid],
                        weights  = weights,
                        hasconst = True).fit()
                    result[i,j,k,0]   = fit.params
                    result[i,j,k,1]   = fit.bse
                    result[i,j,k,2]   = fit.tvalues
                    result[i,j,k,3,0] = fit.mse_resid
                    result[i,j,k,3,1] = fit.df_resid

def fit_sm_wm(result, coordinates, mask, data, design, r, s):
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
                            exog     = design[valid],
                            weights  = weights,
                            hasconst = True).fit()
                        result[i,j,k,0]   = fit.params
                        result[i,j,k,1]   = fit.bse
                        result[i,j,k,2]   = fit.tvalues
                        result[i,j,k,3,0] = fit.mse_resid
                        result[i,j,k,3,1] = fit.df_resid

@jit(nopython=True, fastmath=True, parallel=True)
def fit_nm(result, coordinates, data, design, r, s):
    ni, nj, nk, _ = coordinates.shape
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                squared_distances = ((data[...,:3] - coordinates[i,j,k])**2).sum(axis=1)
                valid = np.where(squared_distances < r)
                if valid[0].size > 120:
                    weights = np.exp(squared_distances[valid] / s)
                    endog   = data[valid][...,3]
                    exog    = design[valid]

                    W = np.diag(weights)
                    V_inverse = exog.T.dot(W).dot(exog)
                    params    = solve(V_inverse, exog.T.dot(W).dot(endog))   # regression parameters
                    residuals = endog - exog.dot(params)                     # residuals
                    df_resid  = exog.shape[0] - exog.shape[1]                # degrees of freedom
                    mse_resid = residuals.T.dot(W).dot(residuals) / df_resid # mean squared error
                    bse       = np.sqrt(mse_resid*np.diag(inv(V_inverse)))   # standard error

                    result[i,j,k,0]   = params
                    result[i,j,k,1]   = bse
                    result[i,j,k,2]   = params / bse
                    result[i,j,k,3,0] = mse_resid
                    result[i,j,k,3,1] = df_resid

@jit(nopython=True, fastmath=True, parallel=True)
def fit_wm(result, coordinates, mask, data, design, r, s):
    ni, nj, nk, _ = coordinates.shape
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                if mask[i,j,k]:
                    squared_distances = ((data[...,:3] - coordinates[i,j,k])**2).sum(axis=1)
                    valid = np.where(squared_distances < r)
                    if valid[0].size > 120:
                        weights = np.exp(squared_distances[valid] / s)
                        endog   = data[valid][...,3]
                        exog    = design[valid]

                        W = np.diag(weights)
                        V_inverse = exog.T.dot(W).dot(exog)
                        params    = solve(V_inverse, exog.T.dot(W).dot(endog))   # regression parameters
                        residuals = endog - exog.dot(params)                     # residuals
                        df_resid  = exog.shape[0] - exog.shape[1]                # degrees of freedom
                        mse_resid = residuals.T.dot(W).dot(residuals) / df_resid # mean squared error
                        bse       = np.sqrt(mse_resid*np.diag(inv(V_inverse)))   # standard error

                        result[i,j,k,0]   = params
                        result[i,j,k,1]   = bse
                        result[i,j,k,2]   = params / bse
                        result[i,j,k,3,0] = mse_resid
                        result[i,j,k,3,1] = df_resid
