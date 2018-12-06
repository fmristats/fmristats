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

"""

Calculates the population effect field from fMRI data of the brain.

"""

import statsmodels.api as sm

import numpy as np

from numpy.linalg import solve, svd

def _kernel_matrix(h, d, x):
    _, s, v = svd(x.T)
    K = v[len(s):].T
    return K.dot(solve(K.T.dot(np.diag(h+d)).dot(K), K.T))

def _qprofile_analysis(h, y, d):
    """
    In case of no covariates, i.e., with only an intercept in the
    regression, we are in the case of an meta analysis. The q-profile
    funtion q_τ simplifies:

    .. math:

    q_δ(τ) = \sum_j \frac 1{δ_i + τ} ⋅ (y - \bar y)^2

    """
    return ((y - y.mean())**2 / (h+d)).sum()

def _qprofile_regression(h, y, d, x):
    return y.dot( _kernel_matrix(h=h,d=d,x=x) ).dot(y)

def hedge_estimator(y, d):
    """
    Hedge estimate for the heterogeneity

    Parameters
    ----------
    y : ndarray, shape (k,)
        The observations.
    d : ndarray, shape (k,)
        The heteroscedasticity vector.
    x : ndarray, shape (k,p)
        The covariates.

    Return
    ------
    λ : float
        The estimated heterogeneity.
    """
    n = len (y)
    return float( max(0, (((y - y.mean())**2).sum() - (d -  d/n).sum()) / (n-1)))

def hedge_type_estimator(y, d, x):
    """
    Hedge estimate for the heterogeneity

    Parameters
    ----------
    y : ndarray, shape (k,)
        The observations.
    d : ndarray, shape (k,)
        The heteroscedasticity vector.
    x : ndarray, shape (k,p)
        The covariates.

    Return
    ------
    λ : float
        The estimated heterogeneity.
    """
    # TODO: do this with MP-pseudo inverse instead
    H  = x.dot(solve(x.T.dot(x), x.T))
    E  = np.eye(x.shape[0]) - H
    resid_df = -np.diff(x.shape)
    return float( max(0, (y.dot(E).dot(y) - np.trace(E.dot(np.diag(d)))) / resid_df) )

def meta_analysis(y, d):
    """
    Random effect meta regression

    Will regress y onto x using a random effects meta regression model
    with heteroscedasticity d.  Heterogeneity will be estimated by Hedge.

    Parameters
    ----------
    y : ndarray, shape (k,)
        The observations.
    d : ndarray, shape (k,)
        The heteroscedasticity vector.

    Return
    ------
    β : ndarray, shape (p)
        The location parameters.
    t : ndarray, shape (p)
        Knapp-Hartung adjusted t-test statistic
    λ : float
        The heterogeneity.
    adj_stderr : float
        The adjusted standard error of the location estimator
    rdf : int
        Residual degrees of freedom

    Notes
    -----
    In case of a meta analysis, many formula simplify.

    .. math:

        σ(\hat μ) = \frac 1 {\sum_j (1/ δ_j + τ)}
        \hat μ    = \frac {\sum_j (1/ δ_j + τ) * y_j}{\sum_j (1/ δ_j + τ)}
    """
    h   = hedge_estimator(y,d)
    w   = 1/(h+d)
    v   = 1 / w.sum()
    b   = (w * y).sum() / w.sum()
    rdf = len(y) - 1
    adj = _qprofile_analysis(h=h,y=y,d=d) / rdf
    adj_stderr = np.sqrt(v * adj)
    t = b / adj_stderr
    return b, t, h, adj_stderr, rdf

def meta_regression(y, d, x):
    """
    Random effect meta regression

    Will regress y onto x using a random effects meta regression model
    with heteroscedasticity d.  Heterogeneity will be estimated by Hedge.

    Parameters
    ----------
    y : ndarray, shape (k,)
        The observations.
    d : ndarray, shape (k,)
        The heteroscedasticity vector.
    x : ndarray, shape (k,p)
        The covariates.

    Return
    ------
    β : ndarray, shape (p)
        The location parameters.
    t : ndarray, shape (p)
        Knapp-Hartung adjusted t-test statistic
    λ : float
        The heterogeneity.
    adj_stderr : float
        The adjusted standard error of the location estimator
    rdf : int
        Residual degrees of freedom
    """
    h   = hedge_type_estimator(y,d,x)
    fit = sm.WLS(y, x, weights=1/(h+d)).fit()
    v   = fit.bse
    b   = fit.params
    rdf = fit.df_resid
    adj = _qprofile_regression(h=h,y=y,d=d,x=x) / rdf
    adj_stderr = v * np.sqrt(adj)
    t   = b / adj_stderr
    return b, t, h, adj_stderr, rdf

def fit_field(statistics, design=None, mask=None):
    """
    Random effect meta regression

    Will regress each point of an effect and certainty field onto the
    covariates in x using a random effects meta regression model

    In the following, k is the number of subjects and p the number of
    regression parameters in the meta regression model.

    Parameters
    ----------
    statistics : ndarray, shape (x,y,z,3,k)
        The observations. In 0: BOLD effect, 1: BOLD stderr, 2: NA
    mask : ndarray, shape 3D-image
        A mask where to fit the field
    design : ndarray, shape (k,p)
        The design matrix. If None, a meta analysis will be performed.

    Returns
    -----
    ndarray, shape (x,y,z,3,p+1)
    """

    ###################################################################
    # Parameter names in the fit
    ###################################################################

    if design is None:
        p  = 1
    else:
        p  = design.shape[1]

    assert statistics.shape[-1] > p+1, 'model not identifiable, too many parameters to fit'

    ###################################################################
    # Check identifiability and create mask
    ###################################################################

    sample_size = np.isfinite(statistics).all(axis=-2).sum(axis=-1)

    # pixels are 'valid' if there exists enough data such that the meta
    # regression / analysis model is identifiable.
    valid = sample_size > p+1

    # pixels are 'fully valid' if there are no missing values along the
    # subject axis, i.e, for all subjects in the sample there is an
    # estimate available for this pixel.
    fully_valid = np.isfinite(statistics).all(axis=(-1,-2))
    fully_valid = fully_valid & valid

    if mask is None:
        mask = valid
        print("  … setting the mask to default.")
    else:
        assert mask.dtype is np.dtype(bool), 'mask must be dtype: bool or None'
        assert mask.shape == valid.shape, 'shape of mask is wrong'

        if (~mask | valid).all():
            print('  … all points in the mask are identifiable.')
        else:
            print('  … not all points in the mask are identifiable.')

        if (~mask | fully_valid).all():
            print('  … no points with missing data along subject dimension.')
        else:
            print('  … points with missing data along subject dimension.')
        mask = valid & mask

    assert mask.any(), 'no identifiable points in the mask'

    ###################################################################
    # Need the squared standard error of the BOLD estimator
    ###################################################################

    stats = statistics.copy()
    stats[:,:,:,1] **= 2

    ###################################################################
    # Create result array
    ###################################################################

    result = np.zeros(mask.shape + (3, p+1), dtype=float)
    result[...] = np.nan

    ###################################################################
    # Reshape result, statistics, and mask
    ###################################################################

    rresult = result.reshape((-1,result.shape[-2],result.shape[-1]))
    rstats = stats.reshape((-1,stats.shape[-2],stats.shape[-1]))
    to_fit  = mask.reshape((-1,))

    print("  … number of points to estimate: {}".format(len(to_fit)))

    ###################################################################
    # Fit the model
    ###################################################################

    if design is None:
        print("  … perform a meta analysis")
        fit_meta_analysis(rresult, to_fit, rstats)
    else:
        print("  … perform a meta regression")
        fit_meta_regression(rresult, to_fit, rstats, design)

    return result

###################################################################
# Backends
###################################################################

def fit_meta_analysis(result, to_fit, statistics):
    for i in range(result.shape[0]):
        if to_fit[i]:
            b, t, h, adj_stderr, r = meta_analysis(
                    y = statistics[i,0][np.where(np.isfinite(statistics[i,0]))],
                    d = statistics[i,1][np.where(np.isfinite(statistics[i,1]))])
            result[i,0,0] = b
            result[i,1,0] = adj_stderr
            result[i,2,0] = t
            result[i,0,1] = h
            result[i,1,1] = r
            result[i,2,1] = 0

def fit_meta_regression(result, to_fit, statistics, design):
    for i in range(result.shape[0]):
        if to_fit[i]:
            valid = np.where(np.isfinite(statistics[i,0]))
            b, t, h, adj_stderr, r = meta_regression(
                    y = statistics[i,0][valid],
                    d = statistics[i,1][valid],
                    x = design[valid])
            result[i,0,:-1] = b
            result[i,1,:-1] = adj_stderr
            result[i,2,:-1] = t
            result[i,0,-1]  = h
            result[i,1,-1]  = r
            result[i,2,-1]  = 0
