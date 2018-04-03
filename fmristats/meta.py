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
    t = b / np.sqrt(v * adj)
    return b, t, h, rdf

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
    """
    h   = hedge_type_estimator(y,d,x)
    fit = sm.WLS(y, x, weights=1/(h+d)).fit()
    v   = fit.bse
    b   = fit.params
    rdf = fit.df_resid
    adj = _qprofile_regression(h=h,y=y,d=d,x=x) / rdf
    t   = b / (v * np.sqrt(adj))
    return b, t, h, rdf

def fit_field(obs, design=None, mask=None):
    """
    Random effect meta regression

    Will regress each point of an effect and certainty field onto the
    covariates in x using a random effects meta regression model

    Parameters
    ----------
    obs : ndarray, shape 5D-image
        The observations.
    mask : ndarray, shape 3D-image
        A mask where to fit the field
    x : ndarray, shape (k,p)
        The covariates.

    Notes
    -----
    Results will be written into obs
    """
    if design is None:
        df = obs.shape[-1] - 1
        p  = 1
    else:
        df = -np.diff(design.shape)
        p  = design.shape[1]

    assert obs.shape[-1] > p+1, 'model not identifiable, too many parameters to fit'
    assert obs.shape[-1] > df, 'model not identifiable, too many parameters to fit'

    # pixels are 'fully valid' if there are no missing values along the
    # subject axis, i.e, for all subjects in the sample there is an
    # estimate available for this pixel.
    fully_valid = np.isfinite(obs).all(axis=(-1,-2))

    # pixels are 'valid' if there exists enough data such that the meta
    # regression / analysis model is identifiable.
    valid = np.isfinite(obs).all(axis=-2).sum(axis=-1)
    valid = valid > max(p+1,df)

    if mask is None:
        mask = valid
    else:
        assert mask.dtype is np.dtype(bool), 'mask must be of dtype bool or None'
        assert mask.shape == valid.shape, 'shape of mask does not fit'
        mask = mask & valid

    print("  … mask has shape: {}".format(mask.shape))

    assert mask.any(), 'mask is empty, i.e., there exists no valid pixels in this mask'

    # if the statement 'if valid, then fully_valid' holds, the following
    # loop don't need checks.  This is identical to 'not valid or
    # fully_valid'.
    skip_checks = (~mask | fully_valid).all()

    res = obs[mask]
    print("  … number of pixels to be fitted: {}".format(len(res)))

    it = iter(res)

    if skip_checks:
        print("  … all pixels in the mask are fully valid")
        if design is None:
            print("  … a meta analysis will be performed")
            for v in it:
                b, t, h, r = meta_analysis(y=v[0], d=v[1])
                v[0,0] = b
                v[1,0] = t
                v[0,1] = h
                v[1,1] = r
        else:
            print("  … a meta regression will be performed")
            for v in it:
                b, t, h, r = meta_regression(y=v[0], d=v[1], x=design)
                v[0,:p] = b
                v[1,:p] = t
                v[0,p] = h
                v[1,p] = r
    else:
        if design is None:
            print("  … a meta analysis will be performed")
            for v in it:
                sampled = ~np.isnan(v).any(axis=0)

                y=v[0,sampled]
                d=v[1,sampled]

                if (sampled).sum() > max(p+1,df):
                    b, t, h, r = meta_analysis(y=v[0], d=v[1])
                    v[0,0] = b
                    v[1,0] = t
                    v[0,1] = h
                    v[1,1] = r
                else:
                    v[...] = np.nan
        else:
            print("  … a meta regression will be performed")
            for v in it:
                sampled = ~np.isnan(v).any(axis=0)

                y=v[0,sampled]
                d=v[1,sampled]
                x=design[sampled]

                if (sampled).sum() > max(p+1,df):
                    b, t, h, r = meta_regression(y=v[0], d=v[1], x=x)
                    v[0,:p] = b
                    v[1,:p] = t
                    v[0,p] = h
                    v[1,p] = r
                else:
                    v[...] = np.nan

    result = np.zeros(mask.shape + (2,p+1,))
    result[...] = np.nan
    result[mask] = res[...,:p+1]

    return result
