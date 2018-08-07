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

Defines a class for the FMRI population model and its fits.

"""

from .sample import Sample

from .meta import fit_field

from .diffeomorphisms import Image

from patsy import dmatrix

import numpy as np

from scipy.stats.distributions import t

import pandas as pd

import pickle

class PopulationModel:
    """
    The FMRI population model

    Parameters
    ----------
    sample : Sample
    formula_like : str
        A formula like object that is understood by patsy.
    exog : ndarray
        If None, will be created from formula_like
        Directly specifying a design matrix, i.e., providing exog
        will take precedence from the formula interface.
    mask : None, bool, ndarray, or str
        If False or None, the population model will be fitted only
        at points, at which the population model is identifiable.
        If True or 'template', the population model will,
        additionally, only be fitted at points, at which the
        template in the population has valid intensities (i.e. > 0
        and not NAN). If 'sample', the population model will be
        fitted only at points at which the sample provides valid
        summary statistics for *all* fields in the sample.
    """

    def __init__(self, sample, formula_like=None, exog=None, mask=True):
        assert type(sample) is Sample, 'sample must be of type Sample'

        self.sample     = sample
        self.covariates = sample.covariates
        self.statistics = sample.statistics

        if (formula_like is not None) and (exog is None):
            exog = np.asarray(dmatrix(
                formula_like=formula_like,
                data=self.covariates))

        datamask = Image(
            self.sample.vb.reference,
            np.isfinite(self.sample.statistics).all(axis=-2).sum(axis=-1))

        datamask.mask()

        datamask = datamask.get_mask()

        if mask is True:
            mask = 'vb'

        if type(mask) is str:
            if mask == 'vb':
                mask = self.sample.vb.get_mask()
                massname = 'template (vb)'
            elif mask == 'vb_background':
                mask = self.sample.vb_background.get_mask()
                maskname = 'template background (vb_background)'
            elif mask == 'vb_estimate':
                mask = self.sample.vb_estimate.get_mask()
                maskname = 'template estimate (vb_estimate)'
            else:
                mask = datamask
                maskname = 'data mask (foreground/background)'

        assert mask.any(), 'mask is empty, i.e., there exist no valid pixels in this mask'

        if mask is not None:
            assert type(mask) is np.ndarray, 'mask must be an ndarray'
            assert mask.dtype == bool, 'mask must be of dtype bool'
            assert mask.shape == datamask.shape, 'shape of mask and image must match'
            mask = mask & datamask

        self.formula_like = formula_like
        self.exog = exog
        self.mask = mask

    def fit(self):
        """
        Fit the population model to data
        """
        statistics, p, parameter_names = fit_field(
                obs=self.statistics,
                design=self.exog,
                mask=self.mask)

        return MetaResult(statistics=statistics, model=self, p=p,
                parameter_names=parameter_names)

    ####################################################################
    # Save instance to and from disk
    ####################################################################

    def save(self, file, **kwargs):
        """
        Save model instance to disk

        This will save the current model instance to disk for later use.

        Parameters
        ----------
        file : str
            File name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)

class MetaResult:

    def __init__(self, statistics, model, p, parameter_names):
        self.statistics = statistics
        self.model      = model
        self.p          = p
        self.parameter_names = parameter_names

    def get_parameter(self):
        f = np.moveaxis(self.statistics[...,0,:-1], -1, 0)
        return Image(data=f[p], reference=self.model.sample.vb.reference)

    def get_tstatistic(self, p=0):
        f = np.moveaxis(self.statistics[...,2,:-1], -1, 0)
        return Image(data=f[p], reference=self.model.sample.vb.reference)

    def get_heterogeneity(self):
        f = self.statistics[...,0,-1]
        return Image(data=f, reference=self.model.sample.vb.reference)

    def get_degree_of_freedom(self):
        f = self.statistics[...,1,-1]
        return Image(data=f, reference=self.model.sample.vb.reference)

    def at_index(self, index):
        # TODO: also add a stderr to the h-estimate (Knapp-Hartung!)

        x  = self.statistics[index]
        tstatistics = x[2,:self.p]
        df = x[1,-1]
        pvalues = 1 - t.cdf(tstatistics, df=df)

        df = pd.DataFrame({
            'parameter' : self.parameter_names + ['heterogeneity'],
            'point'     : x[0],
            'stderr'    : np.hstack((x[1,:self.p], np.nan)),
            'tstatistic': np.hstack((tstatistics, np.nan)),
            'df'        : df,
            'pvalue'    : np.hstack((pvalues,np.nan))})

        return df

    ####################################################################
    # Save instance to and from disk
    ####################################################################

    def save(self, file, **kwargs):
        """
        Save model instance to disk

        This will save the current model instance to disk for later use.

        Parameters
        ----------
        file : str
            File name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)
