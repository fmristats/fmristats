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
    formula : str
        A formula like object that is understood by patsy.
    design : ndarray
        If None, will be created from formula
        Directly specifying a design matrix, i.e., providing design
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

    def __init__(self, sample, formula=None, design=None):
        assert type(sample) is Sample, 'sample must be of type Sample'

        self.sample     = sample
        self.covariates = sample.covariates
        self.statistics = sample.statistics

        if (formula is not None) and (design is None):
            dmat = dmatrix(formula, self.covariates, eval_env=-1)
            parameter_names = dmat.design_info.column_names
            design = np.asarray(dmat)
        else:
            parameter_names = ['Intercept']

        self.formula = formula
        self.parameter_names = parameter_names
        self.design = design

    def fit(self, mask=True):
        """
        Fit the model to the data

        Returns
        -------
        PopulationResult
        """
        if mask is True:
            mask = 'vb'

        if mask is False:
            mask = None

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
                mask = None
                maskname = 'default'

        if mask is not None:
            assert type(mask) is np.ndarray, 'mask must be an ndarray'
            assert mask.any(), 'no valid points in mask'
            assert mask.dtype == bool, 'mask must be of dtype bool'

        self.mask = mask

        result = fit_field(
                statistics=self.statistics,
                design=self.design,
                mask=mask)

        return PopulationResult(statistics=result, model=self)

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

class PopulationResult:

    def __init__(self, statistics, model):
        self.statistics = statistics
        self.model      = model
        self.parameter_names = model.parameter_names

    def get_parameter(self, p=0):
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
        tstatistics = x[2,:-1]
        df = x[1,-1]
        pvalues = t.sf(tstatistics, df=df)

        df = pd.DataFrame({
            'parameter' : self.parameter_names + ['heterogeneity'],
            'point'     : x[0],
            'stderr'    : np.hstack((x[1,:-1], np.nan)),
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
