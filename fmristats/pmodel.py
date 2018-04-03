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

from patsy import dmatrix

import numpy as np

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
        will take presendence from the formula interface.
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

        self.covariates = sample.covariates
        self.statistics = sample.statistics
        self.population_space = sample.population_space

        if (formula_like is not None) and (exog is None):
            exog = np.asarray(dmatrix(
                formula_like=formula_like,
                data=self.covariates))

        if mask is True:
            mask = 'template'

        if type(mask) is str:
            if mask == 'template':
                mask = self.population_space.get_mask()
            elif mask == 'sample':
                mask = np.isfinite(self.statistics).all(axis=(-1,-2))

        assert mask.any(), 'mask is empty, i.e., there exists no valid pixels in this mask'

        self.formula_like = formula_like
        self.exog = exog
        self.mask = mask

    def fit(self):
        """
        Fit the population model to data
        """
        statistics = fit_field(
                obs=self.statistics,
                design=self.exog,
                mask=self.mask)

        return MetaResult(
                statistics=statistics,
                formula_like=self.formula_like,
                exog=self.exog,
                population_space=self.population_space)

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
    def __init__(self, statistics, formula_like, exog, population_space):
        self.statistics = statistics
        self.formula = formula_like
        self.exog = exog
        self.population_space = population_space
        self.reference = population_space.reference

    def get_activation_patterns(self):
        return np.moveaxis(self.statistics[...,0,:-1], -1, 0)

    def get_tvalues(self):
        return np.moveaxis(self.statistics[...,1,:-1], -1, 0)

    def get_heterogeneities(self):
        return self.statistics[...,0,-1]

    def get_degrees_of_freedom(self):
        return self.statistics[...,1,-1]

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
