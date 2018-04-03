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

Sample

"""

from .load import load

from .diffeomorphisms import Image

import numpy as np

import pickle

from pandas import DataFrame

class Sample:
    """
    Sampled activation fields of a FMRI study
    """
    def __init__(self, population_space, covariates, statistics):
        """
        Parameters
        ----------
        population_space : Image
        covariates : DataFrame
        statistics : ndarray, shape (…,3)
        """
        assert type(population_space) is Image, 'population_space must be of type Image'
        self.population_space = population_space

        assert type(covariates) is DataFrame
        covariates.reset_index(inplace=True, drop=True)
        self.covariates = covariates

        assert statistics.shape == population_space.shape + (2,len(covariates))
        self.statistics = statistics

    def filter(self, b=None):
        """
        Notes
        -----
        Here, b should be a slice object, you cannot work with the index
        of the covariate data frame, but must use integer location
        indices instead.
        """
        if b is None:
            b = self.covariates.valid

        if b.dtype == np.dtype(bool):
            covariates = self.covariates.ix[b]
        else:
            covariates = self.covariates.iloc[b]

        return Sample(
                population_space = self.population_space,
                covariates       = covariates.ix[b],
                statistics       = self.statistics[...,b])

    ###################################################################
    # Description
    ###################################################################

    def describe(self):
        description = """
No. of samples: {:d}
{}
        """
        valid = self.covariates.groupby(
                ['cohort','paradigm','valid']).id.agg(['count'])
        return description.format(
                len(self.covariates), valid)

    ####################################################################
    # Save nstance to and from disk
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
