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

Sample

"""

from .load import load

from .diffeomorphisms import Image

from .study import Study

import numpy as np

import pickle

from pandas import DataFrame

class Sample:
    """
    Sampled activation fields of a FMRI study
    """
    def __init__(self, covariates, statistics, study):
        """
        Parameters
        ----------
        covariates : DataFrame
            Should have an integer index (with start=0,
            stop=len(covariates)) of the same length as the number of
            statistics that are saved in statistics.
        statistics : ndarray, shape (…,3)
            The statistics field.
        study : Study
            The study.
        """
        assert type(study) is Study, 'study must be of type Study'

        self.study = study
        self.vb = study.vb
        self.vb_ati = study.vb_ati

        assert type(covariates) is DataFrame
        self.covariates = covariates

        assert statistics.shape == self.vb.shape + (3,len(covariates))
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
                covariates = covariates.ix[b],
                statistics = self.statistics[...,b],
                study      = self.study)

    def at_index(self, index):
        """
        Returns the summary statistics at an index

        Returns the effect and respective standard error of all subjects
        in the sample.

        Parameters
        ----------
        index : tuple(int)
            The index. A 3-tuple of integers

        Returns
        -------
        DataFrame
        """
        df = self.covariates.copy()
        df['task'] = self.statistics[index[0],index[1],index[2],0]
        df['stderr'] = self.statistics[index[0],index[1],index[2],1]
        return df

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

    def __str__(self):
        return self.describe()
