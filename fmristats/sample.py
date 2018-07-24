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
    def __init__(self, vb, vb_background, vb_ati, covariates, statistics):
        """
        Parameters
        ----------
        vb : Image
        vb_background : Image
        vb_ati : Image
        covariates : DataFrame
        statistics : ndarray, shape (…,3)
        """
        assert type(vb) is Image, 'vb must be of type Image'
        self.vb = vb

        assert type(vb_background) is Image, 'vb_background must be of type Image'
        self.vb_background = vb_background

        assert type(vb_ati) is Image, 'vb_ati must be of type Image'
        self.vb_ati = vb_ati

        assert type(covariates) is DataFrame
        covariates.reset_index(inplace=True, drop=True)
        self.covariates = covariates

        assert statistics.shape == vb.shape + (3,len(covariates))
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
                vb               = self.vb,
                vb_background    = self.vb_background,
                vb_ati           = self.vb_ati,
                covariates       = covariates.ix[b],
                statistics       = self.statistics[...,b])

    def at_index(self, index):

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
