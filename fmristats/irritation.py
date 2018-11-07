# Copyright 2015-2018 Thomas W. D. Möbius
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

Stimuluss, paradigms.

"""

from .name import Identifier

import numpy as np

import pandas as pd

from pandas import DataFrame

import pickle

class Stimulus:
    """
    Stimulus design

    Parameters
    ----------
    name : Identifier
        identifier of the session
    """

    def __init__(self, name):
        assert type(name) is Identifier, 'name must be of type Identifier'
        self.name = name

    def design(self, slice_time):
        return np.zeros_like(slice_time)

    def describe(self):
        """
        Give a description of the instance
        """
        return self.name.describe()

    def save(self, file, **kwargs):
        """
        Save stimulus instance to disk

        This will save the stimulus instance to disk for later use.

        Parameters
        ----------
        file : str
            File name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)

class Block(Stimulus):
    """
    Block stimulus of the subject

    Parameters
    ----------
    name : Identifier
        identifier of the session
    names : list
    onsets : dict
    durations : dict
    """

    def __init__(self, name, names, onsets, durations):
        assert type(name) is Identifier, 'name must be of type Identifier'

        number = len(names)

        assert len(names) == number, 'number of names does not match number of blocks'
        assert len(durations) == number, 'number of durations does not match number of blocks'
        assert len(onsets) == number, 'number of onsets does not match number of blocks'

        self.name = name
        self.number = number
        self.names = names
        self.durations = durations
        self.onsets = onsets

    def to_dataframe(self):
        # TODO: currently only works in balanced designs (equal number
        # of blocks in each condition!)
        df = DataFrame(self.onsets)
        df.columns.name = 'block'
        df = df.stack()
        df.name = 'onset'
        df = df.reset_index()
        del df['level_0']
        df.sort_values(by='onset', inplace=True)
        df['duration'] = 0.
        df = df.reset_index()
        del df['index']

        for c,d in df.groupby('block'):
            df.loc[df.block == c, 'duration'] = self.durations[c]

        return df

    def design(self, slice_timing, s, c, offset=0, preset=0):
        """
        Defines the stimulus design vector

        Parameters
        ----------
        slice_timing : ndarray, shape (n,) or (n,m)
            times at which slices have been measured during the
            experiment.  First dimension must be the same as the number
            of full scans
        s : str
            stimulus block of interest
        c : str
            name of the control block to be used
        offset : float
            offset to apply at the beginning of an stimulus phase
        preset : float
            offset to apply at the end of an stimulus phase

        Returns
        -------
        ndarray, shape_like(slice_timing), dtype: float
            The stimulus vector indicating the kind of stimulus under
            which the subject was exposed while the particular slice has
            been measured.
        """
        assert s != c, 'block names: {} and {c} must be different'.format(s,c)
        assert s in self.names, 'name of stimulus block: {} does not exist'.format(s)
        assert c in self.names, 'name of control block: {} does not exist'.format(c)

        design = np.empty(slice_timing.shape + (2,) )
        design[...] = np.nan

        onsets = self.onsets[c]
        durations = self.durations[c]
        srt = onsets + offset
        end = onsets + durations - preset
        mat = np.vstack((srt, end)).T

        group = 0
        for x in iter(mat):
            block = (slice_timing > x[0]) & (slice_timing < x[1])
            design[block,0] = 0
            design[block,1] = group
            group += 1

        onsets = self.onsets[s]
        durations = self.durations[s]
        srt = onsets + offset
        end = onsets + durations - preset
        mat = np.vstack((srt, end)).T

        group = 0
        for x in iter(mat):
            block = (slice_timing > x[0]) & (slice_timing < x[1])
            design[block,0] = 1
            design[block,1] = group
            group += 1

        return design

    def number_of_onsets(self):
        return {name : len(self.onsets[name]) for name in self.names}

    def irritated_scans(self, slice_time, offset=0, preset=0):
        """
        Test whether blocks have been aquired during blocks of stimulus

        Parameters
        ----------
        slice_time : ndarray, shape (n,) or (n,m), dtype: float

        Returns
        -------
        ndarray, shape_like(slice_time), dtype: bool
            If True, this scan has been aquired during a block of subject
            stimulus.
        """
        onsets = self.onsets.items()

        srt = np.hstack([v + offset for k,v in onsets])
        end = np.hstack([v + irr.durations[k] - preset for k,v in onsets])
        mat = np.vstack((srt, end)).T

        irritated = np.zeros_like(slice_time).astype(bool)
        for x in iter(mat):
            irritated = irritated | (slice_time > x[0]) & (slice_time < x[1])

        return irritated

    def describe(self):
        """
        Give a description of the instance
        """
        description = """
        Type of stimulus: block design
        Block number: {}
        Block names:  {}
        Number of onsets per block: {}"""
        return description.format(
                self.number,
                self.names,
                self.number_of_onsets(),
                )

