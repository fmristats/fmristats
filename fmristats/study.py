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

Study Layout

"""

from .name import Identifier

from .load import load, load_block_irritation, load_session, load_refmaps, \
    load_population_map, load_result

import pandas as pd

from pandas import Series, DataFrame

import os

from os.path import isfile, isdir, join

def load_verbose(f, verbose=False, name=None):
    try:
        instance = load(f)
        if verbose:
            print('{}: Read {}'.format(name.name(), f))
        return instance
    except Exception as e:
        print('{}: Unable to read {}, {}'.format(name.name(), f, e))
        return None

class StudyIterator:
    def __init__(self, df, keys, verbose=True):
        assert type(df) is DataFrame, 'df must be DataFrame'

        self.df = df
        self.keys = keys
        self.verbose = verbose
        iter(self)

    def __iter__(self):
        self.it = self.df.itertuples()
        return self

    def __next__(self):
        r = next(self.it)
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)
        return {k : load_verbose(getattr(r, k), self.verbose, name) for k in self.keys}

class Study:
    """
    Study Layout
    """

    def __init__(self,
            protocol,
            covariates,
            vb=None,
            vb_background=None,
            vb_ati=None,
            layout=None,
            strftime=None,

            ):
        """
        Parameters
        ----------
        covariates : DataFrame
        statistics : ndarray, shape (…,3)
        vb : Image
        vb_background : Image
        vb_ati : Image
        irritation : str
        session : str
        reference_maps : str
        population_map : str
        result : str
        strftime : str
        """
        self.protocol        = protocol
        self.covariates      = covariates
        self.vb              = vb
        self.vb_background   = vb_background
        self.vb_ati          = vb_ati

        self.file_layout = {
            'irritation' : '../data/irr/{2}/{0}-{1:04d}-{2}-{3}.irr',
            'session' : '../data/ses/{2}/{0}-{1:04d}-{2}-{3}.ses',
            'reference_maps' : '../data/ref/{2}/{0}-{1:04d}-{2}-{3}.ref',
            'population_map' : '../data/pop/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}.pop',
            'result' : '../data/fit/{2}/{4}/{5}/{6}/{0}-{1:04d}-{2}-{3}-{4}.fit',
            'strftime' : '%Y-%m-%d-%H%M'}

        if layout is not None:
            self.layout = self.file_layout.update(layout)

        if strftime == 'short':
            strftime = '%Y-%m-%d'

        if strftime is not None:
            self.layout = self.file_layout.update({'strftime' : strftime})

        #self.irritation      = irritation
        #self.session         = session
        #self.reference_maps  = reference_maps
        #self.population_map  = population_map
        #self.result          = result
        #self.strftime        = strftime

    def iterate(self, *keys,
            vb_name=None, diffeomorphism_name=None, scale_type=None):
        df = self.protocol.copy()

        for key in keys:
            df [key] = Series(
                    data = [self.file_layout[key].format(
                        r.cohort,
                        r.id,
                        r.paradigm,
                        r.date.strftime(self.file_layout['strftime']),
                        vb_name,
                        diffeomorphism_name,
                        scale_type)
                        for r in df.itertuples()],
                    index = df.index)

        return StudyIterator(df, keys)

