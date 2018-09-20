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

import pickle

def load_verbose(f, verbose=False, name=None):
    try:
        instance = load(f)
        if verbose:
            print('{}: Read {}'.format(name.name(), f))
        return instance
    except Exception as e:
        if verbose > 1:
            print('{}: Unable to read {}, {}'.format(name.name(), f, e))
        return None

class StudyIterator:
    def __init__(self, df, keys, new=None, verbose=True,
            integer_index=False):
        assert type(df) is DataFrame, 'df must be DataFrame'

        self.df = df
        self.keys = keys
        self.new = new
        self.verbose = verbose
        self.integer_index = integer_index

    def __iter__(self):
        df = self.df.reset_index(drop=True)
        self.it = df.itertuples()
        return self

    def __next__(self):
        r = next(self.it)
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)
        if self.new is None:
            if self.integer_index is False:
                return name, \
                    {k : load_verbose(getattr(r, k), self.verbose, name) for k in self.keys}
            else:
                return r.Index, name, \
                    {k : load_verbose(getattr(r, k), self.verbose, name) for k in self.keys}
        else:
            if self.integer_index is False:
                return name, \
                    {k : getattr(r, k) for k in self.new}, \
                    {k : load_verbose(getattr(r, k), self.verbose, name) for k in self.keys}
            else:
                return r.Index, name, \
                    {k : getattr(r, k) for k in self.new}, \
                    {k : load_verbose(getattr(r, k), self.verbose, name) for k in self.keys}

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

        self.layout = {
            'irritation' : '../data/irr/{2}/{0}-{1:04d}-{2}-{3}.irr',
            'session' : '../data/ses/{2}/{0}-{1:04d}-{2}-{3}.ses',
            'reference_maps' : '../data/ref/{2}/{0}-{1:04d}-{2}-{3}.ref',
            'population_map' : '../data/pop/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}.pop',
            'result' : '../data/fit/{2}/{4}/{5}/{6}/{0}-{1:04d}-{2}-{3}-{4}.fit',
            'strftime' : '%Y-%m-%d-%H%M'}

        if layout is not None:
            self.update_layout(layout)

        if strftime == 'short':
            strftime = '%Y-%m-%d'

        if strftime is not None:
            self.layout.update({'strftime' : strftime})

    def update_layout(self, layout):
        self.layout.update( (k,v) for k,v in layout.items() if v is not None)

    def iterate(self, *keys, new=None,
            vb_name=None, diffeomorphism_name=None, scale_type=None,
            verbose=True, integer_index=False):
        df = self.protocol.copy()

        for key in keys:
            df [key] = Series(
                    data = [self.layout[key].format(
                        r.cohort,
                        r.id,
                        r.paradigm,
                        r.date.strftime(self.layout['strftime']),
                        vb_name,
                        diffeomorphism_name,
                        scale_type)
                        for r in df.itertuples()],
                    index = df.index)

        if new is not None:
            for n in new:
                if n not in keys:
                    df [n] = Series(
                            data = [self.layout[n].format(
                                r.cohort,
                                r.id,
                                r.paradigm,
                                r.date.strftime(self.layout['strftime']),
                                vb_name,
                                diffeomorphism_name,
                                scale_type)
                                for r in df.itertuples()],
                            index = df.index)

        return StudyIterator(df, keys, new, verbose, integer_index)

    def save(self, file, **kwargs):
        """
        Save instance to disk

        This will save the current instance to disk for later use.

        Parameters
        ----------
        file : str
            File name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)

def add_study_parser(parser):

    parser.add_argument('--protocol',
            help="""A protocol file""")

    parser.add_argument('--covariates',
            help="""A covariate file""")

    parser.add_argument('--irritation',
            help="""Path to a irritation file or template for such a
            file""")

    parser.add_argument('--session',
            help="""Path to a session file or template for such a
            file""")

    parser.add_argument('--reference-maps',
            help="""Path to a reference maps or template for such a
            file""")

    parser.add_argument('--population-map',
            help="""Path to a population map file or template for such a
            file""")

    parser.add_argument('--fit',
            help="""Path to a result file or template for such a
            file""")

    parser.add_argument('--strftime',
            help="""Format of date and time""")

    parser.add_argument('--cohort',
            help="""Cohort""")

    parser.add_argument('--id',
            type=int,
            nargs='+',
            help="""id""")

    parser.add_argument('--datetime',
            help="""Datetime""")

    parser.add_argument('--paradigm',
            help="""Paradigm""")

    parser.add_argument('--vb-name',
            default='self',
            help="""Name of the population space""")

    parser.add_argument('--diffeomorphism-name',
            default='identity',
            help="""Name of the diffeomorphism between population space
            and subject space""")

    parser.add_argument('--scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help="""Scale type""")
