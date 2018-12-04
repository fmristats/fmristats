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

Extract information from instance files and print information to
standard output

"""

########################################################################
#
# Command line program
#
########################################################################

from ..epilog import epilog

import argparse

def create_argument_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=epilog)

    parser.add_argument('files', nargs='+', help='input files')

    ####################################################################
    # Verbosity
    ####################################################################

    control_verbosity  = parser.add_argument_group(
        """Control the level of verbosity""")

    control_verbosity.add_argument('-v', '--verbose',
        action='count',
        default=0,
        help="""Increase output verbosity""")

    return parser

def cmd():
    parser = create_argument_parser()
    args = parser.parse_args()
    call(args)

cmd.__doc__ = __doc__

########################################################################
#
# Load libraries
#
########################################################################

from ..lock import Lock

from ..affines import Affine

from ..diffeomorphisms import Image, Diffeomorphism

from ..stimulus import Stimulus, Block

from ..session import Session

from ..reference import ReferenceMaps

from ..pmap import PopulationMap

from ..smodel import SignalModel, Result

from ..study import Study

from .. import load

from ..sample import Sample

from ..study import Study

from ..pmodel import PopulationModel, PopulationResult

import pandas as pd

from pandas import DataFrame

import numpy as np

########################################################################

def call(args):
    for f in args.files:
        try:
            print_info(load(f), f, args.verbose)
        except FileNotFoundError as e:
            print(e)

def print_info(x, f, verbose=False):
    if type(x) is Image:
        print('{}: image file'.format(f))
        print(x.describe())

    if type(x) is Session:
        print('{}: session file'.format(f))
        print(x.name.describe())
        print(x.describe())
        print(x.stimulus.describe())
        print(x.reference.describe())

    if type(x) is Stimulus:
        print('{}: stimulus file'.format(f))
        print(x.describe())

    if type(x) is Block:
        print('{}: block stimulus file'.format(f))
        print(x.describe())

    if type(x) is ReferenceMaps:
        print('{}: reference maps'.format(f))
        print(x.name.describe())
        print(x.describe())

    if type(x) is PopulationMap:
        print('{}: population map'.format(f))
        print(x.diffeomorphism.describe())
        print(x.describe())

    if type(x) is Result:
        print('{}: fit of a signal model'.format(f))
        print(x.describe())
        print(x.population_map.diffeomorphism.describe())
        print(x.population_map.describe())

    if type(x) is Lock:
        print('{}: lock'.format(f))
        print(x.describe())

    if type(x) is Sample:
        print('{}: sample file'.format(f))
        print(x.describe())

    if type(x) is PopulationModel:
        print('{}: population model file'.format(f))

    if type(x) is PopulationResult:
        print('{}: population result file'.format(f))

    def data_frame_information(x):
        print('Number of entries: {:d}'.format(len(x)))
        if len(x) <= 12:
            print('Entries:\n{}\n'.format(x))
        else:
            print('First twelve entries:\n{}\n'.format(x.head(12)))

        try:
            valid = x.groupby(['cohort','paradigm']).valid.agg(
                    ['sum', 'mean', 'count'])
        except:
            try:
                valid = x.groupby(['cohort', 'sex']).valid.agg(
                        ['sum', 'mean', 'count'])
            except:
                valid = x.groupby(['cohort']).valid.agg(
                        ['sum', 'mean', 'count'])

        valid['sum'] = valid['sum'].astype(int)
        valid['count'] = valid['count'].astype(int)

        #valid['mean'] = valid['mean'].astype(float)
        #valid['mean'] = ['{:.2f}'.format(m) for m in valid['mean']]

        valid.rename(
                columns={'sum':'valid', 'count':'total'},
                inplace=True)

        if x.valid.all():
            print('All entries are valid.\n')
        else:
            print('Number of valid entries:\n{}\n'.format(
                valid[['valid', 'total']]))

    if type(x) is Study:
        print('{}: study file\n'.format(f))
        print('Protocol:')
        print('---------')
        data_frame_information(x.protocol)
        if x.covariates is not None:
            print('Covariates:')
            print('-----------')
            data_frame_information(x.covariates)

        if verbose > 0:
            print()
            print('File layout:')
            print('------------')
            for k,v in x.file_layout.items():
                print('{:<15}: {}'.format(k, v))

    if type(x) is DataFrame:
        print('{}: protocol or covariates file:\n'.format(f))
        data_frame_information(x)
