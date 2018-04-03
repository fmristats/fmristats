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

import fmristats.cmd.hp as hp

import argparse

def create_argument_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=hp.epilog)

    parser.add_argument('files', nargs='+', help='input files')

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

from ...lock import Lock

from ...affines import Affine

from ...diffeomorphisms import Image, Diffeomorphism

from ...irritation import Irritation, Block

from ...session import Session

from ...reference import ReferenceMaps

from ...pmap import PopulationMap

from ...smodel import SignalModel, Result

from ... import load

from ...sample import Sample

#from ...pmodel import PopulationModel, MetaResult

import pandas as pd

from pandas import DataFrame

import numpy as np

########################################################################

def call(args):
    for f in args.files:
        print_info(load(f), f)

def print_info(x, f):
    if type(x) is Image:
        print('{}: image file'.format(f))
        print(x.describe())

    if type(x) is Session:
        print('{}: session file'.format(f))
        print(x.name.describe())
        print(x.describe())
        print(x.irritation.describe())
        print(x.reference.describe())

    if type(x) is Irritation:
        print('{}: irritation file'.format(f))
        print(x.describe())

    if type(x) is Block:
        print('{}: block irritation file'.format(f))
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

    if type(x) is DataFrame:
        if 'cohort' in x.columns:
            x.cohort.cat.remove_unused_categories(inplace=True)
            if 'valid' in x.columns:
                if 'paradigm' in x.columns:
                    x.paradigm.cat.remove_unused_categories(inplace=True)
                    valid = x.groupby(['cohort','paradigm','valid']).id.agg(['count'])
                else:
                    valid = x.groupby(['cohort','valid']).id.agg(['count'])
        else:
            valid = None

        print("""{}: protocol file

No. of entries: {:d}

First five entries:
{}

Data types:
{}

Valid entries:
{}""".format(f, len(x), x.head(), x.dtypes, valid))

