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

Extract summary field statistics of the effect field of a fitted
Signal Model and save the statistics as a DataFrame.

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

    parser.add_argument('df',
            help="""output file""")

########################################################################
# Input arguments
########################################################################

    parser.add_argument('--fit',
            default='../data/fit/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}-{5}.fit',
            help='input file;' + hp.sfit)

########################################################################
# Output arguments
########################################################################

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-fit2df.pkl',
            help=hp.protocol_log)

########################################################################
# Arguments specific for using the protocol API
########################################################################

    parser.add_argument('--protocol',
            help=hp.protocol)

    parser.add_argument('--cohort',
            help=hp.cohort)

    parser.add_argument('--id',
            type=int,
            nargs='+',
            help=hp.j)

    parser.add_argument('--datetime',
            help=hp.datetime)

    parser.add_argument('--paradigm',
            help=hp.paradigm)

    parser.add_argument('--strftime',
            default='%Y-%m-%d-%H%M',
            help=hp.strftime)

    parser.add_argument('--population-space',
            default='reference',
            help=hp.population_space)

    parser.add_argument('--scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help=hp.scale_type)

########################################################################
# Miscellaneous
########################################################################

    parser.add_argument('-f', '--force',
            action='store_true',
            help=hp.force.format('result'))

    parser.add_argument('-v', '--verbose',
            action='store_true',
            help=hp.verbose)

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

from ..df import get_df

from ...lock import Lock

from ...load import load_result

from ...name import Identifier

from ...protocol import layout_dummy, layout_sdummy

from ...smodel import Result

import pandas as pd

from pandas import DataFrame

import datetime

import sys

import os

from os.path import isfile, isdir, join

import numpy as np

########################################################################

def call(args):
    if isfile(args.df) and not args.force:
        print('File {} already exists, use -f/--force to overwrite'.format(args.df))
        sys.exit()

    output = args.protocol_log.format(
            datetime.datetime.now().strftime('%Y-%m-%d-%H%M'))

    if args.strftime == 'short':
        args.strftime = '%Y-%m-%d'

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args, fall_back=args.fit)

    if df is None:
        sys.exit()

    ####################################################################
    # Add file layout
    ####################################################################

    df_layout = df.copy()

    layout_sdummy(df_layout, 'file',
            template=args.fit,
            urname=args.population_space,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    ####################################################################
    # Apply wrapper
    ####################################################################

    summary = []

    if args.verbose:
        print('Extract statistics'.format(output))

    for r in df_layout.itertuples():
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)

        tmp_summary = wrapper(
                name     = name,
                df       = df,
                index    = r.Index,
                filename = r.file,
                verbose  = args.verbose,
                vb       = args.population_space,
                )

        if tmp_summary is not None:
            summary += [tmp_summary]

    ####################################################################
    # Write summaries
    ####################################################################

    if args.verbose:
        print('Concatenate statistics'.format(args.df))

    summary = pd.concat(summary, ignore_index=True)

    if args.verbose:
        print('Save: {}'.format(args.df))

    dfile = os.path.dirname(args.df)
    if dfile and not isdir(dfile):
       os.makedirs(dfile)

    summary.to_pickle(args.df)

    ####################################################################
    # Write protocol
    ####################################################################

    if args.verbose:
        print('Save: {}'.format(output))

    dfile = os.path.dirname(output)
    if dfile and not isdir(dfile):
       os.makedirs(dfile)

    df.to_pickle(output)

########################################################################

def wrapper(name, df, index, filename, verbose, vb):

    ####################################################################
    # Load fit from disk
    ####################################################################

    result = load_result(filename, name, df, index, vb, verbose)
    if df.ix[index,'valid'] == False:
        return

    if verbose > 1:
        print('{}: Description of the fit:'.format(result.name.name()))
        print(result.describe())

    result.mask(True, verbose)

    tsimg = result.get_field('task', 'tstatistic')
    dwimg = result.get_field('durbin_watson')

    df = DataFrame({
        'cohort'   : result.name.cohort,
        'id'       : result.name.j,
        'date'     : result.name.datetime,
        'paradigm' : result.name.paradigm,
        'name'     : result.name.name(),
        't'        : tsimg.data[dwimg.get_mask()],
        'dw'       : dwimg.data[dwimg.get_mask()]
        })

    return df
