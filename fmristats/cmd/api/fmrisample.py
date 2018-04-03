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

Command line tool to create a sample of fits of the FMRI signal model

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

########################################################################
# Input arguments
########################################################################

# TODO: make population space optional, save ati in Sample and also the
# cfactor field!!

    parser.add_argument('sample',
            help="""path where to save the sample""")

    parser.add_argument('population_space',
            help="""path to a population space""")

    parser.add_argument('population_ati',
            help="""path to a population ATI reference""")

    parser.add_argument('--fit',
            default='../data/fit/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}-{5}.fit',
            help='output file;' + hp.sfit)

    parser.add_argument('--scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help=hp.scale_type)

########################################################################
# Output arguments
########################################################################

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-fmrisample.pkl',
            help=hp.protocol_log)

########################################################################
# Arguments specific for using the protocol API
########################################################################

    parser.add_argument('protocol',
            help=hp.protocol)

    parser.add_argument('--cohort',
            help=hp.cohort)

    parser.add_argument('--id',
            type=int,
            nargs='+',
            help=hp.j)

    parser.add_argument('--paradigm',
            help=hp.paradigm)

    parser.add_argument('--strftime',
            default='%Y-%m-%d-%H%M',
            help=hp.strftime)

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

from ...smodel import Result

from ...pmap import PopulationMap

from ...sample import Sample

from ...load import load, load_result

from ...name import Identifier

from ...protocol import layout_sdummy

import numpy as np

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

########################################################################

def call(args):

    try:
        population_space = load(args.population_space)
    except Exception as e:
        print('Unable to read: {}'.format(args.population_space))
        print('Exception: {}'.format(e))
        exit()

    try:
        population_ati = load(args.population_ati)
    except Exception as e:
        print('Unable to read: {}'.format(args.population_space))
        print('Exception: {}'.format(e))
        exit()

    ####################################################################

    output = args.protocol_log.format(
            datetime.datetime.now().strftime('%Y-%m-%d-%H%M'))

    if args.strftime == 'short':
        args.strftime = '%Y-%m-%d'

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args)

    if df is None:
        sys.exit()

    ####################################################################
    # Add file layout
    ####################################################################

    df.reset_index(inplace=True, drop=True)

    df_layout = df.copy()

    layout_sdummy(df_layout, 'file',
            template=args.fit,
            urname=population_space.name,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    ####################################################################
    # Apply wrapper
    ####################################################################

    if not isfile(args.sample) or args.force:
        if args.verbose:
            print('Create population sample…'.format(args.sample))

        statistics = np.empty(population_space.shape + (2,len(df),))
        statistics [ ... ] = np.nan

        for r in df_layout.itertuples():
            name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)

            wrapper(
                    name             = name,
                    df               = df,
                    index            = r.Index,
                    file             = r.file,
                    population_space = population_space,
                    population_ati   = population_ati,
                    statistics       = statistics,
                    verbose          = args.verbose,
                    )

        sample = Sample(
                population_space = population_space,
                covariates = df,
                statistics = statistics)

        sample = sample.filter()

        if args.verbose:
            print('Save: {}'.format(args.sample))

        sample.save(args.sample)

    else:
        if args.verbose:
            print('Parse: {}'.format(args.sample))
        sample = load(args.sample)

    if args.verbose:
        print('Description of the population space:')
        print(sample.population_space.describe())
        print('Description of the sample:')
        print(sample.describe())

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

def wrapper(name, df, index, file, population_space, population_ati,
        statistics, verbose):

    result = load_result(file, name, df, index, population_space.name, verbose)

    if not df.ix[index,'valid']:
        return

    intercept = result.get_field('intercept','point')
    c = population_ati.data / intercept.data

    beta = result.get_field('task','point')
    beta_stderr = result.get_field('task','stderr')

    statistics[...,0,index] = c*beta.data
    statistics[...,1,index] = (c*beta_stderr.data)**2

    #statistics[...,2,index] = c
