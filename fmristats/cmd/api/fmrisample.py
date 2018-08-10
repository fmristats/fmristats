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

# TODO: make population space optional, save ati in Sample and also
# save the cfactor field!!

    parser.add_argument('sample',
            help="""path where to save the sample""")

    parser.add_argument('vb',
            help="""path to an image in population space""")

    parser.add_argument('vb_background',
            help="""path to a background image in population space""")

    parser.add_argument('vb_ati',
            help="""path to an ATI reference field in population space""")

    parser.add_argument('--vb-name',
            help="""name if different to the name saved in vb""")

    parser.add_argument('--fit',
            default='../data/fit/{2}/{4}/{5}/{6}/{0}-{1:04d}-{2}-{3}-{4}.fit',
            help=hp.sfit)

    parser.add_argument('--diffeomorphism',
            default='ants',
            help="""Name of the fitted diffeomorphisms.""")

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

from ...protocol import layout_sdummy, layout_fdummy

import numpy as np

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

########################################################################

def call(args):

    try:
        vb = load(args.vb)
    except Exception as e:
        print('Unable to read: {}'.format(args.vb))
        print('Exception: {}'.format(e))
        exit()

    try:
        vb_background = load(args.vb_background)
    except Exception as e:
        print('Unable to read: {}'.format(args.vb_background))
        print('Exception: {}'.format(e))
        exit()

    try:
        vb_ati = load(args.vb_ati)
    except Exception as e:
        print('Unable to read: {}'.format(args.vb_ati))
        print('Exception: {}'.format(e))
        exit()

    if args.vb_name is None:
        vb_name = vb.name
    else:
        vb_name = args.vb_name

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

    layout_fdummy(df_layout, 'file',
            template=args.fit,
            vb=vb_name,
            diffeo=args.diffeomorphism,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    ####################################################################
    # Apply wrapper
    ####################################################################

    if not isfile(args.sample) or args.force:
        if args.verbose:
            print('Create population sample…'.format(args.sample))

        statistics = np.empty(vb.shape + (3,len(df),))
        statistics [ ... ] = np.nan

        for r in df_layout.itertuples():
            name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)

            wrapper(
                    name             = name,
                    df               = df,
                    index            = r.Index,
                    file             = r.file,
                    vb               = vb,
                    vb_ati           = vb_ati,
                    vb_name          = vb_name,
                    statistics       = statistics,
                    verbose          = args.verbose,
                    )

        sample = Sample(
                vb            = vb,
                vb_background = vb_background,
                vb_ati        = vb_ati,
                covariates    = df,
                statistics    = statistics)

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
        print(sample.vb.describe())
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

def wrapper(name, df, index, file, vb, vb_ati, vb_name,
        statistics, verbose):

    result = load_result(file, name, df, index, vb_name, verbose)

    if not df.ix[index,'valid']:
        return

    intercept = result.get_field('intercept','point')
    c = vb_ati.data / intercept.data

    beta = result.get_field('task','point')
    beta_stderr = result.get_field('task','stderr')

    statistics[...,0,index] = c*beta.data
    statistics[...,1,index] = c*beta_stderr.data
    statistics[...,2,index] = c
