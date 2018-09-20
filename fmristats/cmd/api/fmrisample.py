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

from ...study import add_study_parser

def create_argument_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=hp.epilog)

    add_study_parser(parser)

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

from ...study import Study

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
    # Parse protocol
    ####################################################################

    df = get_df(args)

    if df is None:
        sys.exit()

    ####################################################################
    # Create study
    ####################################################################

    layout = {
        'irritation':args.irritation,
        'session':args.session,
        'reference_maps':args.reference_maps,
        'result':args.fit,
        'population_map':args.population_map}

    study = Study(df, df, layout=layout, strftime=args.strftime)

    study_iterator = study.iterate('result',
            vb_name=args.vb_name,
            diffeomorphism_name=args.diffeomorphism_name,
            scale_type=args.scale_type,
            integer_index = True)

    ####################################################################
    # Wrapper
    ####################################################################

    if not isfile(args.sample) or args.force:
        if args.verbose:
            print('Create population sample…'.format(args.sample))

        statistics = np.empty(vb.shape + (3,len(df),))
        statistics [ ... ] = np.nan

        for index, name, instances in study_iterator:
            result = instances['result']
            if result is not None:
                intercept = result.get_field('intercept','point')
                c = vb_ati.data / intercept.data

                beta = result.get_field('task','point')
                beta_stderr = result.get_field('task','stderr')

                statistics[...,0,index] = c*beta.data
                statistics[...,1,index] = c*beta_stderr.data
                statistics[...,2,index] = c
            else:
                study_iterator.df.ix[index,'valid'] = False

        df = study_iterator.df.copy()
        del df['result']

        print(df.head())

        sample = Sample(
                vb            = vb,
                vb_background = vb_background,
                vb_ati        = vb_ati,
                covariates    = study_iterator.df,
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
