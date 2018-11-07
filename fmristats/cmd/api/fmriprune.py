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

Prune

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
# Arguments specific for the application of masks
########################################################################

    handling_df_cutoff = parser.add_mutually_exclusive_group()

    handling_df_cutoff.add_argument('-p', '--proportion',
            type=float,
            default=.6,
            help="""estimates which degrees of freedom are below the
            proportional threshold of the degrees of freedom in the
            effect field estimate are set to null.""")

    handling_df_cutoff.add_argument('-t', '--threshold',
            type=int,
            help="""estimates which degrees of freedom are below the
            threshold of the degrees of freedom in the effect field
            estimate are set to null.""")

########################################################################
# Miscellaneous
########################################################################

    parser.add_argument('-f', '--force',
            action='store_true',
            help="""Overwrite mask if it already exits""")

    parser.add_argument('-v', '--verbose',
            action='count',
            default=0,
            help=hp.verbose)

########################################################################
# Multiprocessing
########################################################################

    parser.add_argument('-j', '--cores',
            type=int,
            help=hp.cores)

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

from ...load import load_result

from ...name import Identifier

from ...protocol import layout_dummy, layout_sdummy

from ...smodel import Result

from ...study import Study

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

import nibabel as ni

########################################################################

def call(args):

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args, fall_back=args.fit)

    if df is None:
        sys.exit()

    ####################################################################
    # Add file layout
    ####################################################################

    layout = {
        'stimulus':args.stimulus,
        'session':args.session,
        'reference_maps':args.reference_maps,
        'result':args.fit,
        'population_map':args.population_map,
        'diffeomorphism_name':args.diffeomorphism_name,
        'scale_type':args.scale_type,
        }

    study = Study(df, df, layout=layout, strftime=args.strftime)

    study_iterator = study.iterate('result', new=['result'],
            vb_name=args.vb_name,
            diffeomorphism_name=args.diffeomorphism_name,
            scale_type=args.scale_type)

    ####################################################################
    # wrapper
    ####################################################################

    def wm(result, filename, name):
            verbose        = args.verbose
            threshold_df   = args.threshold
            proportion_df  = args.proportion

            if verbose > 1:
                print('{}: Description of the fit:'.format(name.name()))
                print(result.describe())
                print(result.population_map.describe())

            gf = result.get_field('degrees_of_freedom')

            if proportion_df:
                threshold_df = int(proportion_df * np.nanmax(gf.data))

            if verbose:
                print('{}: Lower df threshold: {:d}'.format(name.name(), threshold_df))

            inside = (gf.data >= threshold_df)

            result.population_map.set_vb_mask(inside)

            result.population_map.vb_mask.name = 'pruned_intercept'

            if verbose:
                print('{}: Save: {}'.format(name.name(), filename))

            try:
                result.save(filename)
            except Exception as e:
                print('{}: Unable to write: {}'.format(name.name(), filename))
                print('{}: Exception: {}'.format(name.name(), e))

    ###################################################################

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for name, files, instances in study_iterator:
                result   = instances['result']
                if result is not None:
                    filename = files['result']
                    pool.apply_async(wm, args=(result, filename, name))

            pool.close()
            pool.join()
        except Exception as e:
            pool.close()
            pool.terminate()
            print('Pool execution has been terminated')
            print(e)
        finally:
            pass
    else:
        try:
            print('Process protocol entries sequentially')
            for name, files, instances in study_iterator:
                result   = instances['result']
                if result is not None:
                    filename = files['result']
                    wm(result, filename, name)
        finally:
            pass
