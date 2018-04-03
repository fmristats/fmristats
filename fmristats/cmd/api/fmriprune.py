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

def create_argument_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=hp.epilog)

########################################################################
# Input arguments
########################################################################

    parser.add_argument('--fit',
            default='../data/fit/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}-{5}.fit',
            help='input file;' + hp.sfit)

########################################################################
# Output arguments
########################################################################

    parser.add_argument('--pruned',
            default='../data/fit/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}-{5}.fit',
            help='input file;' + hp.sfit)

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-fsl4prune.pkl',
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
# Arguments specific for the application of masks
########################################################################

    handling_df_cutoff = parser.add_mutually_exclusive_group()

    handling_df_cutoff.add_argument('-p', '--proportion',
            type=float,
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

    layout_sdummy(df_layout, 'filename',
            template=args.fit,
            urname=args.population_space,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    ####################################################################
    # Apply wrapper
    ####################################################################

    def wm(r):
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)

        wrapper(name                  = name,
                df                    = df,
                index                 = r.Index,
                filename              = r.filename,
                verbose               = args.verbose,
                force                 = args.force,
                vb                    = args.population_space,
                threshold_df          = args.threshold,
                proportion_df         = args.proportion,
                )

    it =  df_layout.itertuples()

    if len(df_layout) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            results = pool.map(wm, it)
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
            for r in it:
                wm(r)
        finally:
            pass

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

def wrapper(name, df, index, filename, verbose, force, vb, threshold_df, proportion_df):

    result = load_result(filename, name, df, index, vb, verbose)
    if df.ix[index,'valid'] == False:
        return

    if verbose > 1:
        print('{}: Description of the fit:'.format(result.name.name()))
        print(result.describe())

    if hasattr(result.population_map, 'template_mask') and \
        (result.population_map.template_mask is not None) and not force:
        if verbose:
            print('{}: Template mask already exits. Use -f/--force to overwrite'.format(name.name()))
        return

    gf = result.get_field('df', 'other')

    if proportion_df:
        threshold_df = int(proportion_df * np.nanmax(gf.data))

    if verbose:
        print('{}: Lower df threshold: {:d}'.format(name.name(), threshold_df))

    result.population_map.set_template_mask(gf.data >= threshold_df)

    if verbose:
        print('{}: Save: {}'.format(name.name(), filename))

    try:
        result.save(filename)
        df.ix[index,'locked'] = False
    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to write: {}'.format(name.name(), filename))
        print('{}: Exception: {}'.format(name.name(), e))
