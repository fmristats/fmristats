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

from ..api.fmriprune import create_argument_parser as cap

import argparse

def create_argument_parser():
    parser = cap()

    parser.add_argument('--cmd-bet',
            default='fsl5.0-bet',
            help="""FSL bet command. Must be in your path.""")

    parser.add_argument('--variante',
            default='R',
            help="""FSL bet command. Must be in your path.""")

    parser.add_argument('--vb-file',
            default='pruning/{0}-{1:04d}-{2}-{3}-{4}-{5}-intercept.nii.gz',
            help="""brain mask in image space.""")

    parser.add_argument('--vb-mask',
            default='pruning/{0}-{1:04d}-{2}-{3}-{4}-{5}-mask.nii.gz',
            help="""brain mask in image space.""")

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

from ...nifti import image2nii

from ...fsl import bet

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

    layout_sdummy(df_layout, 'vb_file',
            template=args.vb_file,
            urname=args.population_space,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'vb_mask',
            template=args.vb_mask,
            urname=args.population_space,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

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
                intercept_file        = r.vb_file,
                mask_file             = r.vb_mask,
                cmd_bet               = args.cmd_bet,
                variante              = args.variante
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

def wrapper(name, df, index, filename, verbose, force, vb, threshold_df,
        proportion_df, intercept_file, mask_file, cmd_bet, variante):

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

    gf = result.get_field('degrees_of_freedom')

    if proportion_df:
        threshold_df = int(proportion_df * np.nanmax(gf.data))

    if verbose:
        print('{}: Lower df threshold: {:d}'.format(name.name(), threshold_df))

    inside = (gf.data >= threshold_df)

    intercept = result.get_field('intercept', 'point')
    intercept = intercept.mask(inside)
    intercept = intercept.round()

    if verbose:
        print('{}: Fit brain mask'.format(name.name()))

    template = bet(
            intercept = intercept,
            intercept_file = intercept_file,
            mask_file = mask_file,
            cmd = cmd_bet,
            variante = variante,
            verbose = verbose)

    if template is None:
        df.ix[index,'valid'] = False
        print('{}: Unable to bet'.format(name.name()))
        lock.conditional_unlock(df, index, verbose, True)
        return

    result.population_map.set_template_mask(template)

    if verbose:
        print('{}: Save: {}'.format(name.name(), filename))

    try:
        result.save(filename)
        df.ix[index,'locked'] = False
    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to write: {}'.format(name.name(), filename))
        print('{}: Exception: {}'.format(name.name(), e))
