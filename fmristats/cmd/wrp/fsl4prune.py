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

from ...load import load_result

from ...name import Identifier

from ...protocol import layout_dummy, layout_sdummy

from ...smodel import Result

from ...study import Study

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

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args, fall_back=args.fit)

    if df is None:
        sys.exit()

    ####################################################################
    # Add file layout
    ####################################################################

    study = Study(df,df,strftime=args.strftime)

    study.layout.update({'vb_file':args.vb_file, 'vb_mask':args.vb_mask})

    study_iterator = study.iterate('result',
            new=['result', 'vb_file', 'vb_mask'],
            vb_name=args.vb_name,
            diffeomorphism_name=args.diffeomorphism_name,
            scale_type=args.scale_type)

    ####################################################################
    # Wrapper
    ####################################################################

    def wm(result, vb_file, vb_mask, filename, name):
            verbose        = args.verbose
            threshold_df   = args.threshold
            proportion_df  = args.proportion
            cmd_bet        = args.cmd_bet
            variante       = args.variante

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

            intercept = result.get_field('intercept', 'point')
            intercept.mask(inside)
            intercept = intercept.round()

            if verbose:
                print('{}: Fit brain mask'.format(name.name()))

            template = bet(
                    intercept = intercept,
                    intercept_file = vb_file,
                    mask_file = vb_mask,
                    cmd = cmd_bet,
                    variante = variante,
                    verbose = verbose)

            if template is None:
                print('{}: Unable to bet'.format(name.name()))
                return

            result.population_map.set_vb_mask(gf.data >= threshold_df)

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
                    vb_file  = files['vb_file']
                    vb_mask  = files['vb_mask']
                    filename = files['result']
                    pool.apply_async(wm, args=(result, vb_file, vb_mask, filename, name))

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
                    vb_file  = files['vb_file']
                    vb_mask  = files['vb_mask']
                    filename = files['result']
                    wm(result, vb_file, vb_mask, filename, name)
        finally:
            pass
