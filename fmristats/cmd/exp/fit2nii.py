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

Export fit to Nifti1

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

    parser.add_argument('--nii',
            default='../results/nii/{2}/{0}-{1:04d}-{2}-{3}-{4}-{5}-{{}}-{{}}.nii.gz',
            help="""output file""")

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-fit2nii.pkl',
            help=hp.protocol_log)

########################################################################
# Arguments specific for the RSM Signal Model: where to
########################################################################

    where_to_fit = parser.add_argument_group(
            """where to""",
            """By default, fit2nii will respect the brain mask which is
            saved in the population map. This default behaviour can, of
            course, be changed.""")

    where_to_fit.add_argument('--ignore-mask',
            action='store_true',
            help="""Ignore any brain masks saved in the respective
            population map.""")

    # TODO: mask: not implemented yet
    #where_to_fit.add_argument('--mask',
    #        action='store_true',
    #        help=hp.population_mask)

########################################################################
# Arguments specific for using the protocol API
########################################################################

    to_process = parser.add_argument_group(
            """specifying the protocol entries to process""",
            """Arguments which give control which protocol entries to
            process. If no protocol file is given, it will be checked
            if files being processed comply to the given information.""")

    to_process.add_argument('--protocol',
            help=hp.protocol)

    to_process.add_argument('--cohort',
            help=hp.cohort)

    to_process.add_argument('--id',
            type=int,
            nargs='+',
            help=hp.j)

    to_process.add_argument('--datetime',
            help=hp.datetime)

    to_process.add_argument('--paradigm',
            help=hp.paradigm)

    to_process.add_argument('--strftime',
            default='%Y-%m-%d-%H%M',
            help=hp.strftime)

    which_to_process = parser.add_argument_group(
            """specifying fit and population space""",
            """If in-  and output files are specified as templates, set
            some of the fixed terms in the template.""")

    which_to_process.add_argument('--scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help=hp.scale_type)

    which_to_process.add_argument('--population-space',
            default='reference',
            help=hp.population_space)

########################################################################
# Arguments which field to extract
########################################################################

    parser.add_argument('--parameter',
            type=str,
            nargs='+',
            default=['intercept', 'activation'],
            help="""which parameter field or which parameter fields to
            plot.  Available options are `intercept`, `activation`, or
            `time`. If you want to plot more than one parameter field,
            then separate these by a space.""")

    parser.add_argument('--value',
            type=str,
            nargs='+',
            default=['point', 'tstatistic'],
            #default=['point', 'stderr', 'tstatistic'],
            help="""which summary statistic to choose for the parameters
            in PARAMETER. Available options are `point` (the point
            estimate), `stderr` (standard deviation of the point
            estimator), or `tstatistic` (the t-test statistic testing
            whether the respected parameter is non zero). If you want to
            plot more than one parameter field, then separate these by a
            space.""")

########################################################################
# Arguments specific for the application of masks
########################################################################

    parser.add_argument('--round-to-integer',
            action='store_true',
            help="""if this flag is set and the output field is the
            intercept of a fit, then round the grey scales of the
            intercept field to the next nearest integer, and cast the
            array to integer data type.""")

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

from ...plot import picture

from ...nifti import image2nii

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

    layout_sdummy(df_layout, 'nii',
            template=args.nii,
            urname=args.population_space,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'file',
            template=args.fit,
            urname=args.population_space,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    ####################################################################
    # Respect the mask mask
    ####################################################################

    mask = True

    if args.ignore_mask:
        mask = False

    ####################################################################
    # Apply wrapper
    ####################################################################

    def wm(r):
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)

        wrapper(name                  = name,
                df                    = df,
                index                 = r.Index,
                filename              = r.file,
                template              = r.nii,
                verbose               = args.verbose,
                vb                    = args.population_space,
                params                = args.parameter,
                values                = args.value,
                mask                  = mask,
                round_to_integer      = args.round_to_integer)

    it =  df_layout.itertuples()

    for r in it:
        wm(r)

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

def wrapper(name, df, index, filename, template, verbose,
        vb, params, values, mask, round_to_integer):

    ####################################################################
    # Load fit from disk
    ####################################################################

    result = load_result(filename, name, df, index, vb, verbose)
    if df.ix[index,'valid'] == False:
        return

    if verbose > 1:
        print('{}: Description of the fit:'.format(result.name.name()))
        print(result.describe())

    result.mask(mask, verbose)

    for param in params:
        if param in result.parameter_dict.keys():
            for value in values:

                fname = template.format(param, value)
                dfile = os.path.dirname(fname)
                if dfile and not isdir(dfile):
                   os.makedirs(dfile)
                field = result.get_field(param=param, value=value)

                if value == 'intercept' and round_to_integer:
                    field = field.round()

                if verbose:
                    print('{}: Save ({}) to {}'.format(
                        name.name(), field.name, fname))
                ni.save(image2nii(field), fname)
