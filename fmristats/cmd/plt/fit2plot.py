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

Plot summary field statistics of the effect field of a fitted Signal
Model, e.g. the estimated activation pattern fields, activation
certainty fields, t-test statistic fields, or intercept fields.

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

    parser.add_argument('--figure',
            default='../results/figures/{2}/fields/{0}-{1:04d}-{2}-{3}-{4}-{5}-{{}}-{{}}.pdf',
            help="""output file""")

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-fit2plot.pkl',
            help=hp.protocol_log)

########################################################################
# Arguments specific for the RSM Signal Model: where to
########################################################################

    where_to_fit = parser.add_argument_group(
            """where to fit""",
            """By default, fmrifit will respect the brain mask which is
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

    parser.add_argument('--protocol',
            help=hp.protocol)

    parser.add_argument('--cohort',
            help=hp.cohort)

    parser.add_argument('--id',
            type=int,
            nargs='+',
            help=hp.j)

    # TODO: test if I need this when not using the protocol interface
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
# Arguments specific for creating pictures
########################################################################

    parser.add_argument('--parameter',
            type=str,
            nargs='+',
            default=['intercept', 'task'],
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
            estimator), or `tstatistic` (t-test statistic testing
            whether the respected parameter is non zero). If you want to
            plot more than one parameter field, then separate these by a
            space.""")

    parser.add_argument('--epi',
            type=int,
            nargs='+',
            default=[3],
            #default=[1,2,3],
            help="""plot along EPI-CODE axis. You may provide more than
            one. Default is: 1 2 3.""")

    parser.add_argument('--nx',
            type=int, default=6,
            help="""number of images on the x-axis.""")

    parser.add_argument('--ny',
            type=int, default=6,
            help="""number of images on the y-axis.""")

    parser.add_argument('--interpolation',
            default='none',
            help="""anything understood by pyplot.imshow.""")

    parser.add_argument('--cmap',
            type=str,
            nargs='+',
            default=['magma', 'inferno', 'seismic'],
            help="""anything that is understood by matplotlib. For a
            list see `http://matplotlib.org/users/colormaps.html`.
            Colour maps are allocated to parameter-value fields in the
            same order as they appear in VALUE.  E.g., if VALUE is the
            (default) list: point std teststat, then the first entry is
            used for the point estimate, second for standard deviation
            field, third for test fields.""")

    parser.add_argument('--dpi',
            type=int,
            default=1200,
            help="""dpi.""")

########################################################################
# Arguments specific
########################################################################

    parser.add_argument('-m', '--mark-peak',
            action='store_true',
            help=hp.verbose)

########################################################################
# Miscellaneous
########################################################################

    misc = parser.add_argument_group(
            """miscellaneous""")

    misc.add_argument('-f', '--force',
            action='store_true',
            help=hp.force.format('result'))

    misc.add_argument('-v', '--verbose',
            action='count',
            default=0,
            help=hp.verbose)

#    misc.add_argument('-j', '--cores',
#            type=int,
#            help=hp.cores)

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

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

#from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

import matplotlib.pyplot as pt

import matplotlib.cm as cm

pt.ioff()

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

    layout_sdummy(df_layout, 'template',
            template=args.figure,
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

    cmaps = {v:c for (v,c) in zip(args.value, args.cmap)}

    def wm(r):
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)

        wrapper(name                  = name,
                df                    = df,
                index                 = r.Index,
                filename              = r.file,
                template              = r.template,
                verbose               = args.verbose,
                vb                    = args.population_space,

                params                = args.parameter,
                values                = args.value,

                epi_code              = args.epi,
                nx                    = args.nx,
                ny                    = args.ny,
                cmaps                 = cmaps,
                dpi                   = args.dpi,

                mask                  = mask,
                mark_peak             = args.mark_peak,
                )

    # TODO: it = study.iter('results')
    it =  df_layout.itertuples()

    for r in it:
        wm(r)

    #if len(df_layout) > 1 and ((args.cores is None) or (args.cores > 1)):
    #    try:
    #        pool = ThreadPool(args.cores)
    #        results = pool.map(wm, it)
    #        pool.close()
    #        pool.join()
    #    except Exception as e:
    #        pool.close()
    #        pool.terminate()
    #        print('Pool execution has been terminated')
    #        print(e)
    #    finally:
    #        pass
    #else:
    #    try:
    #        print('Process protocol entries sequentially')
    #        for r in it:
    #            wm(r)
    #    finally:
    #        pass

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
        vb, params, values, epi_code, nx, ny, cmaps,
        dpi, mask, mark_peak):

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
                if param == 'intercept':
                    cmap = cm.magma
                elif value == 'stderr':
                    cmap = cm.viridis
                else:
                    cmap = cm.bwr

                for epi in epi_code:

                    fname = template.format(param, value, epi)
                    dfile = os.path.dirname(fname)
                    if dfile and not isdir(dfile):
                       os.makedirs(dfile)
                    field = result.get_field(param, value)

                    if verbose:
                        print('{}: Save ({}) to {}'.format(
                            name.name(), field.name, fname))

                    pt.figure(figsize=(10,10))
                    picture(field,epi,nx,ny,cmap=cmap,
                            mark_peak=mark_peak)
                    pt.savefig(fname, dpi=dpi)
                    pt.close()
