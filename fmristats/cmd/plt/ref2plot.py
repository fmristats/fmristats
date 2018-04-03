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

Quality assessment statistics for the ability to track head movements.

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

    parser.add_argument('--reference-maps',
            default='../data/ref/{2}/{0}-{1:04d}-{2}-{3}.ref',
            help='input file;' + hp.reference_maps)

########################################################################
# Output arguments
########################################################################

    parser.add_argument('-d', '--directory',
            default='figures/{2}/assessments',
            help="""output directory""")

    parser.add_argument('-e', '--extension',
            default='pdf',
            help="""anything that matplotlib understands""")

    parser.add_argument('--dpi',
            type=int,
            default=1200,
            help="""dpi.""")

    parser.add_argument('-t', '--template',
            default='{0}-{1:04d}-{2}-{3}',
            help="""template for the file name""")

#    parser.add_argument('-o', '--protocol-log',
#            default='logs/{}-ref2plot.pkl',
#            help=hp.protocol_log)

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

########################################################################
# Arguments specific for plot
########################################################################

    parser.add_argument('--grubbs',
            type=float,
            help=hp.grubbs)

########################################################################
# Miscellaneous
########################################################################

    #parser.add_argument('-f', '--force',
    #        action='store_true',
    #        help=hp.force.format('result'))

    parser.add_argument('-v', '--verbose',
            action='store_true',
            help=hp.verbose)

########################################################################
# Multiprocessing
########################################################################

    #parser.add_argument('-j', '--cores',
    #        type=int,
    #        help=hp.cores)

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

from ...load import load_refmaps

from ...name import Identifier

from ...protocol import layout_dummy

from ...session import Session

from ...reference import ReferenceMaps

from ...euler import rotation_matrix_to_euler_angles

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

#from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

import scipy.stats.distributions as dist

import matplotlib.pyplot as pt

pt.ioff()

########################################################################

def call(args):
#    output = args.protocol_log.format(
#            datetime.datetime.now().strftime('%Y-%m-%d-%H%M'))

    if args.strftime == 'short':
        args.strftime = '%Y-%m-%d'

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args, fall_back=args.reference_maps)

    if df is None:
        sys.exit()

    ####################################################################
    # Add file layout
    ####################################################################

    df_layout = df.copy()

    layout_dummy(df_layout, 'file_semiaxis',
            template=args.template + '-semi-axis.' + args.extension,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'file_euler',
            template=args.template + '-euler-angles.' + args.extension,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'file_euler_studenised',
            template=args.template + '-euler-angles-studenised.' + args.extension,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'file_bary',
            template=args.template + '-bary-centre.' + args.extension,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'file_bary_studenised',
            template=args.template + '-bary-centre-studenised.' + args.extension,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'dir',
            template=args.directory,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'file',
            template=args.reference_maps,
            strftime=args.strftime
            )

    ####################################################################
    # Apply wrapper
    ####################################################################

    def wm(r):
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)

        try:
            dfile = r.dir
            if dfile and not isdir(dfile):
                os.makedirs(dfile)
        except Exception as e:
            print('{}: {}'.format(name.name(), e))

        wrapper(
                name                  = name,
                df                    = df,
                index                 = r.Index,
                filename              = r.file,
                directory             = r.dir,

                file_semiaxis         = r.file_semiaxis,
                file_euler            = r.file_euler,
                file_euler_studenised = r.file_euler_studenised,
                file_bary             = r.file_bary,
                file_bary_studenised  = r.file_bary_studenised,

                sgnf                  = args.grubbs,
                dpi                   = args.dpi,
                verbose               = args.verbose)

    it =  df_layout.itertuples()

    for r in it:
        wm(r)

#######################################################################
# Plot quality measures of head movement
#######################################################################

def qc_plot(x, o, t, outlying_cycles, ax, studenise=False):
    """
    Parameters
    ----------
    x : ndarray, shape (3,n)
        Values to plot
    o : ndarray, shape (3,n)
        Outlying per axis
    t : ndarray, shape (n,)
        Time vector for the x-axis
    outlying_cycles : ndarray, shape (n,), dtype: bool
        Outlying scans
    ax : Axes
        Axes to which to plot
    studenise : bool
        Whether to studenise x first
    """
    #palette = sb.color_palette("muted", n_colors=3)
    if studenise:
        x = (x.T - x[...,~outlying_cycles].mean(axis=1)) / x[...,~outlying_cycles].std(axis=1)
        x = x.T
    ax.plot(t[~o[0]], x[0,~o[0]], '-', label='1st coordinate')#, c=palette[0])
    ax.plot(t[~o[1]], x[1,~o[1]], '-', label='2nd coordinate')#, c=palette[1])
    ax.plot(t[~o[2]], x[2,~o[2]], '-', label='3rd coordinate')#, c=palette[2])
    if outlying_cycles.any():
        ax.plot(t[o[0]], x[0,o[0]], '*', label='1st outlying')#, c=palette[0])
        ax.plot(t[o[1]], x[1,o[1]], '*', label='2nd outlying')#, c=palette[1])
        ax.plot(t[o[2]], x[2,o[2]], '*', label='3rd outlying')#, c=palette[2])
    if studenise:
        ax.axhline(0, c='k', ls='-', lw=0.5)
        ax.axhline(dist.norm.ppf(0.025), c='g', ls='--', lw=0.5)
        ax.axhline(dist.norm.ppf(1-0.025), c='g', ls='--', lw=0.5)
        ax.axhline(dist.norm.ppf(0.01), c='r', ls='-.', lw=0.5)
        ax.axhline(dist.norm.ppf(1-0.01), c='r', ls='-.', lw=0.5)

########################################################################

def wrapper(name, df, index, filename, directory, file_semiaxis,
        file_euler, file_euler_studenised,
        file_bary, file_bary_studenised, sgnf, dpi,
        verbose):

    ####################################################################
    # Load reference maps from disk
    ####################################################################

    reference_maps = load_refmaps(filename, name, df, index, verbose)
    if df.ix[index,'valid'] == False:
        return

    if (sgnf is not None) and (not np.isclose(sgnf, 1)):
        if verbose:
            print('{}: Detect outlying scans'.format(name.name()))
        reference_maps.detect_outlying_scans(sgnf)

    ####################################################################
    # Get slice times and outlying
    ####################################################################

    slice_times     = reference_maps.slice_time[:,0]
    outlying        = reference_maps.outlying
    values          = reference_maps.values
    outlying_cycles = reference_maps.outlying_cycles

    #######################################################################
    # Studenised Semi axis norms
    #######################################################################

    fig, ax = pt.subplots()
    qc_plot(
            values[:3],
            outlying[:3],
            slice_times,
            outlying_cycles,
            ax, True)
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Studenised semi axis length')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
              ncol=2, fancybox=True, shadow=True)
    ax.set_ylim((-6,6))

    figurename = join(directory, file_semiaxis)
    if verbose:
        print('{}: Save figure to: {}'.format(name.name(), figurename))
    pt.savefig(figurename, dpi=dpi)
    pt.close()

    #######################################################################
    # Calculate Euler angles for all head movements
    #######################################################################

    fig, ax = pt.subplots()
    qc_plot(
            np.degrees(values[6:]),
            outlying[6:],
            slice_times,
            outlying_cycles,
            ax, False)
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Euler angels (degrees)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
              ncol=2, fancybox=True, shadow=True)
    ax.axhline(0, c='k', ls='--', lw=0.5)

    figurename = join(directory, file_euler)
    if verbose:
        print('{}: Save figure to: {}'.format(name.name(), figurename))
    pt.savefig(figurename, dpi=dpi)
    pt.close()

    fig, ax = pt.subplots()
    qc_plot(
            values[6:],
            outlying[6:],
            slice_times,
            outlying_cycles,
            ax, True)
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Euler angels (studenised)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
              ncol=2, fancybox=True, shadow=True)
    ax.axhline(0, c='k', ls='--', lw=0.5)
    ax.set_ylim((-6,6))

    figurename = join(directory, file_euler_studenised)
    if verbose:
        print('{}: Save figure to: {}'.format(name.name(), figurename))
    pt.savefig(figurename, dpi=dpi)
    pt.close()

    #######################################################################
    # Plot barycentre of the head with respect to different references
    #######################################################################

    fig, ax = pt.subplots()
    qc_plot(
            values[3:6],
            outlying[3:6],
            slice_times,
            outlying_cycles,
            ax, False)
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Coordinates of bary centres (mm)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
              ncol=2, fancybox=True, shadow=True)

    figurename = join(directory, file_bary)
    if verbose:
        print('{}: Save figure to: {}'.format(name.name(), figurename))
    pt.savefig(figurename, dpi=dpi)
    pt.close()

    fig, ax = pt.subplots()
    qc_plot(
            values[3:6],
            outlying[3:6],
            slice_times,
            outlying_cycles,
            ax, True)
    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Coordinates of bary centres (studenised)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
              ncol=2, fancybox=True, shadow=True)
    ax.set_ylim((-6,6))

    figurename = join(directory, file_bary_studenised)
    if verbose:
        print('{}: Save figure to: {}'.format(name.name(), figurename))
    pt.savefig(figurename, dpi=dpi)
    pt.close()
