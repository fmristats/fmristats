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

Quality assessment statistics for evaluating the ability to track the head
movements of a subject within an FMRI

"""

########################################################################
#
# Command line program
#
########################################################################

from ...epilog import epilog

import argparse

def define_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=epilog)

    ####################################################################
    # Specific arguments
    ####################################################################

    specific = parser.add_argument_group(
        """Arguments controlling the plots""")


    specific.add_argument('-d', '--directory',
        default='figures/head-movements',
        help="""Directory where to save the plots.""")

    specific.add_argument('--extension',
        default='pdf',
        help="""Format of the plots. May be anything that matplotlib
        understands.""")

    specific.add_argument('--dpi',
        type=int,
        default=1200,
        help="""The dpi to use (whenever relevant).""")

    specific.add_argument('--grubbs',
        type=float,
        help="""An outlier detection is performed to identify scans
        which may have been acquired during severe head movements. More
        precisely, a Grubbs' outlying test will be performed on the set
        of estimated principle semi axis for each full scan cycle on the
        given level of significance. When using fmririgid to create
        ReferenceMaps, the default is 0.1, and the information of
        outlying scans is saved to disk together with the estimated
        rigid body transformations. Then, when running fmrifit, this
        information is used. When setting --grubbs in fmrifit, outlier
        estimation is performed again.""")

    ####################################################################
    # Verbosity
    ####################################################################

    control_verbosity  = parser.add_argument_group(
        """Control the level of verbosity""")

    control_verbosity.add_argument('-v', '--verbose',
        action='count',
        default=0,
        help="""Increase output verbosity""")

    ####################################################################
    # Multiprocessing
    ####################################################################

    control_multiprocessing  = parser.add_argument_group(
        """Multiprocessing""")

    control_multiprocessing.add_argument('-j', '--cores',
        type=int,
        help="""Number of cores to use. Default is the number of cores
        on the machine.""")

    return parser

from ..api.fmristudy import add_study_arguments

def cmd():
    parser = define_parser()
    add_study_arguments(parser)
    args = parser.parse_args()
    call(args)

cmd.__doc__ = __doc__

########################################################################
#
# Load libraries
#
########################################################################

import sys

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

from ..api.fmristudy import get_study

from ...lock import Lock

from ...study import Study

from ...reference import ReferenceMaps

from ...euler import rotation_matrix_to_euler_angles

import scipy.stats.distributions as dist

import matplotlib.pyplot as pt

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

def call(args):

    study = get_study(args)

    if study is None:
        sys.exit()

    study_iterator = study.iterate('reference_maps')

    df = study_iterator.df.copy()

    verbose   = args.verbose
    grubbs    = args.grubbs
    dpi       = args.dpi
    directory = args.directory

    if directory and not isdir(directory):
       os.makedirs(directory)

    pt.ioff()

    def wm(name, reference_maps):

        if (grubbs is not None) and (not np.isclose(grubbs, 1)):
            if verbose:
                print('{}: Detect outlying scans'.format(name.name()))
            reference_maps.detect_outlying_scans(grubbs)

        file1 = join(directory, name.name() + '-radii.' + args.extension)
        file2 = join(directory, name.name() + '-euler.' + args.extension)
        file3 = join(directory, name.name() + '-euler-studenised.' + args.extension)
        file4 = join(directory, name.name() + '-bary.' + args.extension)
        file5 = join(directory, name.name() + '-bary-studenised.' + args.extension)

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

        if verbose:
            print('{}: Save figure to: {}'.format(name.name(), file1))
        pt.savefig(file1, dpi=dpi)
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

        if verbose:
            print('{}: Save figure to: {}'.format(name.name(), file2))
        pt.savefig(file2, dpi=dpi)
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

        if verbose:
            print('{}: Save figure to: {}'.format(name.name(), file3))
        pt.savefig(file3, dpi=dpi)
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

        if verbose:
            print('{}: Save figure to: {}'.format(name.name(), file4))
        pt.savefig(file4, dpi=dpi)
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

        if verbose:
            print('{}: Save figure to: {}'.format(name.name(), file5))
        pt.savefig(file5, dpi=dpi)
        pt.close()

    ####################################################################

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for name, files, instances in study_iterator:
                reference_maps  = instances['reference_maps']
                wm
                pool.apply_async(wm, args=(name, reference_maps))

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
            for name, instances in study_iterator:
                reference_maps = instances['reference_maps']
                wm(name, reference_maps)
        finally:
            pass

    ####################################################################
    # Write study to disk
    ####################################################################

    if args.out is not None:
        if args.epi_code is None:
            print('Warning: study protocol has not been equipped with a valid EPI code')

        if args.verbose:
            print('Save: {}'.format(args.out))

        dfile = os.path.dirname(args.out)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        study.save(args.out)
