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

Prune statistic field from non-brain areas

"""

    # TODO: also give the option --reference-maps None or none,
    # which will assume no movement. This is useful for analysing
    # phantom data.

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
    # Handling cut-offs
    ####################################################################

    handling_df_cutoff = parser.add_mutually_exclusive_group()

    handling_df_cutoff.add_argument('-p', '--proportion',
            type=float,
            default=.6842,
            help="""estimates which degrees of freedom are below the
            proportional threshold of the degrees of freedom in the
            effect field estimate are set to null.""")

    handling_df_cutoff.add_argument('-t', '--threshold',
            type=int,
            help="""estimates which degrees of freedom are below the
            threshold of the degrees of freedom in the effect field
            estimate are set to null.""")

    ####################################################################
    # Miscellaneous
    ####################################################################

    parser.add_argument('-f', '--force-mask-overwrite',
            action='store_true',
            help="""Overwrite mask if it already exits""")

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

########################################################################

def call(args):

    study = get_study(args)

    if study is None:
        print('No study found. Nothing to do.')
        sys.exit()

    ####################################################################
    # Options
    ####################################################################

    force         = args.force_mask_overwrite
    verbose       = args.verbose
    threshold_df  = args.threshold
    proportion_df = args.proportion

    ####################################################################
    # Create the iterator
    ####################################################################

    study_iterator = study.iterate('result', new=['result'])
            #vb_name=args.vb_name,
            #diffeomorphism_name=args.diffeomorphism_name,
            #scale_type=args.scale_type)

    df = study_iterator.df.copy()

    ####################################################################
    # Wrapper
    ####################################################################

    def wm(result, filename, name):

        if type(result) is Lock:
            print('{}: Result file is locked. Still fitting?'.format(
                name.name()))
            return

        if verbose > 2:
            print("""{}: Description of the fit:
                {}
                {}""".format(name.name(), result.describe(),
                    result.population_map.describe()))

        if hasattr(result.population_map, 'vb_mask') and not force:
            print("""{}:
            VB mask already present, force overwrite with --force""".format(
                name.name()))
            return

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

    ####################################################################
    # Write study to disk
    ####################################################################

    if args.out is not None:
        if args.verbose:
            print('Save: {}'.format(args.out))

        dfile = os.path.dirname(args.out)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        study.save(args.out)
