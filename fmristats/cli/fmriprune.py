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

########################################################################
#
# Command line program
#
########################################################################

from ..epilog import epilog

import argparse

def add_arguments(parser):

    ####################################################################
    # Handling cut-offs
    ####################################################################

    handling_df_cutoff = parser.add_mutually_exclusive_group()

    handling_df_cutoff.add_argument('--fraction',
            type=float,
            default=.6842,
            help="""estimates which degrees of freedom are below the
            proportional threshold of the degrees of freedom in the
            effect field estimate are set to null.""")

    handling_df_cutoff.add_argument('--threshold',
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
    # Push
    ####################################################################

    control_verbosity  = parser.add_argument_group(
        """Control whether to save the modified (thus overwrite the
        existing) study instance.""")

    control_verbosity.add_argument('-p', '--push',
        action='store_true',
        help="""Will save the modified (and thus overwrite the existing)
        study instance.""")

    ####################################################################
    # Multiprocessing
    ####################################################################

    control_multiprocessing  = parser.add_argument_group(
        """Multiprocessing""")

    control_multiprocessing.add_argument('-j', '--cores',
        type=int,
        default=1,
        help="""Number of threads to use. The implementation will
        usually try to run as many calculations and loops as possible in
        parallel -- this may suggest that it may be adventurous to
        process all entries in the study protocol sequentially (and this
        is the default). It is possible, however, to generate a thread
        for each protocol entry. Note that this may generate a lot of
        I/O-operations. If you set CORES to 0, then the number of cores
        on the machine will be used.""")

    return parser

def define_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=epilog)
    add_arguments(parser)
    return parser

from .fmristudy import add_study_arguments

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

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

from .fmristudy import get_study

from ..lock import Lock

from ..study import Study

########################################################################

def call(args):

    study = get_study(args)

    if study is None:
        print('Nothing to do.')
        return

    ####################################################################
    # Options
    ####################################################################

    force         = args.force_mask_overwrite
    verbose       = args.verbose
    threshold_df  = args.threshold
    proportion_df = args.fraction

    ####################################################################
    # Create the iterator
    ####################################################################

    study_iterator = study.iterate('result', new=['result'])

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
            VB mask already present, use
            -f/--force-mask-overwrite
            to overwrite existing mask""".format(
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

        try:
            if verbose:
                print('{}: Save: {}'.format(name.name(), filename))
            result.save(filename)
        except Exception as e:
            print('{}: Unable to write: {}, {}'.format(name.name(),
                filename, e))

    ###################################################################

    if args.cores == 0:
        args.cores = None

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for name, files, instances in study_iterator:
                result = instances['result']
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
                result = instances['result']
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
