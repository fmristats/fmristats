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

Create a block stimulus instances from MATLAB coded logfiles

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
            """Setup of a two-block stimulus design""")

    specific.add_argument('--mat',
            default='../raw/mat/{2}/{0}-{1:04d}-{2}-{3}.mat',
            help = """A Matlab coded stimulus design""")

    ####################################################################
    # File handling
    ####################################################################

    file_handling = parser.add_argument_group(
            """File handling""")

    lock_handling = file_handling.add_mutually_exclusive_group()

    lock_handling.add_argument('-r', '--remove-lock',
            action='store_true',
            help="""Remove lock, if file is locked. This is useful, if
            used together with -s/--skip to remove orphan locks.""")

    lock_handling.add_argument('-i', '--ignore-lock',
            action='store_true',
            help="""Ignore lock, if file is locked. Together with
            -s/--skip this will also remove orphan locks.""")

    skip_force = file_handling.add_mutually_exclusive_group()

    skip_force.add_argument('-f', '--force',
            action='store_true',
            help="""Force re-writing any files""")

    skip_force.add_argument('-s', '--skip',
            action='store_true',
            help="""Do not perform any calculations.""")

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
            help="""Number of cores to use. Default is the number of
            cores on the machine.""")

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

import scipy.io

from ..api.fmristudy import get_study

from ...lock import Lock

from ...name import Identifier

from ...study import Study

from ...stimulus import Block

from ...matlab import mat2block

########################################################################

def call(args):

    study = get_study(args)

    if study is None:
        print('No study found. Nothing to do.')
        return

    ####################################################################
    # Options
    ####################################################################

    remove_lock       = args.remove_lock
    ignore_lock       = args.ignore_lock
    force             = args.force
    skip              = args.skip
    verbose           = args.verbose

    ####################################################################
    # Create the iterator
    ####################################################################

    study.update_layout({'mat':args.mat})

    study_iterator = study.iterate('stimulus',
            new=['stimulus', 'mat'],
            integer_index=True)

    df = study_iterator.df.copy()

    df['locked'] = False

    ####################################################################
    # Wrapper
    ####################################################################

    def wm(index, stimulus, file_stimulus, file_mat, name):

        if type(stimulus) is Lock:
            if remove_lock or ignore_lock:
                if verbose:
                    print('{}: Remove lock'.format(name.name()))
                stimulus.unlock()
                if remove_lock:
                    return
            else:
                if verbose:
                    print('{}: Locked'.format(name.name()))
                return

        if stimulus is not None and not force:
            if verbose:
                print('{}: Stimulus already exists. Use -f/--force to overwrite'.format(name.name()))
            return

        if skip:
            return

        if verbose:
            print('{}: Lock: {}'.format(name.name(), file_stimulus))

        lock = Lock(name, 'fmriblock', file_stimulus)
        df.ix[index, 'locked'] = True

        dfile = os.path.dirname(file_stimulus)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        lock.save(file_stimulus)

        ####################################################################
        # Load MATLAB instance from disk
        ####################################################################

        try:
            mat = scipy.io.loadmat(file_mat)
            if verbose:
                print('{}: Read: {}'.format(name.name(), file_mat))
        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to read: {}, {}'.format(name.name(),
                file_mat, e))

        if lock.conditional_unlock(df, index, verbose):
            df.ix[index,'valid'] = False
            return

        ####################################################################
        # Create stimulus instance
        ####################################################################

        try:
            stimulus = mat2block(mat, name=name)

            if verbose:
                print('{}: Save {}'.format(name.name(), file_stimulus))

            stimulus.save(file_stimulus)
            df.ix[index,'locked'] = False

        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to create: {}, {}'.format(name.name(),
                file_stimulus, e))
            lock.conditional_unlock(df, index, verbose, True)
            return

        return

    ####################################################################

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for index, name, files, instances in study_iterator:
                stimulus = instances['stimulus']
                file_stimulus = files['stimulus']
                file_mat = files['mat']
                pool.apply_async(wm, args=(index, stimulus,
                    file_stimulus, file_mat, name))

            pool.close()
            pool.join()
        except Exception as e:
            pool.close()
            pool.terminate()
            print('Pool execution has been terminated')
            print(e)
        finally:
            files = df.ix[df.locked, 'stimulus'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)
    else:
        try:
            print('Process protocol entries sequentially')
            for index, name, files, instances in study_iterator:
                stimulus = instances['stimulus']
                file_stimulus = files['stimulus']
                file_mat = files['mat']
                wm(index, stimulus, file_stimulus, file_mat, name)
        finally:
            files = df.ix[df.locked, 'stimulus'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)

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
