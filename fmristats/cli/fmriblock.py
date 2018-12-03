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

Create a block stimulus instances

"""

########################################################################

onsetx = """Times of onsets of condition x, which will be interpreted
as the control stimulus of the subject.  Should, e.g., be given in the
form:

                --onsetx 2 84 166 248 330

to be parsed correctly.  If onsets are not given in order, they will be
sorted automatically.  Must be specified if --stimulus is not."""

onsety = """Times of onset of condition y, which will be interpreted as
the stimulus stimulus of the subject.  Should, e.g., be given in the
form:

                --onsety 37 119 201 283

to be parsed correctly.  If onsets are not given in order, they will be
sorted automatically.  Must be specified if --stimulus is not."""

durationx = """Duration of the stimulus of the subject under condition
x.  Must be specified if --stimulus is not."""

durationy = """Duration of the stimulus of the subject under condition
y.  Must be specified if --stimulus is not."""

########################################################################
#
# Command line program
#
########################################################################

from ..epilog import epilog

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

    specific.add_argument('--onsetsx',
            type=float,
            nargs='+',
            help=onsetx)

    specific.add_argument('--onsetsy',
            type=float,
            nargs='+',
            help=onsety)

    specific.add_argument('--durationsx',
            type=float,
            help=durationx)

    specific.add_argument('--durationsy',
            type=float,
            help=durationy)

    specific.add_argument('--namex',
            default='control',
            help="""name of block x""")

    specific.add_argument('--namey',
            default='stimulus',
            help="""name of block y""")

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

from ..name import Identifier

from ..study import Study

from ..stimulus import Block

########################################################################

def call(args):

    study = get_study(args)

    if study is None:
        print('Nothing to do.')
        return

    study_iterator = study.iterate('stimulus', new=['stimulus'],
            integer_index=True)

    df = study_iterator.df.copy()

    df['locked'] = False

    remove_lock       = args.remove_lock
    ignore_lock       = args.ignore_lock
    force             = args.force
    skip              = args.skip
    verbose           = args.verbose

    namex    = args.namex
    namey    = args.namey
    onsetsx  = np.asarray(args.onsetsx)
    onsetsy  = np.asarray(args.onsetsy)
    durationsx = np.asarray(args.durationsx)
    durationsy = np.asarray(args.durationsy)

    def wm(index, instance, filename, name):

        if type(instance) is Lock:
            if remove_lock or ignore_lock:
                if verbose:
                    print('{}: Remove lock'.format(name.name()))
                instance.unlock()
                if remove_lock:
                    return
            else:
                if verbose:
                    print('{}: Locked'.format(name.name()))
                return

        if instance is not None and not force:
            if verbose:
                print('{}: Stimulus already exists. Use -f/--force to overwrite'.format(name.name()))
            return

        if skip:
            return

        if verbose:
            print('{}: Lock: {}'.format(name.name(), filename))

        lock = Lock(name, 'fmriblock', filename)
        df.ix[index, 'locked'] = True

        dfile = os.path.dirname(filename)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        lock.save(filename)

        ####################################################################
        # Create stimulus instance
        ####################################################################

        try:
            stimulus = Block(name=name,
                    names=[namex, namey],
                    onsets={namex:onsetsx,namey:onsetsy},
                    durations={namex:durationsx,namey:durationsy})

            if verbose:
                print('{}: Save: {}'.format(name.name(), filename))

            stimulus.save(filename)
            df.ix[index,'locked'] = False

        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to create: {}'.format(name.name(), filename))
            print('{}: Exception: {}'.format(name.name(), e))
            lock.conditional_unlock(df, index, verbose, True)
            return

    ####################################################################

    if args.cores == 0:
        args.cores = None

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for index, name, files, instances in study_iterator:
                stimulus = instances['stimulus']
                filename = files['stimulus']
                pool.apply_async(wm, args=(index, stimulus, filename, name))

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
                filename = files['stimulus']
                wm(index, stimulus, filename, name)
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
        if args.verbose:
            print('Save: {}'.format(args.out))

        dfile = os.path.dirname(args.out)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        study.save(args.out)

    if args.push:
        if args.verbose:
            print('Save: {}'.format(args.study))
        study.save(args.study)
