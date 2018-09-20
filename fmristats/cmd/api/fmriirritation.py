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

Create a block irritation instances

"""

########################################################################
#
# Command line program
#
########################################################################

import fmristats.cmd.hp as hp

import argparse

from ...study import add_study_parser

def create_argument_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=hp.epilog)

    add_study_parser(parser)

    parser.add_argument('-o', '--out',
            help="""Save (new) study instance to OUT""")

########################################################################
# Arguments specific for the setup of an irritation instance
########################################################################

    parser.add_argument('--onsetsx',
            type=float,
            nargs='+',
            help=hp.onsetx)

    parser.add_argument('--onsetsy',
            type=float,
            nargs='+',
            help=hp.onsety)

    parser.add_argument('--durationsx',
            type=float,
            help=hp.durationx)

    parser.add_argument('--durationsy',
            type=float,
            help=hp.durationy)

    parser.add_argument('--namex',
            default='control',
            help="""name of block x""")

    parser.add_argument('--namey',
            default='stimulus',
            help="""name of block y""")

########################################################################
# If not using the protocol API, you need this
########################################################################

    parser.add_argument('--epi-code',
            type=int,
            help=hp.epi_code)

########################################################################
# Miscellaneous
########################################################################

    lock_handling = parser.add_mutually_exclusive_group()

    lock_handling.add_argument('-r', '--remove-lock',
            action='store_true',
            help=hp.remove_lock.format('irritation'))

    lock_handling.add_argument('-i', '--ignore-lock',
            action='store_true',
            help=hp.ignore_lock.format('irritation'))

    file_handling = parser.add_mutually_exclusive_group()

    file_handling.add_argument('-f', '--force',
            action='store_true',
            help=hp.force.format('irritation'))

    file_handling.add_argument('-s', '--skip',
            action='store_true',
            help=hp.skip.format('irritation'))

    parser.add_argument('-v', '--verbose',
            action='count',
            default=0,
            help=hp.verbose)

########################################################################
# Multiprocessing
########################################################################

    parser.add_argument('-j', '--cores',
            type=int,
            help=hp.cores)

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

from ...load import load_block_irritation

from ...name import Identifier

from ...study import Study

#from ...protocol import layout_dummy

from ...irritation import Block

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

########################################################################

def call(args):

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args)

    if df is None:
        sys.exit()

    ####################################################################
    # Add file layout
    ####################################################################

    layout = {
        'irritation':args.irritation,
        'session':args.session,
        'reference_maps':args.reference_maps,
        'result':args.fit,
        'population_map':args.population_map,
        'diffeomorphism_name':args.diffeomorphism_name,
        'scale_type':args.scale_type}

    study = Study(df, df, layout=layout, strftime=args.strftime)

    study_iterator = study.iterate('irritation', new=['irritation'])

    ####################################################################
    # Apply wrapper
    ####################################################################

    df = study_iterator.df.copy()

    df['locked'] = False

    def wm(instance, filename, name):
        index = (name.cohort, name.j, name.paradigm, name.datetime)

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
                print('{}: Irritation already exists. Use -f/--force to overwrite'.format(name.name()))
            return

        if skip:
            return

        if verbose:
            print('{}: Lock: {}'.format(name.name(), filename))

        lock = Lock(name, 'fmriirritation', filename)
        df.ix[index, 'locked'] = True

        dfile = os.path.dirname(filename)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        lock.save(filename)

        ####################################################################
        # Create irritation instance
        ####################################################################

        try:
            irritation = Block(name=name,
                    names=[namex, namey],
                    onsets={namex:onsetsx,namey:onsetsy},
                    durations={namex:durationsx,namey:durationsy})

            if verbose:
                print('{}: Save: {}'.format(name.name(), filename))

            irritation.save(filename)
            df.ix[index,'locked'] = False

        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to create: {}'.format(name.name(), filename))
            print('{}: Exception: {}'.format(name.name(), e))
            lock.conditional_unlock(df, index, verbose, True)
            return

    ####################################################################

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for name, files, instances in study_iterator:
                irritation = instances['irritation']
                filename = files['irritation']
                pool.apply_async(wm, args=(irritation, filename, name))

            pool.close()
            pool.join()
        except Exception as e:
            pool.close()
            pool.terminate()
            print('Pool execution has been terminated')
            print(e)
        finally:
            files = df.ix[df.locked, 'irritation'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)
    else:
        try:
            print('Process protocol entries sequentially')
            for name, files, instances in study_iterator:
                irritation = instances['irritation']
                filename = files['irritation']
                wm(irritation, filename, name)
        finally:
            files = df.ix[df.locked, 'irritation'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)

    ####################################################################
    # Write study to disk
    ####################################################################

    # TODO: filter study protocol entries form invalid entries

    if args.out is not None:
        if args.epi_code is None:
            print('Warning: study protocol has not been equipped with a valid EPI code')

        if args.verbose:
            print('Save: {}'.format(args.out))

        dfile = os.path.dirname(args.out)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        study.save(args.out)
