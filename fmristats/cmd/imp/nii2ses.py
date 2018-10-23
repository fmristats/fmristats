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

Command line tool to create a fmristats session instance

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
# Input arguments
########################################################################

    parser.add_argument('--nii',
            default='../raw/nii/{2}/{0}-{1:04d}-{2}-{3}.nii',
            help='input file;' + hp.nii)

########################################################################
# Arguments specific for the setup of a session instance
########################################################################

    parser.add_argument('--epi-code',
            type=int,
            help=hp.epi_code)

    foreground_handling = parser.add_mutually_exclusive_group()

    foreground_handling.add_argument('--detect-foreground',
            action='store_true',
            help=hp.detect_foreground)

    foreground_handling.add_argument('--set-foreground',
            action='store_true',
            help=hp.set_foreground)

    parser.add_argument('--foreground',
            default='../raw/foreground/{2}/{0}-{1:04d}-{2}-{3}-foreground.nii.gz',
            help=hp.set_foreground)

########################################################################
# Miscellaneous
########################################################################

    lock_handling = parser.add_mutually_exclusive_group()

    lock_handling.add_argument('-r', '--remove-lock',
            action='store_true',
            help=hp.remove_lock.format('session'))

    lock_handling.add_argument('-i', '--ignore-lock',
            action='store_true',
            help=hp.ignore_lock.format('session'))

    file_handling = parser.add_mutually_exclusive_group()

    file_handling.add_argument('-f', '--force',
            action='store_true',
            help=hp.force.format('session'))

    file_handling.add_argument('-s', '--skip',
            action='store_true',
            help=hp.skip.format('session'))

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

from ...load import load_block_irritation, load_session

from ...name import Identifier

from ...study import Study

from ...irritation import Block

from ...session import Session, fmrisetup

from ...nifti import nii2session

import nibabel as ni

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

########################################################################

def call(args):
    output = args.protocol_log.format(
            datetime.datetime.now().strftime('%Y-%m-%d-%H%M'))

    if args.strftime == 'short':
        args.strftime = '%Y-%m-%d'

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args, fall_back=args.irritation)

    if df is None:
        sys.exit()

    if not 'epi' in df.columns:
        print("""
No epi column found. You need to provide an EPI code via --epi-code,
it must be integer, within [-3,3], and not null.""")

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
        'scale_type':args.scale_type,
        'foreground':args.foreground,
        'nii':args.nii,
        }

    study = Study(df, df, layout=layout, strftime=args.strftime)

    study_iterator = study.iterate('irritation', 'session',
            new=['session', 'nii', 'foreground'],
            vb_name=args.vb_name,
            diffeomorphism_name=args.diffeomorphism_name,
            scale_type=args.scale_type)

    ####################################################################
    # Apply wrapper
    ####################################################################

    df['locked'] = False

    def wm(r):
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)

        try:
            dfile = os.path.dirname(r.file)
            if dfile and not isdir(dfile):
                os.makedirs(dfile)
        except Exception as e:
            print('{}: {}'.format(name.name(), e))

        wrapper(name              = name,
                df                = df,
                index             = r.Index,
                remove_lock       = args.remove_lock,
                ignore_lock       = args.ignore_lock,
                force             = args.force,
                skip              = args.skip,
                verbose           = args.verbose,
                file              = r.file,

                file_nii = r.nii,
                file_irr = r.irr,

                epi_code = r.epi,
                detect_foreground = args.detect_foreground,
                set_foreground = args.set_foreground,
                foreground = r.foreground,
                )

    it =  df_layout.itertuples()

    if len(df_layout) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            results = pool.map(wm, it)
            pool.close()
            pool.join()
        except Exception as e:
            pool.close()
            pool.terminate()
            print('Pool execution has been terminated')
            print(e)
        finally:
            files = df_layout.ix[df.locked, 'file'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)
            del df['locked']
    else:
        try:
            print('Process protocol entries sequentially')
            for r in it:
                wm(r)
        finally:
            files = df_layout.ix[df.locked, 'file'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)
            del df['locked']

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

def wrapper(name, df, index, remove_lock, ignore_lock, force, skip,
        verbose, file, file_nii , file_irr, epi_code,
        detect_foreground, set_foreground, foreground):

    if isfile(file):
        instance = load_session(file, name, df, index, verbose)
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
        else:
            if df.ix[index,'valid'] and not force:
                if verbose:
                    print('{}: Valid'.format(name.name()))
                return
            else:
                if skip:
                    if verbose:
                        print('{}: Invalid'.format(name.name()))
                    return

    if verbose:
        print('{}: Lock: {}'.format(name.name(), file))

    lock = Lock(name, 'fmrifit', file)
    df.ix[index, 'locked'] = True
    lock.save(file)
    df.ix[index,'valid'] = True

    ####################################################################
    # Load irritation instance from disk
    ####################################################################

    irritation = load_block_irritation(file_irr, name, df, index, verbose)
    if lock.conditional_unlock(df, index, verbose):
        return

    ################################################################
    # Load image data
    ################################################################

    try:
        img = ni.load(file_nii)
        if verbose:
            print('{}: Read: {}'.format(name.name(), file_nii))
    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to read: {}'.format(name.name(), file_nii))
        print('{}: Exception: {}'.format(name.name(), e))

    if lock.conditional_unlock(df, index, verbose):
        return

    ####################################################################
    # Create session instance
    ####################################################################

    try:
        session = nii2session(
                name=name,
                nii=img,
                epi_code=epi_code)

        fmrisetup(session = session, irritation = irritation)

        if detect_foreground:
            if verbose:
                print('{}: Detect foreground'.format(name.name()))
            session.fit_foreground()
        elif set_foreground:
            foreground = ni.load(foreground)
            if verbose:
                print('{}: Set foreground: {}'.format(name.name(), foreground))
            session.set_foreground(foreground.data)

        if verbose:
            print('{}: Save: {}'.format(name.name(), file))

        session.save(file)
        df.ix[index,'locked'] = False

    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to create: {}'.format(name.name(), file))
        print('{}: Exception: {}'.format(name.name(), e))
        lock.conditional_unlock(df, index, verbose, True)
        return

    return
