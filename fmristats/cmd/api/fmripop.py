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

Fit a signal model to FMRI data

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

    parser.add_argument('--session',
            default='../data/ses/{2}/{0}-{1:04d}-{2}-{3}.ses',
            help="""only needed if no fit is provided; """ + hp.session)

########################################################################
# Additional input arguments when a previous fit can be provided
########################################################################

    parser.add_argument('--nb-fit',
            default='../data/fit/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}-{5}.fit',
            help="""input file; if warp coefficient files are provided
            you don't need this""" + hp.sfit)

    parser.add_argument('--nb-scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help="""only needed if part of the template for --nb-fit""" + hp.scale_type)

    parser.add_argument('--nb',
            default='reference',
            help="""name of the population space that was originally
            used for the fit"""  + hp.population_space)

    parser.add_argument('--ignore-fit',
            action='store_true',
            help="""ignore fit""")

########################################################################
# Output arguments
########################################################################

    parser.add_argument('--population-space',
            default='reference',
            help='output name;' + hp.population_space_name)

    parser.add_argument('--population-map',
            default='../data/pop/{2}/{4}/{0}-{1:04d}-{2}-{3}-{4}.pop',
            help='output file;' + hp.population_map)

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-fmripop.pkl',
            help=hp.protocol_log)

########################################################################
# Configuration
########################################################################

    parser.add_argument('--resolution',
            default=2.,
            type=float,
            help="""(optional) only applicable when population space is
            reference.""")

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
# Miscellaneous
########################################################################

    lock_handling = parser.add_mutually_exclusive_group()

    lock_handling.add_argument('-r', '--remove-lock',
            action='store_true',
            help=hp.remove_lock.format('population map'))

    lock_handling.add_argument('-i', '--ignore-lock',
            action='store_true',
            help=hp.ignore_lock.format('population map'))

    file_handling = parser.add_mutually_exclusive_group()

    file_handling.add_argument('-f', '--force',
            action='store_true',
            help=hp.force.format('population map'))

    file_handling.add_argument('-s', '--skip',
            action='store_true',
            help=hp.skip.format('population map'))

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

from ...load import load_result, load_session, load_population_map

from ...name import Identifier

from ...protocol import layout_dummy, layout_sdummy

from ...session import Session

from ...reference import ReferenceMaps

from ...smodel import SignalModel, Result

from ...pmap import PopulationMap, pmap_scanner, pmap_reference

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

########################################################################

def call(args):
    output = args.protocol_log.format(
            datetime.datetime.now().strftime('%Y-%m-%d-%H%M'))

    if args.strftime == 'short':
        args.strftime = '%Y-%m-%d'

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args, fall_back=[args.nb_fit, args.session])

    if df is None:
        sys.exit()

    ####################################################################
    # Add file layout
    ####################################################################

    df_layout = df.copy()

    layout_dummy(df_layout, 'ses',
            template=args.session,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'fit',
            template=args.nb_fit,
            urname=args.nb,
            scale_type=args.nb_scale_type,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'file',
            template=args.population_map,
            urname=args.population_space,
            scale_type=None,
            strftime=args.strftime
            )

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

                file_fit          = r.fit,
                file_ses          = r.ses,
                vb                = args.population_space,
                resolution        = args.resolution,
                ignore_fit        = args.ignore_fit,
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

#######################################################################

def wrapper(name, df, index, remove_lock, ignore_lock, force, skip,
        verbose, file, file_fit, file_ses, vb, resolution, ignore_fit):

    if isfile(file):
        instance = load_population_map(file, name, df, index, vb, verbose)
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

    if skip:
        return

    if verbose:
        print('{}: Lock: {}'.format(name.name(), file))

    lock = Lock(name, 'fmrifit', file)
    df.ix[index, 'locked'] = True
    lock.save(file)
    df.ix[index,'valid'] = True

    ####################################################################
    # Load fit from disk
    ####################################################################

    if isfile(file_fit) and not ignore_fit:
        result = load_result(file_fit, name, df, index, vb, verbose)
        if lock.conditional_unlock(df, index, verbose):
            return

        if verbose > 1:
            print('{}: Description of the fit:'.format(result.name.name()))
            print(result.describe())

        population_map = result.population_map
        population_map.set_template(template=result.get_field('intercept', 'point'))

    else:
        session = load_session(file_ses, name, df, index, verbose)
        if lock.conditional_unlock(df, index, verbose):
            return

        if vb == 'scanner':
            population_map = pmap_scanner(session=session)
        if vb == 'reference':
            population_map = pmap_reference(session=session, resolution=resolution)

    if verbose:
        print('{}: Save: {}'.format(name.name(), file))

    try:
        population_map.save(file)
        df.ix[index,'locked'] = False
    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to write: {}'.format(name.name(), file))
        print('{}: Exception: {}'.format(name.name(), e))
        lock.conditional_unlock(df, index, verbose, True)

    if verbose > 1:
        print('{}: Description'.format(name.name()))
        print(result.describe())
