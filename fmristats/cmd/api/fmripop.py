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

Define a subject-specific standard space for, well, a subject.

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

    parser.add_argument('--vb-type',
            default='reference',
            choices=['reference', 'scanner', 'scan', 'fit'],
            #choices=['reference', 'scan', 'scanner'],
            help=hp.vb_name)

########################################################################
# Input arguments when a provided with reference maps
########################################################################

    parser.add_argument('--session',
            default='../data/ses/{2}/{0}-{1:04d}-{2}-{3}.ses',
            help="""needed if --vb-name is set to reference or scanner.""")

########################################################################
# Input arguments when a provided with reference maps
########################################################################

    parser.add_argument('--reference-maps',
            default='../data/ref/{2}/{0}-{1:04d}-{2}-{3}.ref',
            help="""needed if --vb-name is set to scan.""")

    parser.add_argument('--cycle',
            type=int,
            help="""cycle to pick as reference. Needed if --vb-name is
            set to scan""")

########################################################################
# Input arguments when provided with a previous fit
########################################################################

    parser.add_argument('--fit',
            default='../data/fit/{2}/{4}/{5}/{6}/{0}-{1:04d}-{2}-{3}-{4}.fit',
            help="""needed if --vb-name is set to fit.""" + hp.sfit)

    parser.add_argument('--scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help="""only needed as part of the template for --fit.""" + hp.scale_type)

    parser.add_argument('--nb-name',
            default='self',
            help="""name of the population space that was originally
            used for the fit. Only needed as part of the template for
            --fit."""  + hp.nb_name)

########################################################################
# Output arguments
########################################################################

    parser.add_argument('--population-map',
            default='../data/pop/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}.pop',
            help=hp.population_map)

    parser.add_argument('--vb',
            default='self',
            help=hp.vb)

    parser.add_argument('--diffeomorphism',
            help="""Name of the diffeomorphism. Default is VB_TYPE.""")

########################################################################
# Log file
########################################################################

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

from ...load import load_result, load_session, load_refmaps, load_population_map

from ...name import Identifier

from ...protocol import layout_dummy, layout_sdummy, layout_fdummy

from ...session import Session

from ...reference import ReferenceMaps

from ...smodel import SignalModel, Result

from ...pmap import PopulationMap, pmap_scanner, pmap_reference, pmap_scan

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

    df = get_df(args, fall_back=[args.fit, args.session])

    if df is None:
        sys.exit()

    ####################################################################
    # Add file layout
    ####################################################################

    if args.diffeomorphism is None:
        args.diffeomorphism = args.vb_type

    df_layout = df.copy()

    layout_dummy(df_layout, 'ses',
            template=args.session,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'ref',
            template=args.reference_maps,
            strftime=args.strftime
            )

    layout_fdummy(df_layout, 'fit',
            template=args.fit,
            vb=args.nb_name,
            diffeo=args.diffeomorphism,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'file',
            template=args.population_map,
            urname=args.vb,
            scale_type=args.diffeomorphism,
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

                file_ses          = r.ses,
                file_ref          = r.ref,
                file_fit          = r.fit,

                vb_type           = args.vb_type,
                vb                = args.vb,
                resolution        = args.resolution,
                scan_cycle        = args.cycle,
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
        verbose, file, file_ses, file_ref, file_fit, vb_type, vb,
        scan_cycle, resolution):

    if isfile(file):
        instance = load_population_map(file, name, df, index, vb, vb_type, verbose)
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

    if (vb_type == 'scanner') or (vb_type == 'reference'):

        session = load_session(file_ses, name, df, index, verbose)
        if lock.conditional_unlock(df, index, verbose):
            return

        if vb_type == 'scanner':
            population_map = pmap_scanner(session=session)
            if verbose:
                print('{}: Population space equals to scanner space (with native resolution)'.format(
                    name.name()))

        if vb_type == 'reference':
            population_map = pmap_reference(session=session, resolution=resolution)
            if verbose:
                print('{}: Population space equals scanner space (with resolution ({} mm)**3)'.format(
                    name.name(), resolution))

    if (vb_type == 'scan'):

        session = load_session(file_ses, name, df, index, verbose)
        if lock.conditional_unlock(df, index, verbose):
            return

        reference_maps = load_refmaps(file_ref, name, df, index, verbose)
        if lock.conditional_unlock(df, index, verbose):
            return

        population_map = pmap_scan(
                reference_maps=reference_maps,
                session=session,
                scan_cycle=scan_cycle)
        if verbose:
            print('{}: Population space equals subject position during scan cycle: {:d}'.format(
                name.name(), scan_cycle))

    if (vb_type == 'fit'):

       result = load_result(file_fit, name, df, index, vb, verbose)
       if lock.conditional_unlock(df, index, verbose):
           return

       if verbose > 1:
           print('{}: Description of the fit:'.format(result.name.name()))
           print(result.describe())

       population_map = result.population_map
       population_map.set_template(template=result.get_field('intercept', 'point'))

    try:
        if verbose:
            print('{}: Save: {}'.format(name.name(), file))

        population_map.save(file)
        df.ix[index,'locked'] = False

    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to create: {}'.format(name.name(), file))
        print('{}: Exception: {}'.format(name.name(), e))
        lock.conditional_unlock(df, index, verbose, True)
        return

    return
