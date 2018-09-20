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

    parser.add_argument('--vb-type',
            default='reference',
            choices=['reference', 'scanner', 'scan', 'fit'],
            #choices=['reference', 'scan', 'scanner'],
            help=hp.vb_name)

########################################################################
# Input arguments when provided with reference maps
########################################################################

    parser.add_argument('--cycle',
            type=int,
            help="""cycle to pick as reference. Needed if --vb-name is
            set to scan""")

########################################################################
# Configuration
########################################################################

    parser.add_argument('--resolution',
            default=2.,
            type=float,
            help="""(optional) applicable when vb is reference.""")

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

from ...study import Study

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

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args, fall_back=[args.fit, args.session])

    if df is None:
        sys.exit()

    ####################################################################
    # Add file layout
    ####################################################################

    if args.diffeomorphism_name is None:
        args.diffeomorphism_name = args.vb_type

    layout = {
        'irritation':args.irritation,
        'session':args.session,
        'reference_maps':args.reference_maps,
        'result':args.fit,
        'population_map':args.population_map,
        'diffeomorphism_name':args.diffeomorphism_name,
        'scale_type':args.scale_type}

    study = Study(df, df, layout=layout, strftime=args.strftime)

    study_iterator = study.iterate(
            'session',
            'reference_maps',
            'result',
            'population_map', new=['population_map'])

    ####################################################################
    # Apply wrapper
    ####################################################################

    df = study_iterator.df.copy()

    df['locked'] = False

    def wm(instance, session, reference_maps, result, filename, name):
        index = (name.cohort, name.j, name.paradigm, name.datetime)

        remove_lock       = args.remove_lock
        ignore_lock       = args.ignore_lock
        force             = args.force
        skip              = args.skip
        verbose           = args.verbose

        vb_type           = args.vb_type
        resolution        = args.resolution
        scan_cycle        = args.cycle

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
                print('{}: PopulationMap already exists. Use -f/--force to overwrite'.format(name.name()))
            return

        if skip:
            return

        if verbose:
            print('{}: Lock: {}'.format(name.name(), filename))

        lock = Lock(name, 'fmripop', filename)
        df.ix[index, 'locked'] = True

        dfile = os.path.dirname(filename)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        lock.save(filename)

        ####################################################################
        # Create population map instance from a session instance
        ####################################################################

        if (vb_type == 'scanner') or (vb_type == 'reference'):

            if session is None:
                print('{}: No session provided'.format(name.name()))
                df.ix[index,'valid'] = False
                lock.conditional_unlock(df, index, verbose)
                return

            if vb_type == 'scanner':
                population_map = pmap_scanner(session=session)
                if verbose:
                    print('{}: VB space equals scanner space (with native resolution)'.format(
                        name.name()))

            if vb_type == 'reference':
                population_map = pmap_reference(session=session, resolution=resolution)
                if verbose:
                    print('{}: VB space equals scanner space (with resolution ({} mm)**3)'.format(
                        name.name(), resolution))

        ####################################################################
        # Create population map instance from a session and reference
        # map instance
        ####################################################################

        if (vb_type == 'scan'):

            if (session is None) or (reference_maps is None) or (scan_cycle is None):
                print('{}: No session or reference maps provided or CYCLE not defined'.format(name.name()))
                df.ix[index,'valid'] = False
                lock.conditional_unlock(df, index, verbose)
                return

            population_map = pmap_scan(
                    reference_maps=reference_maps,
                    session=session,
                    scan_cycle=scan_cycle)

            if verbose:
                print('{}: VB space now equals subject position during scan cycle: {:d}'.format(
                    name.name(), scan_cycle))

        ####################################################################
        # Create population map instance from a result instance
        ####################################################################

        if (vb_type == 'fit'):

            if result is None:
                print('{}: No fit provided'.format(name.name()))
                df.ix[index,'valid'] = False
                lock.conditional_unlock(df, index, verbose)
                return

            population_map = result.population_map
            population_map.set_template(template=result.get_field('intercept', 'point'))

        try:
            if verbose:
                print('{}: Save: {}'.format(name.name(), filename))

            population_map.save(filename)
            df.ix[index,'locked'] = False

        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to create: {}'.format(name.name(), filename))
            print('{}: Exception: {}'.format(name.name(), e))
            lock.conditional_unlock(df, index, verbose, True)
            return

        return

    ####################################################################

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for name, files, instances in study_iterator:
                population_map = instances['population_map']
                session = instances['session']
                reference_maps = instances['reference_maps']
                result = instances['result']
                filename = files['population_map']
                pool.apply_async(wm, args=(population_map, session,
                    reference_maps, result, filename, name))

            pool.close()
            pool.join()
        except Exception as e:
            pool.close()
            pool.terminate()
            print('Pool execution has been terminated')
            print(e)
        finally:
            files = df.ix[df.locked, 'population_map'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)
    else:
        try:
            print('Process protocol entries sequentially')
            for name, files, instances in study_iterator:
                population_map = instances['population_map']
                session = instances['session']
                reference_maps = instances['reference_maps']
                result = instances['result']
                filename = files['population_map']
                wm(population_map, session, reference_maps, result, filename, name)
        finally:
            files = df.ix[df.locked, 'population_map'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)

    ####################################################################
    # Write study to disk
    ####################################################################

    # TODO: filter study protocol entries form invalid entries

    if args.out is not None:
        if args.verbose:
            print('Save: {}'.format(args.out))

        dfile = os.path.dirname(args.out)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        study.save(args.out)
