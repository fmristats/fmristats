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

Create a standard space that is isometric to the reference space of a
subject

"""

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
        """Creating a standard space isometric to the reference space.""")

    specific.add_argument('--diffeomorphism-nb',
        default='scanner',
        choices=['scanner', 'identity', 'fit'],
        help="""Image space of diffeomorphism.""")

    specific.add_argument('--new-diffeomorphism',
        default='scanner',
        help="""Name to use for the fitted diffeomorphisms.""")

    specific.add_argument('--resolution',
        default=2.,
        type=float,
        help="""The resolution of the template in standard space. If set
        to 0, the resolution of the scanner will be used.""")

    specific.add_argument('--cycle',
        type=int,
        nargs='+',
        help="""Use the data in scan cycle CYCLE as a template for
        the subject reference space. Enumeration starts at 0 (i.e. the
        first scan cycle is cycle 0). It will be checked whether the
        specified cycle has been marked as a potential outlier (i.e. a
        scan cycle that shows more head movements that usual). You won't
        be able to set an outlying scan cycle as reference. More than
        one cycle can be specified, though, and the list will be used as
        fall backs. If --cycle is given, DIFFEOMORPHISM_NB will be set
        to scan_cycle.""")

    ####################################################################
    # File handling
    ####################################################################

    file_handling = parser.add_argument_group(
        """File handling""")

    lock_handling = file_handling.add_mutually_exclusive_group()

    lock_handling.add_argument('-r', '--remove-lock',
        action='store_true',
        help="""Remove lock, if file is locked. This is useful, if used
        together with -s/--skip to remove orphan locks.""")

    lock_handling.add_argument('-i', '--ignore-lock',
        action='store_true',
        help="""Ignore lock, if file is locked. Together with -s/--skip
        this will also remove orphan locks.""")

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
        help="""Number of cores to use. Default is the number of cores
        on the machine.""")

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

from ..diffeomorphisms import Image

from ..session import Session

from ..reference import ReferenceMaps

from ..pmap import PopulationMap, pmap_scanner

########################################################################

def call(args):

    ####################################################################
    # Options
    ####################################################################

    remove_lock       = args.remove_lock
    ignore_lock       = args.ignore_lock
    force             = args.force
    skip              = args.skip
    verbose           = args.verbose

    scan_cycle         = args.cycle
    resolution         = args.resolution
    new_diffeomorphism = args.new_diffeomorphism
    diffeomorphism_nb  = args.diffeomorphism_nb

    if (resolution is None) or (resolution == 'native') or np.isclose(resolution, 0):
        resolution = None

    ####################################################################
    # Study
    ####################################################################

    study = get_study(args)

    if study is None:
        print('No study found. Nothing to do.')
        return

    study.set_diffeomorphism(new_diffeomorphism)
    study.set_standard_space('isometric')

    ####################################################################
    # Create iterator
    ####################################################################

    study_iterator = study.iterate(
            'session',
            'reference_maps',
            'result',
            'population_map', new=['population_map'],
            integer_index=True)

    df = study_iterator.df.copy()

    df['locked'] = False

    ####################################################################
    # Wrapper
    ####################################################################

    def wm(index, name, session, reference_maps, population_map, result,
            file_population_map):

        if type(population_map) is Lock:
            if remove_lock or ignore_lock:
                if verbose:
                    print('{}: Remove lock'.format(name.name()))
                population_map.unlock()
                if remove_lock:
                    return
            else:
                if verbose:
                    print('{}: Locked'.format(name.name()))
                return

        elif population_map is not None and not force:
            if verbose:
                print('{}: PopulationMap already exists. Use -f/--force to overwrite'.format(
                    name.name()))
            return

        if skip:
            return

        if verbose:
            print('{}: Lock: {}'.format(name.name(), file_population_map))

        lock = Lock(name, 'fmripop', file_population_map)
        df.ix[index, 'locked'] = True

        dfile = os.path.dirname(file_population_map)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        lock.save(file_population_map)

        ####################################################################
        # Create population map from a session instance
        ####################################################################

        if session is None:
            print('{}: No session found'.format(name.name()))
            df.ix[index,'valid'] = False
            lock.conditional_unlock(df, index, verbose)
            return

        if diffeomorphism_nb == 'identity':
            population_map = pmap_scanner(
                    session=session,
                    resolution=resolution,
                    name=new_diffeomorphism)

            if verbose:
                if resolution:
                    print("""{}:
                    Standard space equals scanner space. Diffeomorphism
                    equals identity. Resolution is ({} mm)**3.""".format(
                        name.name(), resolution))
                else:
                    print("""{}:
                    Standard space equals scanner space. Diffeomorphism
                    equals identity. Resolution is native.""".format(
                        name.name()))

        ####################################################################
        # Create population map from a session and reference instance
        ####################################################################

        elif diffeomorphism_nb == 'scanner':
            if reference_maps is None:
                print('{}: No ReferenceMaps found'.format(name.name()))
                df.ix[index,'valid'] = False
                lock.conditional_unlock(df, index, verbose)
                return

            if scan_cycle is None:
                population_map = pmap_scanner(
                        session=session,
                        reference_maps = reference_maps,
                        resolution=resolution,
                        name=new_diffeomorphism)

                if verbose:
                    if resolution:
                        print("""{}:
                        Standard space equals scanner space. Diffeomorphism
                        maps to the average position of the subject in the
                        scanner. Resolution is ({} mm)**3.""".format(
                            name.name(), resolution))
                    else:
                        print("""{}:
                        Standard space equals scanner space. Diffeomorphism
                        maps to the average position of the subject in the
                        scanner. Resolution is native.""".format(
                            name.name()))

            else:
                try:
                    outlying_cycles = reference_maps.outlying_cycles
                except:
                    if verbose:
                        print("""{}:
                        I have found no information
                        about outlying scan cycles!""".format(name.name()))
                    outlying_cycles = None

                if outlying_cycles is None:
                    scan_cycle_to_use = cycle[0]
                else:
                    if outlying_cycles[scan_cycle].all():
                        df.ix[index,'valid'] = False
                        print("""{}:
                        All suggested reference cycles have been marked as
                        outlying. Unable to proceed. Please specify a
                        different scan cycle (using --cycle) as
                        reference.""".format(name.name()))
                        lock.conditional_unlock(df, index, verbose, True)
                        return
                    elif outlying_cycles[scan_cycle].any():
                        for c, co in zip (scan_cycle, outlying_cycles[scan_cycle]):
                            if co:
                                if verbose:
                                    print("""{}:
                                    Scan cycle {:d} marked as outlying,
                                    using fallback.""".format( name.name(), c))
                            else:
                                scan_cycle_to_use = c
                                break
                    else:
                        scan_cycle_to_use = scan_cycle[0]

                population_map = pmap_scanner(
                        session=session,
                        reference_maps = reference_maps,
                        scan_cycle = scan_cycle_to_use,
                        resolution=resolution,
                        name=new_diffeomorphism)

                population_map.set_nb(Image(
                    reference=session.reference,
                    data=session.data[scan_cycle_to_use],
                    name=session.name.name()+'-{:d}'.format(scan_cycle_to_use)))

                if verbose:
                    if resolution:
                        print("""{}:
                        Standard space equals scanner space,
                        diffeomorphism maps to subject position during
                        scan cycle: {:d}. Resolution is ({} mm)**3.""".format(
                            name.name(), scan_cycle_to_use, resolution))
                    else:
                        print("""{}:
                        Standard space equals scanner space,
                        diffeomorphism maps to subject position during
                        scan cycle: {:d}. Resolution is native.""".format(
                            name.name(), scan_cycle_to_use))

        ####################################################################
        # Create population map instance from a result instance
        ####################################################################

        elif (diffeomorphism_nb == 'fit'):

            if result is None:
                print('{}: No fit found'.format(name.name()))
                df.ix[index,'valid'] = False
                lock.conditional_unlock(df, index, verbose)
                return

            population_map = result.population_map
            population_map.set_vb(template=result.get_field('intercept', 'point'))

            if verbose:
                print('{}: VB space equals reference space as given by fit'.format(
                    name.name()))

        else:
            print('{}: Diffeomorphism type not supported'.format(name.name()))

        if verbose > 2:
            print("""{}:
                {}
                {}""".format(name.name(),
                    population_map.diffeomorphism.describe(),
                    population_map.describe()))

        try:
            if verbose:
                print('{}: Save: {}'.format(name.name(),
                    file_population_map))

            population_map.save(file_population_map)
            df.ix[index,'locked'] = False

        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to create: {}, {}'.format(name.name(),
                file_population_map, e))
            lock.conditional_unlock(df, index, verbose, True)
            return

        return

    ####################################################################

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for index, name, files, instances in study_iterator:
                session         = instances['session']
                reference_maps  = instances['reference_maps']
                population_map  = instances['population_map']
                result          = instances['result']
                file_population_map = files['population_map']
                wm
                pool.apply_async(wm, args=(index, name, session,
                    reference_maps, population_map, result,
                    file_population_map))
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
            for index, name, files, instances in study_iterator:
                session         = instances['session']
                reference_maps  = instances['reference_maps']
                population_map  = instances['population_map']
                result          = instances['result']
                file_population_map = files['population_map']
                wm(index, name, session, reference_maps, population_map,
                        result, file_population_map)
        finally:
            files = df.ix[df.locked, 'population_map'].values
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
