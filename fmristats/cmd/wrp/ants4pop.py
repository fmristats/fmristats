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

Set the standard space to the given image and fit a diffeomorphism from
standard space to subject space using a wrapper to ANTS

"""

#import warnings
#warnings.simplefilter(action='ignore', category=FutureWarning)

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
        """Creating a standard space isometric to the reference space.""")

    specific.add_argument('--diffeomorphism-type',
        default='reference',
        choices=['reference', 'scanner', 'scan_cycle', 'fit'],
        help="""Type of diffeomorphism.""")

    specific.add_argument('--resolution',
        default=2.,
        type=float,
        help="""The resolution of the template in standard
        space.""")

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
        fall backs.""")

    specific.add_argument('--new-diffeomorphism',
            default='ants',
            help="""Name to use for the fitted diffeomorphisms.""")

    specific.add_argument('--ants4pop-prefix',
            default='ants/vb/{}.nii.gz'
            help="""Prefix for temporary files.""")

    specific.add_argument('--ants-prefix',
            default='ants/tmp/{2}/{4}/{1:04d}/{0}-{1:04d}-{2}-{3}-{4}-',
            help="""Prefix for temporary ANTS' files.""")

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

from ...session import Session

from ...reference import ReferenceMaps

from ...pmap import PopulationMap

from ...nifti import image2nii, nii2image

from ...ants import fit_population_map

import nibabel as ni

########################################################################

def call(args):

    study = get_study(args)

    if study is None:
        sys.exit()

    study_iterator = study.iterate(
            'session',
            'reference_maps',
            'result',
            'population_map',
            new=['population_map', 'ants_prefix'],
            integer_index=True)

    df = study_iterator.df.copy()

    df['locked'] = False

    remove_lock       = args.remove_lock
    ignore_lock       = args.ignore_lock
    force             = args.force
    skip              = args.skip
    verbose           = args.verbose

    diffeomorphism_type = args.diffeomorphism_type
    resolution          = args.resolution
    cycle               = args.cycle

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
        # Create population map instance from a result instance
        ####################################################################

        if cycle is None:

            result.mask()
            nb = result.get_field('intercept', 'point')
            nb.name = name.name()

        ####################################################################
        # Create population map instance from a session instance
        ####################################################################

        else:

            if reference_maps is None:
                scan_cycle_to_use = cycle[0]
            else:
                try:
                    outlying_cycles = reference_maps.outlying_cycles
                except:
                    if verbose:
                        print('{}: I have found no information about outlying scan cycles!'.format(
                            name.name()))
                    outlying_cycles = None

                if outlying_cycles is None:
                    scan_cycle_to_use = cycle[0]
                else:
                    if outlying_cycles[cycle].all():
                        df.ix[index,'valid'] = False
                        print("""{}:
                        All suggested reference cycles have been marked as
                        outlying. Unable to proceed. Please specify a
                        different scan cycle (using --cycle) as
                        reference.""".format(name.name()))
                        lock.conditional_unlock(df, index, verbose, True)
                        return
                    elif outlying_cycles[cycle].any():
                        for c, co in zip (cycle, outlying_cycles[cycle]):
                            if co:
                                if verbose:
                                    print("""{}: Suggested cycle {:>4d} marked as outlying, using fallback.""".format(
                                        name.name(), c))
                            else:
                                scan_cycle_to_use = c
                                break
                    else:
                        scan_cycle_to_use = cycle[0]

            if verbose:
                print('{}: NB equals subject position during scan cycle: {:d}'.format(
                    name.name(), scan_cycle_to_use))

            nb = Image(
                reference=session.reference,
                data=session.data[scan_cycle_to_use],
                name=name.name())

        if verbose:
            print('{}: Population space: {}'.format(name.name(), vb.name))
            print('{}: Subject space: {}'.format(name.name(), nb.name))
            print('{}: Fit diffeomorphism'.format(name.name()))

        try:
            population_map = fit_population_map(
                    vb_image=vb,
                    nb_image=nb,
                    output_prefix = output_prefix,
                    vb=vb_token,
                    name=diffeomorphism_new)
        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to fit diffeomorphism'.format(name.name()))
            print('{}: Exception: {}'.format(name.name(), e))
            return

        population_map.set_vb_background(study.vb_background)

        if verbose:
            print('{}: Save: {}'.format(name.name(), file_population_map))

        try:
            population_map.save(file_population_map)
            df.ix[index,'locked'] = False
        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to write: {}'.format(name.name(), file_population_map))
            print('{}: Exception: {}'.format(name.name(), e))
            lock.conditional_unlock(df, index, verbose, True)

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
