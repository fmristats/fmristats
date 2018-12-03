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

Fit head movements of a subject in an FMRI session

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
        """Arguments controlling the fit""")

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

    specific.add_argument('--grubbs',
        type=float,
        default=0.1,
        help="""An outlier detection is performed to identify scans
        which may have been acquired during severe head movements. More
        precisely, a Grubbs' outlying test will be performed on the set
        of estimated principle semi axis for each full scan cycle on the
        given level of significance. When using fmririgid to create
        ReferenceMaps, the default is 0.1, and the information of
        outlying scans is saved to disk together with the estimated
        rigid body transformations. Then, when running fmrifit, this
        information is used. When setting --grubbs in fmrifit, outlier
        estimation is performed again.""")

    specific.add_argument('--window-radius',
        type=int,
        nargs='+',
        default=[0],
        help="""Calculate the mean rigid transformation using the this
        window. A --window-radius of 0 does not do anything (the
        default). A WINDOW_RADIUS of n corresponds to all rigid
        transformation n scan cycle before and after each time point
        (including the boundary). Can also be used iteratively by
        providing more than one argument. This effectively produces a
        weighted averaging of rigid transformations.""")

    specific.add_argument('--new-rigids',
            default='pcm',
            help="""Name of the rigid transformations.""")

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

from ..study import Study

from ..reference import ReferenceMaps

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

    cycle         = args.cycle
    grubbs        = args.grubbs
    window_radius = args.window_radius

    ####################################################################
    # Study
    ####################################################################

    study = get_study(args)

    if study is None:
        print('Nothing to do.')
        return

    study.set_rigids(args.new_rigids)

    ####################################################################
    # Iterator
    ####################################################################

    study_iterator = study.iterate('session', 'reference_maps',
            new=['reference_maps'],
            integer_index=True)

    df = study_iterator.df.copy()

    df['locked'] = False

    ####################################################################
    # Wrapper
    ####################################################################

    def wm(index, name, session, reference_maps, file_reference_maps):
        if session is None:
            df.ix[index,'valid'] = False
            print('{}: Unable to open session.'.format(name.name()))
            return

        if type(reference_maps) is Lock:
            if remove_lock or ignore_lock:
                if verbose:
                    print('{}: Remove lock'.format(name.name()))
                reference_maps.unlock()
                if remove_lock:
                    return
            else:
                if verbose:
                    print('{}: Locked'.format(name.name()))
                return

        elif reference_maps is not None and not force:
            if verbose:
                print('{}: ReferenceMaps already exists. Use -f/--force to overwrite'.format(
                    name.name()))
            return

        if skip:
            return

        if verbose:
            print('{}: Lock: {}'.format(name.name(), file_reference_maps))

        lock = Lock(name, 'fmririgid', file_reference_maps)
        df.ix[index, 'locked'] = True

        dfile = os.path.dirname(file_reference_maps)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        lock.save(file_reference_maps)

        ####################################################################
        # Fit rigid body transformations
        ####################################################################

        if verbose:
            print('{}: Start fit of rigid body transformations'.format(name.name()))

        reference_maps = ReferenceMaps(name)
        reference_maps.fit(session)
        reference_maps.reset_reference_space()

        ####################################################################
        # Detect outlying scan cycles
        ####################################################################

        if not np.isclose(grubbs, 1):
            if verbose:
                print('{}: Detect outlying scans'.format(name.name()))
            reference_maps.detect_outlying_scans(grubbs)
            outlying_cycles = reference_maps.outlying_cycles
        else:
            outlying_cycles = None

        ####################################################################
        # Reset reference space
        ####################################################################

        if (outlying_cycles is None) or (cycle is None):
            reference_maps.reset_reference_space(cycle=cycle)
        else:
            if outlying_cycles[cycle].all():
                df.ix[index,'valid'] = False
                print("""{}: All suggested reference cycles have been
                marked as outlying. Unable to proceed. Please specify a
                different scan cycle (using --cycle) as reference.""".format(name.name()))
                lock.conditional_unlock(df, index, verbose, True)
                return
            elif outlying_cycles[cycle].any():
                for c, co in zip (cycle, outlying_cycles[cycle]):
                    if co:
                        print("""{}: Cycle {:d} marked as outlying, using fallback.""".format(
                            name.name(), c))
                    else:
                        print("""{}: Reference cycle is {:d}.""".format(
                            name.name(), c))
                        reference_maps.reset_reference_space(cycle=c)
                        break
            else:
                print("""{}: Reference cycle is {:d}.""".format(name.name(), cycle[0]))
                reference_maps.reset_reference_space(cycle=cycle[0])

        ####################################################################
        # Average rigid body transformations
        ####################################################################

        for r in window_radius:
            if r > 0:
                if verbose:
                    print('{}: Flatten head locations using a radius of {} scan cycles'.format(
                        name.name(), r))
                reference_maps.mean(r)

        ####################################################################
        # Save to disk
        ####################################################################

        try:
            if verbose:
                print('{}: Save: {}'.format(name.name(), file_reference_maps))

            reference_maps.save(file_reference_maps)
            df.ix[index,'locked'] = False

        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to create: {}'.format(name.name(), file_reference_maps))
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
                session         = instances['session']
                reference_maps  = instances['reference_maps']
                file_reference_maps = files['reference_maps']

                if session is None:
                    print('{}: No Session found'.format(name.name()))
                else:
                    pool.apply_async(wm, args=\
                    (index, name, session, reference_maps, file_reference_maps)
                    )

            pool.close()
            pool.join()
        except Exception as e:
            pool.close()
            pool.terminate()
            print('Pool execution has been terminated')
            print(e)
        finally:
            files = df.ix[df.locked, 'reference_maps'].values
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
                file_reference_maps = files['reference_maps']

                if session is None:
                    print('{}: No Session found'.format(name.name()))
                else:
                    wm(index, name, session, reference_maps, file_reference_maps)
        finally:
            files = df.ix[df.locked, 'reference_maps'].values
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
