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

    specific.add_argument('--nii',
        default='{cohort}-{id:04d}-{paradigm}-{date}.nii',
        help="""Path to a file that should contain a 4D-image of a fMRI
        session in any file format understood by the NiBabel project,
        e.g, any of ANALYZE (plain, SPM99, SPM2 and later), GIFTI,
        NIfTI1, NIfTI2, MINC1, MINC2, MGH and ECAT as well as Philips
        PAR/REC.  For more details see http://nipy.org/nibabel/.  Please
        note that fmristats has only been tested with Nifti1 files.""")

    specific.add_argument('--nii-prefix',
        default='raw/nii/{paradigm}',
        help = """Prefix for the path in --nii.""")

    foreground_handling = specific.add_mutually_exclusive_group()

    foreground_handling.add_argument('--detect-foreground',
        action='store_true',
        help="""Detect the foreground in the FMRI""")

    foreground_handling.add_argument('--set-foreground',
        action='store_true',
        help="""Set the foreground to the image in FOREGROUND""")

    specific.add_argument('--foreground',
        default='{cohort}-{id:04d}-{paradigm}-{date}-foreground.nii.gz',
        help="""A 4D-image that contains a mask for the foreground in
        the FMRI.""")

    specific.add_argument('--foreground-prefix',
        default='raw/foreground/{cohort}',
        help="""Prefix for the path in --foreground.""")

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

from ..session import Session, fmrisetup

from ..nifti import nii2session

import nibabel as ni

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

    detect_foreground = args.detect_foreground
    set_foreground    = args.set_foreground

    ####################################################################
    # Study
    ####################################################################

    study = get_study(args)

    if study is None:
        print('Nothing to do.')
        return

    study.update_layout({
        'nii' : join(args.nii_prefix, args.nii),
        'foreground' : join(args.foreground_prefix, args.foreground),
        })

    if not 'epi' in study.protocol.columns:
        if args.epi_code is not None:
            assert args.epi_code >= -3, 'epi code must be ≥ -3'
            assert args.epi_code >= 3, 'epi code must be ≤ 3'
            assert args.epi_code != 0, 'epi code must be ≠ 0'
            study.protocol['epi'] = args.epi_code
        else:
            print("""
            You need to provide an EPI code via --epi-code. This code
            defines the direction in which the planes in the FMRI have been
            measured. Visit:

            https://fmristats.github.io/tutorials/getting-started.html""")
            return

    ####################################################################
    # Iterator
    ####################################################################

    study_iterator = study.iterate('session', 'stimulus',
            new=['session', 'nii', 'foreground'],
            integer_index=True)

    df = study_iterator.df.copy()

    df['locked'] = False

    ####################################################################
    # Wrapper
    ####################################################################

    def wm(index, name, session, stimulus, file_session, file_nii,
            file_foreground):

        if type(session) is Lock:
            if remove_lock or ignore_lock:
                if verbose:
                    print('{}: Remove lock'.format(name.name()))
                session.unlock()
                if remove_lock:
                    return
            else:
                if verbose:
                    print('{}: Locked'.format(name.name()))
                return

        elif session is not None and not force:
            if verbose:
                print('{}: Session already exists. Use -f/--force to overwrite'.format(
                    name.name()))
            return

        if skip:
            return

        if verbose:
            print('{}: Lock: {}'.format(name.name(), file_session))

        lock = Lock(name, 'nii2session', file_session)
        df.ix[index, 'locked'] = True

        dfile = os.path.dirname(file_session)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        lock.save(file_session)

        ################################################################
        # Load image data
        ################################################################

        try:
            img = ni.load(file_nii)
            if verbose:
                print('{}: Read: {}'.format(name.name(), file_nii))
        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to read: {}, {}'.format(name.name(),
                file_nii, e))

        if lock.conditional_unlock(df, index, verbose):
            return

        ####################################################################
        # Create session instance
        ####################################################################

        session = nii2session(
                name=name,
                nii=img,
                epi_code=df.loc[index, 'epi'])

        fmrisetup(session = session, stimulus = stimulus)

        if detect_foreground:
            if verbose:
                print('{}: Detect foreground'.format(name.name()))
            session.fit_foreground()

        elif set_foreground:

            try:
                foreground_nii = ni.load(file_foreground)
            except Exception as e:
                df.ix[index,'valid'] = False
                print('{}: Unable to read: {}, {}'.format(name.name(),
                    file_nii, e))

            if verbose:
                print('{}: Set foreground: {}'.format(name.name(),
                    foreground_nii))

            session.set_foreground(foreground_nii.get_data())

        if verbose:
            print('{}: Save: {}'.format(name.name(), file_session))

        session.save(file_session)
        df.ix[index,'locked'] = False

    ####################################################################

    if args.cores == 0:
        args.cores = None

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for index, name, files, instances in study_iterator:
                session         = instances['session']
                stimulus        = instances['stimulus']
                file_session    = files['session']
                file_nii        = files['nii']
                file_foreground = files['foreground']

                if stimulus is None:
                    print('{}: No Stimulus found'.format(name.name()))
                else:
                    pool.apply_async(wm, args=\
                    (index, name, session, stimulus, file_session,
                        file_nii, file_foreground)
                    )

            pool.close()
            pool.join()
        except Exception as e:
            pool.close()
            pool.terminate()
            print('Pool execution has been terminated')
            print(e)
        finally:
            files = df.ix[df.locked, 'session'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)
    else:
        try:
            print('Process protocol entries sequentially')
            for index, name, files, instances in study_iterator:
                session         = instances['session']
                stimulus        = instances['stimulus']
                file_session    = files['session']
                file_nii        = files['nii']
                file_foreground = files['foreground']

                if stimulus is None:
                    print('{}: No Stimulus found'.format(name.name()))
                else:
                    wm(index, name, session, stimulus, file_session,
                            file_nii, file_foreground)
        finally:
            files = df.ix[df.locked, 'session'].values
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
