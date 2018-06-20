# Copyright 2016-2017 Thomas W. D. Möbius
#
# This filename is part of fmristats.
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

Wrapper to ANT's RegistrationSyNQuick

"""

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

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

    parser.add_argument('--vb',
            default='mni152',
            help=hp.vb_name)

    parser.add_argument('--vb-image',
            help=hp.vb)

    parser.add_argument('--vb-nii',
            default='/usr/share/data/fsl-mni152-templates/MNI152_T1_2mm_brain.nii.gz',
            help=hp.vb_nii)

    parser.add_argument('--vb-name',
            default='MNI152_T1_2mm_brain',
            help=hp.vb_name)

    parser.add_argument('--vb-path',
            default='../data/vbs/{}.image',
            help=hp.vb_path)

    parser.add_argument('--vb-background-image',
            help=hp.vb_background)

    parser.add_argument('--vb-background-nii',
            default='/usr/share/data/fsl-mni152-templates/MNI152_T1_2mm.nii.gz',
            help=hp.vb_background_nii)

    parser.add_argument('--vb-background-name',
            default='MNI152_T1_2mm',
            help=hp.vb_background_name)

    parser.add_argument('--vb-background-path',
            default='../data/vbs/{}.image',
            help=hp.vb_background_path)

    parser.add_argument('--diffeomorphism',
            default='ants',
            help="""Name to use for the fitted diffeomorphisms.""")

########################################################################
# Output arguments
########################################################################

    parser.add_argument('--population-map',
            default='../data/pop/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}.pop',
            help=hp.population_map)

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-ants4pop.pkl',
            help=hp.protocol_log)

########################################################################
# Additional input arguments when warp coefficient files provided
########################################################################

    parser.add_argument('--output-prefix',
            default='ants/{2}/{4}/{1:04d}/{0}-{1:04d}-{2}-{3}-{4}-',
            help="""Prefix for temporary files from ANTS""")

########################################################################
# Additional input arguments when no warp coefficient files provided
########################################################################

    parser.add_argument('--session',
            default='../data/ses/{2}/{0}-{1:04d}-{2}-{3}.ses',
            help=hp.session)

    parser.add_argument('--cycle',
            type=int,
            help="""cycle to pick as reference""")

    parser.add_argument('--fit',
            default='../data/fit/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}-{5}.fit',
            help="""if working with session files, you don't need this""" + hp.sfit)

    parser.add_argument('--scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help="""only needed if part of the template for --fit""" + hp.scale_type)

    parser.add_argument('--nb-name',
            default='self',
            help="""name of the population space that was originally
            used for the fit"""  + hp.nb_name)

########################################################################
# Arguments specific for using the protocol API
########################################################################

    to_process = parser.add_argument_group(
            """specifying the protocol entries to process""",
            """Arguments which give control which protocol entries to
            process. If no protocol filename is given, it will be checked
            if files being processed comply to the given information.""")

    to_process.add_argument('--protocol',
            help=hp.protocol)

    to_process.add_argument('--cohort',
            help=hp.cohort)

    to_process.add_argument('--id',
            type=int,
            nargs='+',
            help=hp.j)

    to_process.add_argument('--datetime',
            help=hp.datetime)

    to_process.add_argument('--paradigm',
            help=hp.paradigm)

    to_process.add_argument('--strftime',
            default='%Y-%m-%d-%H%M',
            help=hp.strftime)

########################################################################
# Miscellaneous
########################################################################

    lock_handling = parser.add_mutually_exclusive_group()

    lock_handling.add_argument('--remove-lock',
            action='store_true',
            help=hp.remove_lock.format('population space'))

    lock_handling.add_argument('--ignore-lock',
            action='store_true',
            help=hp.ignore_lock.format('population space'))

    file_handling = parser.add_mutually_exclusive_group()

    file_handling.add_argument('-f', '--force',
            action='store_true',
            help=hp.force.format('population space'))

    file_handling.add_argument('-s', '--skip',
            action='store_true',
            help=hp.skip.format('population space'))

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

from ...load import load_session, load_result, load_population_map

from ...name import Identifier

from ...protocol import layout_sdummy, layout_dummy

from ...diffeomorphisms import Image

from ...nifti import image2nii, nii2image

from ...ants import fit_population_map

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

import nibabel as ni

########################################################################

def call(args):

    if args.verbose:
        print('-------------------------------------')

    if args.vb_image:
        try:
            if args.verbose:
                print('Read: {}'.format(args.vb_image))
            vb = load(args.vb_image)
        except Exception as e:
            print('Unable to read vb filename: {}'.format(args.vb_image))
            print('Fall back to --vb-nii: {}'.format(e))
    else:
        try:
            if args.verbose:
                print('Read: {}'.format(args.vb_nii))
            vb = nii2image(ni.load(args.vb_nii), name=args.vb_name)

            fname = args.vb_path.format(args.vb_name)

            dfile = os.path.dirname(fname)
            if dfile and not isdir(dfile):
               os.makedirs(dfile)

            if args.verbose:
                print('Save: {}'.format(fname))
            vb.save(fname)
        except Exception as e:
            print('Unable to read or write vb files')
            print('Exception: {}'.format(e))
            exit()

    if args.vb_background_image:
        try:
            if args.verbose:
                print('Read: {}'.format(args.vb_background_image))
            vb_background = load(args.vb_background_image)
        except Exception as e:
            print('Unable to read vb-background: {}'.format(args.vb_background_image))
            print('Fall back to --vb-background-nii: {}'.format(e))
    else:
        try:
            if args.verbose:
                print('Read: {}'.format(args.vb_background_nii))
            vb_background = nii2image(ni.load(args.vb_background_nii), name=args.vb_background_name)

            fname = args.vb_background_path.format(args.vb_background_name)

            dfile = os.path.dirname(fname)
            if dfile and not isdir(dfile):
               os.makedirs(dfile)

            if args.verbose:
                print('Save: {}'.format(fname))
            vb_background.save(fname)
        except Exception as e:
            print('Unable to read or write vb-background files')
            print('Exception: {}'.format(e))
            exit()

    ####################################################################

    output = args.protocol_log.format(
            datetime.datetime.now().strftime('%Y-%m-%d-%H%M'))

    if args.strftime == 'short':
        args.strftime = '%Y-%m-%d'

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args, fall_back=[args.session, args.fit])

    if df is None:
        sys.exit()

    ####################################################################
    # Add filename layout
    ####################################################################

    df_layout = df.copy()

    layout_dummy(df_layout, 'ses',
            template=args.session,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'output_prefix',
            template=args.output_prefix,
            urname=args.nb_name,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'fit',
            template=args.fit,
            urname=args.nb_name,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'filename',
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
            dfile = os.path.dirname(r.filename)
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
                filename          = r.filename,

                file_res          = r.fit,
                file_ses          = r.ses,

                output_prefix     = r.output_prefix,

                vb                = vb,
                vb_token          = args.vb,
                vb_background     = vb_background,
                nb_name           = args.nb_name,
                cycle             = args.cycle,
                diffeomorphism_name = args.diffeomorphism,
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
            files = df_layout.ix[df.locked, 'filename'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)
            del df['locked']
    else:
        try:
            print('-------------------------------------')
            print('Process protocol entries sequentially')
            for r in it:
                wm(r)
        finally:
            files = df_layout.ix[df.locked, 'filename'].values
            if len(files) > 0:
                for f in files:
                    print('Unlock: {}'.format(f))
                    os.remove(f)
            del df['locked']

    ####################################################################
    # Write protocol
    ####################################################################

    if args.verbose:
        print('-------------------------------------')
        print('Save: {}'.format(output))

    dfile = os.path.dirname(output)
    if dfile and not isdir(dfile):
       os.makedirs(dfile)

    df.to_pickle(output)

########################################################################

def wrapper(name, df, index, remove_lock, ignore_lock, force, skip,
        verbose, filename, file_res, file_ses, output_prefix, vb,
        vb_token, vb_background, nb_name, cycle, diffeomorphism_name):

    if isfile(filename):
        instance = load_population_map(filename, name, df, index, vb.name, verbose)
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
        print('{}: Lock: {}'.format(name.name(), filename))

    lock = Lock(name, 'ants4pop', filename)
    df.ix[index, 'locked'] = True
    lock.save(filename)
    df.ix[index,'valid'] = True

    ####################################################################

    if cycle is None:

        result = load_result(file_res, name, df, index, nb_name, verbose)
        if lock.conditional_unlock(df, index, verbose):
            return

        result.mask()
        nb = result.get_field('intercept', 'point')

    else:

        session = load_session(file_ses, name, df, index, verbose)
        if lock.conditional_unlock(df, index, verbose):
            return

        nb = Image(
            reference=session.reference,
            data=session.data[cycle],
            name=session.name) #.name()+'-{:d}'.format(cycle))

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
                name=diffeomorphism_name)
    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to fit diffeomorphism'.format(name.name()))
        print('{}: Exception: {}'.format(name.name(), e))
        return

    population_map.set_vb_background(vb_background)

    if verbose:
        print('{}: Save: {}'.format(name.name(), filename))

    try:
        population_map.save(filename)
        df.ix[index,'locked'] = False
    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to write: {}'.format(name.name(), filename))
        print('{}: Exception: {}'.format(name.name(), e))
        lock.conditional_unlock(df, index, verbose, True)
