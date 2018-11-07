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

    input_files = parser.add_argument_group(
            """input files""",
            """Paths or templates for the paths to input files.""")

    input_files.add_argument('--session',
            default='../data/ses/{2}/{0}-{1:04d}-{2}-{3}.ses',
            help='input file;' + hp.session)

    # TODO: also give the option --reference-maps None or none,
    # which will assume no movement. This is useful for analysing
    # phantom data.

    input_files.add_argument('--reference-maps',
            default='../data/ref/{2}/{0}-{1:04d}-{2}-{3}.ref',
            help='input file;' + hp.reference_maps)

    input_files.add_argument('--population-map',
            default='../data/pop/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}.pop',
            help=hp.population_map)

    input_files.add_argument('--vb',
            default='self',
            help=hp.vb)

    parser.add_argument('--diffeomorphism',
            default='reference',
            help="""Name of the diffeomorphism.""")

########################################################################
# Output arguments
########################################################################

    output_files = parser.add_argument_group(
            """output files""",
            """Paths or templates for the paths to output field.""")

    output_files.add_argument('--fit',
            default='../data/fit/{2}/{4}/{5}/{6}/{0}-{1:04d}-{2}-{3}-{4}.fit',
            help=hp.sfit)

########################################################################
# Log file
########################################################################

    output_files.add_argument('-o', '--protocol-log',
            default='logs/{}-fmrifit.pkl',
            help=hp.protocol_log)

########################################################################
# Arguments specific for the RSM Signal Model: design matrix
########################################################################

    signal_model = parser.add_argument_group(
            """define the signal model""",
            """Model specifications of the signal model""")

    signal_model.add_argument('--formula',
            default='C(task)/C(block, Sum)',
            help="""formula""")

    signal_model.add_argument('--parameter',
            default=['intercept', 'task'],
            nargs='+',
            help="""parameter""")

########################################################################
# Arguments specific for the RSM Signal Model: stimulus II
########################################################################

    experimental_design = parser.add_argument_group(
            """define the experimental design""",
            """Arguments which define the experimental design of the
            session.""")

    experimental_design.add_argument('--stimulus-block',
            default='stimulus',
            help=hp.stimulus_block)

    experimental_design.add_argument('--control-block',
            default='control',
            help=hp.control_block)

    experimental_design.add_argument('--acquisition-burn-in',
            default=4,
            type=int,
            help=hp.burn_in)

    experimental_design.add_argument('--offset-beginning',
            type=float,
            default=5.242,
            help=hp.offset_beginning)

    experimental_design.add_argument('--offset-end',
            type=float,
            default=1.242,
            help=hp.offset_end)

########################################################################
# Arguments specific for the RSM Signal Model: weighting kernel
########################################################################

    weighting = parser.add_argument_group(
            """control the weighting of observations""",
            """Arguments which control the weighting scheme used in
            estimation.""")

    weighting.add_argument('--scale',
            type=float,
            help=hp.scale)

    weighting.add_argument('--scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help=hp.scale_type)

    weighting.add_argument('--factor',
            type=float,
            default=3,
            help=hp.factor)

    weighting.add_argument('--mass',
            type=float,
            help=hp.mass)

########################################################################
# Arguments specific for the RSM Signal Model: where to fit
########################################################################

    where_to_fit = parser.add_argument_group(
            """where to fit""",
            """By default, fmrifit will respect the brain mask which is
            saved in the population map, and it will fit the model at
            all points within this mask. This default behaviour can, of
            course, be changed.""")

    where_to_fit.add_argument('--mask',
            help="""mask to use""")

    where_to_fit.add_argument('--ignore-mask',
            action='store_true',
            help="""Ignore any brain masks saved in the respective
            population map.""")

    where_to_fit.add_argument('--at-slice',
            type=int,
            nargs='+',
            help="""only fit the model at a slice of the index lattice""")

    # TODO: mask: not implemented yet
    #where_to_fit.add_argument('--mask',
    #        help="""an image file""")

########################################################################
# Arguments specific for using the protocol API
########################################################################

    to_process = parser.add_argument_group(
            """specifying the protocol entries to process""",
            """Arguments which give control which protocol entries to
            process. If no protocol file is given, it will be checked
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

    handling = parser.add_argument_group(
            """file handling""",
            """You can either ignore the lock, remove the lock, or
            respect the lock of a file. If you choose to respect the
            lock (default), then protocol entries that are locked will
            be skipped.""")

    lock_handling = handling.add_mutually_exclusive_group()

    lock_handling.add_argument('-r', '--remove-lock',
            action='store_true',
            help=hp.remove_lock.format('fit'))

    lock_handling.add_argument('-i', '--ignore-lock',
            action='store_true',
            help=hp.ignore_lock.format('fit'))

    file_handling = handling.add_mutually_exclusive_group()

    file_handling.add_argument('-f', '--force',
            action='store_true',
            help=hp.force.format('fit'))

    file_handling.add_argument('-s', '--skip',
            action='store_true',
            help=hp.skip.format('fit'))

########################################################################
# Arguments
########################################################################

    detect = parser.add_argument_group(
            """detections""",
            """If these things have not already happend, you may perform
            these operations know.""")

    detect.add_argument('--grubbs',
            type=float,
            help=hp.grubbs)

    detect.add_argument('--detect-foreground',
            action='store_true',
            help=hp.detect_foreground)

########################################################################
# Multiprocessing
########################################################################

    misc = parser.add_argument_group(
            """miscellaneous""")

    misc.add_argument('-v', '--verbose',
            action='count',
            default=0,
            help=hp.verbose)

    misc.add_argument('-j', '--cores',
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

    df = get_df(args, fall_back=args.session)

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

    layout_dummy(df_layout, 'ref',
            template=args.reference_maps,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'pop',
            template=args.population_map,
            urname=args.vb,
            scale_type=args.diffeomorphism,
            strftime=args.strftime
            )

    layout_fdummy(df_layout, 'file',
            template=args.fit,
            vb=args.vb,
            diffeo=args.diffeomorphism,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    ####################################################################
    # Respect the mask mask
    ####################################################################

    if args.ignore_mask:
        mask = False
    elif args.mask:
        mask = args.mask
    else:
        mask = True

    print('Mask: {}'.format(mask))

    ####################################################################
    # Fit at slice
    ####################################################################

    if args.at_slice:
        fit_at_slice = True
        slice_object = (slice(args.at_slice[0], args.at_slice[1]),
                        slice(args.at_slice[2], args.at_slice[3]),
                        slice(args.at_slice[4], args.at_slice[5]))
    else:
        fit_at_slice = False
        slice_object = None

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

                file_ses = r.ses,
                file_ref = r.ref,
                file_pop = r.pop,

                vb                   = args.vb,
                stimulus_block       = args.stimulus_block,
                control_block        = args.control_block,
                scale_type           = args.scale_type,
                scale                = args.scale,
                factor               = args.factor,
                mass                 = args.mass,
                mask                 = mask,
                offset               = args.offset_beginning,
                preset               = args.offset_end,
                sgnf                 = args.grubbs,
                detect_foreground    = args.detect_foreground,
                burn_in              = args.acquisition_burn_in,
                formula              = args.formula,
                parameter            = args.parameter,

                fit_at_slice         = fit_at_slice,
                slice_object         = slice_object,
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
        verbose, file, file_ses, file_ref, file_pop, vb,
        stimulus_block, control_block, scale_type, scale, factor, mass,
        mask,
        offset, preset, sgnf,
        detect_foreground, burn_in,
        formula, parameter,
        fit_at_slice, slice_object):

    if isfile(file):
        instance = load_result(file, name, df, index, vb, verbose)
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
    # Load session instance from disk
    ####################################################################

    session = load_session(file_ses, name, df, index, verbose)
    if lock.conditional_unlock(df, index, verbose):
        return

    if detect_foreground:
        if verbose:
            print('{}: Detect foreground'.format(name.name()))
        session.fit_foreground()

    ####################################################################
    # Load reference maps from disk
    ####################################################################

    reference_maps = load_refmaps(file_ref, name, df, index, verbose)
    if lock.conditional_unlock(df, index, verbose):
        return

    if (sgnf is not None) and (not np.isclose(sgnf, 1)):
        if verbose:
            print('{}: Detect outlying scans'.format(name.name()))
        reference_maps.detect_outlying_scans(sgnf)

    ####################################################################
    # Load or created population map from disk (as needed)
    ####################################################################

    population_map = load_population_map(file_pop, name, df, index, None, verbose)
    if lock.conditional_unlock(df, index, verbose):
        return

    if verbose:
        print('{}: Standard space is: {}'.format(name.name(),
            population_map.name))
        print('{}: Diffeomorphism is: {}'.format(name.name(),
            population_map.diffeomorphism.name))

    ########################################################################
    # Create signal model instance
    ########################################################################

    if verbose:
        print('{}: Configure signal model'.format(name.name()))

    smodel = SignalModel(
        session=session,
        reference_maps=reference_maps,
        population_map=population_map)

    smodel.create_stimulus_design(
            s=stimulus_block,
            c=control_block,
            offset=offset,
            preset=preset)

    smodel.create_data_matrix(burn_in=burn_in, verbose=verbose)

    smodel.set_hyperparameters(
            scale_type=scale_type,
            scale=scale,
            factor=factor,
            mass=mass)

    if verbose > 1:
        print(smodel.describe())

    if verbose:
        print('{}: Create design matrix'.format(name.name()))

    smodel.create_design_matrix(
            formula=formula,
            parameter=parameter,
            verbose=verbose)

    ########################################################################
    # Create signal model instance
    ########################################################################

    if fit_at_slice:
        if verbose:
            print('{}: Fit at {}'.format(name.name(), slice_object))

        coordinates = smodel.population_map.diffeomorphism.coordinates()
        coordinates = coordinates[slice_object]

        result = smodel.fit_at_subject_coordinates(coordinates,
                mask=mask, verbose=verbose)

        if verbose:
            print('{}: Done fitting'.format(name.name()))

    else:
        result = smodel.fit(mask=mask, verbose=verbose)

        if verbose:
            print('{}: Done fitting'.format(name.name()))

    try:
        if verbose:
            print('{}: Save: {}'.format(name.name(), file))

        result.save(file)
        df.ix[index,'locked'] = False

    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to create: {}'.format(name.name(), file))
        print('{}: Exception: {}'.format(name.name(), e))
        lock.conditional_unlock(df, index, verbose, True)
        return

    return
