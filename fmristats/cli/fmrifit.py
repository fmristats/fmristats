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

Fit FMRI data

"""

    # TODO: also give the option --reference-maps None or none,
    # which will assume no movement. This is useful for analysing
    # phantom data.

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
    # Setting the design matrix
    ####################################################################

    signal_model = parser.add_argument_group(
        """Setting the design matrix of the model""")

    signal_model = parser.add_argument_group(
        """Define the design matrix of the model""",
        """The following two arguments allow you to set the model
        specifications of the signal model.""")

    signal_model.add_argument('--formula',
        default='C(task)/C(block, Sum)',
        help="""A formula defining the design matrix.""")

    signal_model.add_argument('--parameter',
        default=['intercept', 'task'],
        nargs='+',
        help="""The parameters to keep in the result field.""")

    signal_model.add_argument('--window-radius',
        type=int,
        nargs='+',
        default=[0],
        help="""Calculate the mean rigid transformation using the this
        window. A --window-radius of 0 does not do anything (the
        default). A WINDOW_RADIUS of n corresponds to all rigid
        transformation n scan cycle before and after each time point
        (including the boundary).""")

    ####################################################################
    # Setting the stimulus
    ####################################################################

    experimental_design = parser.add_argument_group(
        """Define the experimental design""",
        """Define the experimental design of the model.""")

    experimental_design.add_argument('--stimulus-block',
        default='stimulus',
        help="""Name of the stimulus block""")

    experimental_design.add_argument('--control-block',
        default='control',
        help="""Name of the control block""")

    experimental_design.add_argument('--acquisition-burn-in',
        default=4,
        type=int,
        help="""Acquisition burn in. Defines the *first* scan cycle that
        should be considered during the fit. Prior scan cycles are
        discarded.""")

    experimental_design.add_argument('--offset-beginning',
        type=float,
        default=5.242,
        help = """The haemodynamic response to stimulus is not
        immediate.  It is usually assumed that the HR-function spikes
        approximately five seconds after the first stimulus.  The value
        of OFFSET_BEGINNING is added to the onset times of the stimulus
        phases to allow you to wait till the subject's brain is in the
        respected stimulus modus.""")

    experimental_design.add_argument('--offset-end',
        type=float,
        default=1.242,
        help = """ Similar to OFFSET_BEGINNING the value OFFSET_END is
        removed from the end of each stimulus phase and not considered
        in the fitting.""")

    ####################################################################
    # Weighting kernel
    ####################################################################

    weighting = parser.add_argument_group(
        """Control the spatial weighting of observations""",
        """Defines the spatial weighting scheme used in estimation.""")

    weighting.add_argument('--scale',
        type=float,
        help = """Standard deviation of a Gaussian kernel that defines
        the weighting scheme of the underlying WLS regression.  If not
        given explicitly, SCALE will be set to one half of the length of
        the diagonal of the orthorhombic measure lattice of the
        acquisition grid of the session (and this is recommended). This
        default behaviour corresponds to the SCALE_TYPE ``diagonal`` and
        can be overwritten by setting SCALE_TYPE to a different value.
        The parameter SCALE will determine the final curvature of the
        fitted effect field. The larger SCALE, the flatter the fitted
        effect field will appear.""")

    weighting.add_argument('--factor',
        type=float,
        default=3,
        help = """Only observations within an area of FACTOR standard
        deviations under the Gaussian distribution defined by SCALE will
        find their way into the WLS regression. The rational behind this
        parameter is that weights outside of FACTOR standard deviations
        will be very small and will not contribute to the actual
        estimate at any given point. They will only slow down
        computation.""")

    weighting.add_argument('--mass',
        type=float,
        help = """Only observations within an area of mass `mass` under
        the Gaussian distribution defined by SCALE will find their way
        into the WLS regression.  If specified, must be a value between
        0 and 1, and should be a value close to 1.""")

    ####################################################################
    # Arguments specific for the RSM Signal Model: where to fit
    ####################################################################

    where_to_fit = parser.add_argument_group(
        """Control in which areas to fit the field""",
        """By default, fmrifit will respect the template which is saved
        in the population map, and it will fit the model only at points
        which are non-null in the template. This only makes sense,
        though, if you are sure that the diffeomorphism that maps from
        standard space to the subject reference space is really correct.
        If this is the case it will considerably speed up computation
        time. If not, however, it will make a visual control of the fit
        more difficult. This is why this default behaviour can be be
        changed.""")

    where_to_fit.add_argument('--mask',
        help="""Set the mask to use. If yes, true, apply or not given,
        both vb and vb_mask will apply. If no, false or ignore, neither
        mask will be applied. If vb, vb_background, vb_estimate or
        vb_mask, the respective mask will be applied.""")

    where_to_fit.add_argument('--at-slice',
        type=int,
        nargs='+',
        help="""Only fit the model at a slice of the index lattice""")

    ####################################################################
    # Detections
    ####################################################################

    detect = parser.add_argument_group(
        """Detections""",
        """If these things have not already happend, you may perform
        these operations now. In particular, you may want to set a
        window size for averaging the head movements.""")

    detect.add_argument('--grubbs',
        type=float,
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

    detect.add_argument('--detect-foreground',
        action='store_true',
        help="""Detect the foreground in the FMRI""")

    ####################################################################
    # Backends
    ####################################################################

    backends = parser.add_argument_group(
        """Backends""",
        """Currently, two backends are supported: Numba and
        Statsmodels.""")

    backends.add_argument('--backend',
        default='numba',
        choices=['numba', 'jit', 'statsmodels'],
        help="""Choose your backend. Backends numba and jit are two
        names for the same backend. JIT is fast, statsmodels is slow.
        Statsmodels also calculates a Durbin-Watson type statistics.""")

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

    control_push  = parser.add_argument_group(
        """Control whether to save the modified (thus overwrite the
        existing) study instance.""")

    control_push.add_argument('-p', '--push',
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

from ..session import Session

from ..reference import ReferenceMaps

from ..pmap import PopulationMap

from ..smodel import SignalModel

########################################################################

def call(args):

    study = get_study(args)

    if study is None:
        print('Nothing to do.')
        return

    ####################################################################
    # Respect the mask mask
    ####################################################################

    if (args.mask is None) or (args.mask == 'yes') or \
            (args.mask == 'true') or (args.mask == 'apply'):
        mask = True
    elif (args.mask == 'no') or (args.mask == 'false') or \
            (args.mask == 'none') or (args.mask == 'ignore'):
        mask = False
    else:
        mask = args.mask

    if args.verbose:
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
    # Options
    ####################################################################

    remove_lock       = args.remove_lock
    ignore_lock       = args.ignore_lock
    force             = args.force
    skip              = args.skip
    verbose           = args.verbose

    scale_type           = study.scale_type

    stimulus_block       = args.stimulus_block
    control_block        = args.control_block
    scale                = args.scale
    factor               = args.factor
    mass                 = args.mass
    offset               = args.offset_beginning
    preset               = args.offset_end
    sgnf                 = args.grubbs
    window_radius        = args.window_radius
    detect_foreground    = args.detect_foreground
    burn_in              = args.acquisition_burn_in
    formula              = args.formula
    parameter            = args.parameter

    mask                 = mask
    fit_at_slice         = fit_at_slice
    slice_object         = slice_object

    backend              = args.backend

    ####################################################################
    # Create the iterator
    ####################################################################

    study_iterator = study.iterate(
            'session',
            'reference_maps',
            'result',
            'population_map',
            new=['result'],
            integer_index=True,
            verbose=verbose)

    df = study_iterator.df.copy()

    df['locked'] = False

    ####################################################################
    # Wrapper
    ####################################################################

    def wm(index, name, session, reference_maps, population_map, result,
            file_result):

        if type(result) is Lock:
            if remove_lock or ignore_lock:
                if verbose:
                    print('{}: Remove lock'.format(name.name()))
                result.unlock()
                if remove_lock:
                    return
            else:
                if verbose:
                    print('{}: Locked'.format(name.name()))
                return

        elif result is not None and not force:
            if verbose:
                print('{}: Result already exists. Use -f/--force to overwrite'.format(
                    name.name()))
            return

        if skip:
            return

        if verbose:
            print('{}: Lock: {}'.format(name.name(), file_result))

        lock = Lock(name, 'fmrifit', file_result)
        df.ix[index, 'locked'] = True

        dfile = os.path.dirname(file_result)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        lock.save(file_result)

        ###############################################################
        # Detect foreground
        ###############################################################

        if detect_foreground:
            if verbose:
                print('{}: Detect foreground'.format(name.name()))
            session.fit_foreground()

        ####################################################################
        # Detect outlying scan cycles
        ####################################################################

        if (sgnf is not None) and (not np.isclose(sgnf, 1)):
            if verbose:
                print('{}: Detect outlying scans'.format(name.name()))
            reference_maps.detect_outlying_scans(sgnf)

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
        # Inform about the standard space and diffeomorphism
        ####################################################################

        if verbose:
            print("""{}:
            Standard space: {}
            Diffeomorphism: {}""".format(name.name(),
                population_map.vb.name,
                population_map.diffeomorphism.name))

        ########################################################################
        # Defining the signal model
        ########################################################################

        if verbose:
            print('{}: Configure signal model'.format(name.name()))

        smodel = SignalModel(
            session=session,
            reference_maps=reference_maps,
            population_map=population_map)

        smodel.set_stimulus_design(
                s=stimulus_block,
                c=control_block,
                offset=offset,
                preset=preset)

        smodel.set_hyperparameters(
                scale_type=scale_type,
                scale=scale,
                factor=factor,
                mass=mass)

        smodel.set_data(burn_in=burn_in, verbose=verbose)

        smodel.set_design(
                formula=formula,
                parameter=parameter,
                verbose=verbose)

        if verbose > 2:
            print("""{}: {}""".format(name.name(), smodel.describe()))

        ########################################################################
        # Fit
        ########################################################################

        if fit_at_slice:
            if verbose:
                print('{}: Fit at {}'.format(name.name(), slice_object))

            coordinates = smodel.population_map.diffeomorphism.coordinates()
            coordinates = coordinates[slice_object]

            result = smodel.fit_at_subject_coordinates(
                    coordinates=coordinates,
                    mask=mask,
                    verbose=verbose,
                    backend=backend)

            if verbose:
                print('{}: Done fitting'.format(name.name()))

        else:
            result = smodel.fit(
                    mask=mask,
                    verbose=verbose,
                    backend=backend)

            if verbose:
                print('{}: Done fitting'.format(name.name()))

        ###############################################################
        # Save the result to disk
        ###############################################################

        try:
            if verbose:
                print('{}: Save: {}'.format(name.name(),
                    file_result))

            result.save(file_result)
            df.ix[index,'locked'] = False

        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to create: {}, {}'.format(name.name(),
                file_result, e))
            lock.conditional_unlock(df, index, verbose, True)
            return

        if verbose > 2:
            print("""{}: {}""".format(name.name(), result.describe()))

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
                population_map  = instances['population_map']
                result          = instances['result']
                file_result     = files['result']

                skip = False
                if session is None:
                    print('{}: No Session found'.format(name.name()))
                    skip = True
                if reference_maps is None:
                    print('{}: No ReferenceMaps found'.format(name.name()))
                    skip = True
                if population_map is None:
                    print('{}: No PopulationMap found'.format(name.name()))
                    skip = True
                if not skip:
                    pool.apply_async(wm, args=(index, name, session,
                        reference_maps, population_map, result,
                        file_result))
            pool.close()
            pool.join()
        except Exception as e:
            pool.close()
            pool.terminate()
            print('Pool execution has been terminated')
            print(e)
        finally:
            files = df.ix[df.locked, 'result'].values
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
                file_result     = files['result']

                skip = False
                if session is None:
                    print('{}: No Session found'.format(name.name()))
                    skip = True
                if reference_maps is None:
                    print('{}: No ReferenceMaps found'.format(name.name()))
                    skip = True
                if population_map is None:
                    print('{}: No PopulationMap found'.format(name.name()))
                    skip = True
                if not skip:
                    wm(index, name, session, reference_maps,
                            population_map, result, file_result)
        finally:
            files = df.ix[df.locked, 'result'].values
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
