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
standard space to subject space using a wrapper to FSL_ FNIRT_.

On the one hand, this is a wrapper to ``std2imgcoord`` and converts a
given FSL warp coefficient file produced by FNIRT_ to the population map
of the corresponding session. On the other hand, it can be used as a
wrapper to the FSL_ command line tools FNIRT_ to estimate this warp
coefficient file.

.. _FNIRT: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT

.. _FSL: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/

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
        choices=['scan_cycle', 'fit'],
        help="""Type of diffeomorphism.""")

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
        default='fnirt',
        help="""Name to use for the fitted diffeomorphisms.""")

    specific.add_argument('--fnirt-prefix',
            default='warping-by-fnirt/{cohort}-{id:04d}/{cohort}-{id:04d}-{paradigm}-{date}-{space}-',
            help="""Prefix for FNIRT files.""")

    ########################################################################
    # Additional input arguments when warp coefficient files are
    ########################################################################

    specific.add_argument('--fnirt-subject-reference-space',
        default='nb.nii.gz',
        help="""Image used as template for the subject in subject
        reference space""")

    specific.add_argument('--fnirt-spline-coefficients',
        default='warpcoef.nii.gz',
        help="""Spline or warp coefficient file created by FNIRT.""")

    specific.add_argument('--ignore-existing-spline-coefficients',
        action='store_true',
        help="""By default existing spline coefficient files will be
        used, as it is assumed that you have created them manually
        to fit your needs. When this flag is set, it will re-fit the
        spline coefficients using some reasonable defaults.""")

    specific.add_argument('--cmd-fnirt',
        default='fsl5.0-fnirt',
        help="""FSL fnirt command. Must be in your path.""")

    specific.add_argument('--cmd-config',
        default='T1_2_MNI152_2mm',
        help="""Config file for FSL FNIRT.""")

    specific.add_argument('--cmd-std2imgcoord',
        default='fsl5.0-std2imgcoord',
        help="""FSL std2imgcoord command. Must be in your path.""")

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

from ..diffeomorphisms import Image

from ..session import Session

from ..reference import ReferenceMaps

from ..pmap import PopulationMap

from ..nifti import image2nii, nii2image

from ..fsl import fnirt, splines2warp

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

    if args.diffeomorphism_nb is None:
        if args.cycle is None:
            diffeomorphism_nb = 'fit'
        else:
            diffeomorphism_nb = 'scan_cycle'
    else:
        diffeomorphism_nb = args.diffeomorphism_nb

    new_diffeomorphism = args.new_diffeomorphism
    cycle              = args.cycle

    cmd_fnirt        = args.cmd_fnirt
    cmd_config       = args.cmd_config
    cmd_std2imgcoord = args.cmd_std2imgcoord

    nb_file_template       = args.fnirt_subject_reference_space
    warpcoef_file_template = args.fnirt_spline_coefficients
    vb_estimate_template   = 'vb-estimate.nii.gz'

    cvb = 'vb.txt'
    cnb = 'nb.txt'

    ignore_existing_warpcoef = args.ignore_existing_spline_coefficients

    ####################################################################
    # Study
    ####################################################################

    study = get_study(args)

    if study is None:
        print('No study found. Nothing to do.')
        return

    if study.vb is None:
        print("""
        You need to provide a template in standard space (for example by
        either setting a template in the study or by providing a
        template using --vb-image or --vb-nii).""")
        return

    study.set_rigids(None)
    study.set_diffeomorphism(new_diffeomorphism)
    study.set_standard_space(study.vb.name)

    study.update_layout({'fnirt_prefix':args.fnirt_prefix})

    ####################################################################
    # Iterator
    ####################################################################

    study_iterator = study.iterate(
            'session',
            'reference_maps',
            'result',
            'population_map',
            new=['population_map', 'fnirt_prefix'],
            lookup=['result'],
            integer_index=True,
            verbose=verbose)

    df = study_iterator.df.copy()

    df['locked'] = False

    ####################################################################
    # Save or look for template in standard space
    ####################################################################

    if args.vb_nii is None:
        vb_nii = 'warping-by-fnirt/vb.nii.gz'

        if verbose:
            print('Save template in standard space to: {}'.format(vb_nii))

        dfile = os.path.dirname(vb_nii)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        ni.save(image2nii(study.vb), vb_nii)
    else:
        vb_nii = args.vb_nii

        if verbose:
            print('Use existing template in standard space in: {}'.format(vb_nii))

    ####################################################################
    # Wrapper
    ####################################################################

    def wm(index, name, session, reference_maps, population_map, result,
            file_population_map, fnirt_prefix):

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

        lock = Lock(name, 'fsl4pop', file_population_map)
        df.ix[index, 'locked'] = True

        dfile = os.path.dirname(file_population_map)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        lock.save(file_population_map)

        ####################################################################
        # File names
        ####################################################################

        nb_nii = fnirt_prefix + nb_file_template
        vb_estimate_nii = fnirt_prefix + vb_estimate_template
        coefficients_vb = fnirt_prefix + cvb
        coefficients_nb = fnirt_prefix + cnb

        ####################################################################
        # Get NB from a result instance
        ####################################################################

        if diffeomorphism_nb == 'fit':

            if result is None:
                print('{}: No Result found'.format(name.name()))
                df.ix[index,'valid'] = False
                lock.conditional_unlock(df, index, verbose)
                return

            result.mask()
            nb = result.get_field('intercept', 'point')

        ####################################################################
        # Get NB from a session instance
        ####################################################################

        elif diffeomorphism_nb == 'scan_cycle':

            if (session is None) or (cycle is None):
                print('{}: No Session found or CYCLE defined'.format(name.name()))
                df.ix[index,'valid'] = False
                lock.conditional_unlock(df, index, verbose)
                return

            if reference_maps is None:
                print('{}: No ReferenceMaps found, continue…'.format(name.name()))
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
                                    print("""{}: Cycle {:d} marked as outlying, using fallback.""".format(
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
                name=name.name()+'-{:d}'.format(scan_cycle_to_use))

        else:
            print('{}: Diffeomorphism type not supported'.format(name.name()))
            return

        ####################################################################
        # Spline coefficients file
        ####################################################################

        if verbose:
            print('{}: Save subject reference image in NB to: {}'.format(
                name.name(), nb_nii))

        dfile = os.path.dirname(nb_nii)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        ni.save(image2nii(nb), nb_nii)

        if ignore_existing_warpcoef:
            warpcoef_file = fnirt_prefix + warpcoef_file_template
        else:
            if isfile(warpcoef_file_template):
                warpcoef_file = warpcoef_file_template
            else:
                warpcoef_file = fnirt_prefix + warpcoef_file_template

        if verbose:
            print('{}: Spline coefficients file: {}'.format(
                name.name(), warpcoef_file))

        if not isfile(warpcoef_file) or ignore_existing_warpcoef:
            if verbose:
                print('{}: Fit spline coefficients.'.format(
                    name.name()))

            dfile = os.path.dirname(warpcoef_file)
            if dfile and not isdir(dfile):
               os.makedirs(dfile)

            dfile = os.path.dirname(vb_estimate_nii)
            if dfile and not isdir(dfile):
               os.makedirs(dfile)

            status = fnirt(
                        warpcoef_file    = warpcoef_file,
                        nb_nii           = nb_nii,
                        vb_estimate_nii  = vb_estimate_nii,
                        vb_nii           = vb_nii,
                        nb_mask          = None,
                        vb_mask          = None,
                        config           = cmd_config,
                        cmd              = cmd_fnirt,
                        verbose          = verbose,
                        )

            if not status:
                df.ix[index,'valid'] = False
                print('{}: Unable to fit warp coefficients'.format(name.name()))
                lock.conditional_unlock(df, index, verbose, True)
                return
        else:
            if verbose:
                print('{}: Use existing spline coefficients.'.format(
                    name.name()))

        if verbose:
            print('{}: Parse diffeomorphism and create population map'.format(name.name()))

        diffeomorphism = splines2warp(
                    warpcoef_file      = warpcoef_file,
                    vb                 = study.vb,
                    vb_nii             = vb_nii,
                    nb_nii             = nb_nii,
                    name               = name,
                    new_diffeomorphism = new_diffeomorphism,
                    coefficients_vb    = coefficients_vb,
                    coefficients_nb    = coefficients_nb,
                    cmd                = cmd_std2imgcoord)

        try:
            vb_estimate = nii2image(ni.load(vb_estimate_nii),
                    name='vb_estimate_fnirt')
        except Exception as e:
            print('{}: Unable to read: {}, {}'.format(
                name.name(), vb_estimate_nii, e))
            vb_estimate = None

        if study.vb_background:
            vb_background = study.vb_background
        else:
            vb_background = None

        population_map = PopulationMap(diffeomorphism,
                vb=study.vb,
                nb=nb,
                vb_estimate=vb_estimate,
                vb_background = vb_background
                )

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

        if verbose > 2:
            print("""{}:
                {}
                {}""".format(name.name(),
                    population_map.diffeomorphism.describe(),
                    population_map.describe()))

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
                fnirt_prefix        = files['fnirt_prefix']
                pool.apply_async(wm, args=(index, name, session,
                    reference_maps, population_map, result,
                    file_population_map, fnirt_prefix))
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
                fnirt_prefix        = files['fnirt_prefix']
                wm(index, name, session, reference_maps, population_map,
                        result, file_population_map, fnirt_prefix)
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
