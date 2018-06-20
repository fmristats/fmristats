# Copyright 2016-2017 Thomas W. D. Möbius
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

On the one hand, this tool is a wrapper to ``std2imgcoord`` and converts
a FSL warp coefficient file produced by ``fnirt`` and produces
population map for the corresponding session.  On the other hand, it can
be used as a wrapper to the FSL command line tools ``fnirt`` to estimate
former warp coefficient file from the intercept field of a signal model
and a given template brain.  The resulting population diffeomorphism is
save to disk.

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

    parser.add_argument('--template',
            default='/usr/share/data/fsl-mni152-templates/MNI152_T1_2mm.nii.gz',
            help=hp.template)

    parser.add_argument('--template-mask',
            default='/usr/share/data/fsl-mni152-templates/MNI152_T1_2mm_brain.nii.gz',
            help=hp.template_mask)

########################################################################
# Output arguments
########################################################################

    parser.add_argument('--population-space',
            default='fnirt',
            help='output name;' + hp.population_space_name)

    parser.add_argument('--image',
            default='../data/vbs/{}.image',
            help='output name;' + hp.population_space_directory)

    parser.add_argument('--population-map',
            default='../data/pop/{2}/{4}/{0}-{1:04d}-{2}-{3}-{4}.pop',
            help='output file;' + hp.population_map)

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-fsl4pop.pkl',
            help=hp.protocol_log)

########################################################################
# Additional input arguments when warp coefficient files provided
########################################################################

    parser.add_argument('--reference-space',
            default='warping/{0}-{1:04d}-{2}-{3}-reference-space.nii.gz',
            help='reference image used for the subject.')

    parser.add_argument('--warpcoef',
            default='warping/{0}-{1:04d}-{2}-{3}-warpcoef.nii.gz',
            help='warp coefficient file created by fnirt.')

    parser.add_argument('--ignore-existing-warpcoef',
            action='store_true',
            help="""by default ``fsl4pop`` will use the warp coefficient
            files it finds on disk, as it will assume that you will have
            created them to fit your needs. When this flag is set,
            ``fsl4pop`` will re-fit the warp coefficients using FNIRT
            using some reasonable defaults.""")

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
            help="""if warp coefficient files are provided or orking
            with session files, you don't need this""" + hp.sfit)

    parser.add_argument('--scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help="""only needed if part of the template for --nb-fit""" + hp.scale_type)

    parser.add_argument('--nb',
            default='reference',
            help="""name of the population space that was originally
            used for the fit"""  + hp.population_space)

########################################################################
# (optional) Controlling the behaviour of FNIRT
########################################################################

    parser.add_argument('--cmd-fnirt',
            default='fsl5.0-fnirt',
            help="""FSL fnirt command. Must be in your path.""")

    parser.add_argument('--cmd-config',
            default='T1_2_MNI152_2mm',
            help="""Config file for FSL FNIRT.""")

    parser.add_argument('--cmd-std2imgcoord',
            default='fsl5.0-std2imgcoord',
            help="""FSL std2imgcoord command. Must be in your path.""")

    parser.add_argument('--cmd-bet',
            default='fsl5.0-bet',
            help="""FSL bet command. Must be in your path.""")

    parser.add_argument('--variante',
            default='R',
            help="""FSL bet command. Must be in your path.""")

    # parser.add_argument('--provide-fnirt-with-templates',
    #         action='store_true',
    #         help="""by default FNIRT will use the reference in
    #         CMD_CONFIG. If this flag is set, the TEMPLATE will replace
    #         the default defined therin.""")

########################################################################
# (optional) Intermediate or temporary files
########################################################################

    parser.add_argument('--preimage',
            default='warping/{0}-{1:04d}-{2}-{3}-preimage.nii.gz',
            help='warped INTERCEPT image estimated by fnirt.')

    parser.add_argument('--ccycle',
            default='warping/{0}-{1:04d}-{2}-{3}-cycle.nii.gz',
            help="""coordinates of the template's image grid in
            population space""")

    parser.add_argument('--cpopulation',
            default='warping/{0}-{1:04d}-{2}-{3}-cpopulation.txt',
            help="""coordinates of the template's image grid in
            population space""")

    parser.add_argument('--csubject',
            default='warping/{0}-{1:04d}-{2}-{3}-csubject.txt',
            help="""coordinates of the template's image grid in subject
            space""")

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

from ...fsl import fit_warpcoef, warpcoef2pmap, bet

import pandas as pd

import datetime

import sys

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

import nibabel as ni

########################################################################

def call(args):
    try:
        template = ni.load(args.template)
    except Exception as e:
        print('Unable to read template file: {}'.format(args.template))
        print('Exception: {}'.format(e))
        exit()

    template = nii2image(nii=template, name=args.population_space)
    fname = args.image.format(args.population_space)
    dfile = os.path.dirname(fname)
    if dfile and not isdir(dfile):
       os.makedirs(dfile)

    template.save(fname)

    try:
        template_mask = ni.load(args.template_mask)
    except Exception as e:
        print('Unable to read template_mask file: {}'.format(args.template_mask))
        print('Exception: {}'.format(e))
        exit()

    template_mask = nii2image(nii=template_mask,
            name=args.population_space)
    fname = args.image.format(args.population_space+'-mask')
    dfile = os.path.dirname(fname)
    if dfile and not isdir(dfile):
       os.makedirs(dfile)

    template_mask.save(fname)

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
    # Add file layout
    ####################################################################

    df_layout = df.copy()

    layout_dummy(df_layout, 'ses',
            template=args.session,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'fit',
            template=args.fit,
            urname=args.nb,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'warpcoef',
            template=args.warpcoef,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'intercept',
            template=args.reference_space,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'preimage',
            template=args.preimage,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'ccycle',
            template=args.ccycle,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'cpopulation',
            template=args.cpopulation,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'csubject',
            template=args.csubject,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'filename',
            template=args.population_map,
            urname=args.population_space,
            scale_type=None,
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

        try:
            dfile = os.path.dirname(r.warpcoef)
            if dfile and not isdir(dfile):
                os.makedirs(dfile)
        except Exception as e:
            print('{}: {}'.format(name.name(), e))

        try:
            dfile = os.path.dirname(r.intercept)
            if dfile and not isdir(dfile):
                os.makedirs(dfile)
        except Exception as e:
            print('{}: {}'.format(name.name(), e))

        try:
            dfile = os.path.dirname(r.preimage)
            if dfile and not isdir(dfile):
                os.makedirs(dfile)
        except Exception as e:
            print('{}: {}'.format(name.name(), e))

        try:
            dfile = os.path.dirname(r.ccycle)
            if dfile and not isdir(dfile):
                os.makedirs(dfile)
        except Exception as e:
            print('{}: {}'.format(name.name(), e))

        try:
            dfile = os.path.dirname(r.cpopulation)
            if dfile and not isdir(dfile):
                os.makedirs(dfile)
        except Exception as e:
            print('{}: {}'.format(name.name(), e))

        try:
            dfile = os.path.dirname(r.csubject)
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
                file              = r.filename,

                file_res          = r.fit,
                file_ses          = r.ses,

                warpcoef_file     = r.warpcoef,
                intercept_file    = r.intercept,
                preimage_file     = r.preimage,
                ccycle            = r.ccycle,
                cpopulation       = r.cpopulation,
                csubject          = r.csubject,

                vb                = args.population_space,
                nb                = args.nb,
                vb_file           = args.template,
                vb_mask_file      = args.template_mask,
                template          = template,
                template_mask     = template_mask,

                cycle             = args.cycle,
                ignore_existing_warpcoef = args.ignore_existing_warpcoef,
                #provide_fnirt_with_templates = args.provide_fnirt_with_templates,

                cmd_fnirt        = args.cmd_fnirt,
                cmd_config       = args.cmd_config,
                cmd_std2imgcoord = args.cmd_std2imgcoord,
                cmd_bet          = args.cmd_bet,
                variante         = args.variante
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
        print('Save: {}'.format(output))

    dfile = os.path.dirname(output)
    if dfile and not isdir(dfile):
       os.makedirs(dfile)

    df.to_pickle(output)

########################################################################

def wrapper(name, df, index, remove_lock, ignore_lock, force, skip,
        verbose, file, file_res, file_ses, warpcoef_file, intercept_file,
        preimage_file, ccycle, cpopulation, csubject, vb, nb,
        vb_file, vb_mask_file, template, template_mask, cycle,
        ignore_existing_warpcoef, cmd_fnirt, cmd_config,
        cmd_std2imgcoord, cmd_bet, variante):

    if isfile(file):
        instance = load_population_map(file, name, df, index, vb, verbose)
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

    lock = Lock(name, 'fsl4pop', file)
    df.ix[index, 'locked'] = True
    lock.save(file)
    df.ix[index,'valid'] = True

    ####################################################################
    # Load or create warp coefficient file
    ####################################################################

    if not isfile(warpcoef_file) or ignore_existing_warpcoef:

        if cycle is None:
            result = load_result(file_res, name, df, index, nb, verbose)
            if lock.conditional_unlock(df, index, verbose):
                return

            result.mask()
            reference_space = result.get_field('intercept', 'point').round()
            intercept = image2nii(reference_space)

            try:
                intercept.to_filename(intercept_file)
            except Exception as e:
                df.ix[index,'valid'] = False
                print('{}: Unable to write: {}'.format(name.name(), intercept_file))
                print('{}: Exception: {}'.format(name.name(), e))
                lock.conditional_unlock(df, index, verbose, True)

        else:
            session = load_session(file_ses, name, df, index, verbose)
            if lock.conditional_unlock(df, index, verbose):
                return

            reference_space = Image(
                reference=session.reference,
                data=session.raw[cycle],
                name=session.name.name()+'-{:d}'.format(cycle))

            reference_space = bet(
                    intercept = reference_space,
                    intercept_file = ccycle,
                    mask_file = intercept_file,
                    cmd = cmd_bet,
                    variante = variante,
                    verbose = verbose)

            if reference_space is None:
                df.ix[index,'valid'] = False
                print('{}: Unable to bet'.format(name.name()))
                lock.conditional_unlock(df, index, verbose, True)
                return

        #if provide_fnirt_with_templates:
        #    if verbose:
        #        print('{}: … nb-templates (and nb-masks) defined by user'.format(
        #            name.name()))
        #else:
        #    if verbose:
        #        print('{}: … nb-templates (and nb-masks) as defined in: {}'.format(
        #            name.name(), cmd_config))
        #    vb_file = None
        #    vb_mask = None

        if verbose:
            print('{}: Reference space: {}'.format(name.name(), reference_space.name))
            print('{}: {}'.format(name.name(), intercept_file))
            print('{}: Fit warp coefficients...'.format(name.name()))

        status = fit_warpcoef(
                    nb_file        = intercept_file,
                    vb_file        = None, #vb_file,
                    vb_mask        = None, #vb_mask,
                    warpcoef_file  = warpcoef_file,
                    preimage_file  = preimage_file,
                    config         = cmd_config,
                    cmd            = cmd_fnirt,
                    verbose        = verbose,
                    )

        if not status:
            df.ix[index,'valid'] = False
            print('{}: Unable to fit warp coefficients'.format(name.name()))
            lock.conditional_unlock(df, index, verbose, True)
            return
    else:
        if verbose:
            print('{}: Use existing warp coefficients'.format(name.name()))

    ####################################################################
    # Create population map from warp coefficients file
    ####################################################################

    if verbose:
        print('{}: Create population map'.format(name.name()))

    population_map = warpcoef2pmap(
                warpcoef_file    = warpcoef_file,
                vb_file          = vb_file,
                nb_file          = intercept_file,
                vb_name          = vb,
                nb_name          = name,
                name             = 'fnirt',  # TODO: should be user choice
                vb               = 'mni152', # TODO: should be user choice
                cpopulation_file = cpopulation,
                csubject_file    = csubject,
                cmd              = cmd_std2imgcoord)

    population_map.diffeomorphism.metadata['vb_mask'] = vb_mask_file
    population_map.set_template_mask(template_mask)

    if population_map is None:
        df.ix[index,'valid'] = False
        print('{}: Unable to create: {}'.format(name.name(), file))
        print('{}: Exception: {}'.format(name.name(), e))
        lock.conditional_unlock(df, index, verbose, True)
        return

    if isfile(preimage_file):
        try:
            preimage = ni.load(preimage_file)
            population_map.set_target(nii2image(preimage, name='fnirt'))
        except Exception as e:
            df.ix[index,'valid'] = False
            print('{}: Unable to read: {}'.format(name.name(), preimage_file))
            print('{}: Exception: {}'.format(name.name(), e))

    if verbose:
        print('{}: Save: {}'.format(name.name(), file))

    try:
        population_map.save(file)
        df.ix[index,'locked'] = False
    except Exception as e:
        df.ix[index,'valid'] = False
        print('{}: Unable to write: {}'.format(name.name(), file))
        print('{}: Exception: {}'.format(name.name(), e))
        lock.conditional_unlock(df, index, verbose, True)
