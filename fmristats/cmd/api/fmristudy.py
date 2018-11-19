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

Handling studies

"""

import argparse

# TODO: --protocol_from_csv
# TODO: --covariates_from_csv

def add_study_arguments(parser):

    study_parser = parser.add_argument_group(
            """Create or manage a study""")

    study_parser.add_argument('--study',
        help="""A study.""")

    study_parser.add_argument('--single-subject',
        help="""If this is a study that only contains a single subject,
        you may change the default file layout to name all files
        starting with the string: SINGLE_SUBJECT.""")

    study_parser.add_argument('--root_directory',
        help="""Will overwrite the default (the work-directory) and set
        the root directory of the study to ROOT_DIRECTORY.""")

    study_parser.add_argument('--traverse-upwards',
        type=int,
        default=6,
        help="""When looking for a valid study file, traverse upwards
        TRAVERSE_UPWARDS number of directories. Note that the search
        will never traverse further upwards than your home
        directory.""")

    study_parser.add_argument('--protocol',
        help="""A protocol. The protocol file contains basic information
        about the FMRI sessions in a study.""")

    study_parser.add_argument('--covariates',
        help="""The covariates. The covariates file information about
        the subjects in the study.""")

    parser.add_argument('--protocol-update',
        nargs='+',
        help="""Update field 'valid' in the protocol file with the
        entries in the this file.""")

    parser.add_argument('--covariates-update',
        nargs='+',
        help="""Update field 'valid' in the protocol file with the
        entries in the this file.""")

    study_parser.add_argument('--protocol-query',
            help="""query the protocol.""")

    study_parser.add_argument('--covariates-query',
            help="""query the covariates.""")

    study_parser.add_argument('--stimulus',
        help="""Path to a stimulus file or template for such a
        file""")

    study_parser.add_argument('--session',
        help="""Path to a session file or template for such a file""")

    study_parser.add_argument('--epi-code',
        type=int,
        help="""This code defines the direction in which the planes in
        the FMRI have been measured. If the EPIs have been measured
        inferior to superior, then EP is 2, if they have been measured
        superior to inferior, then EP is -2. if they have been measured
        posterior to anterior, then EP is 1. Values in a protocol file
        will always take precedence.""")

    study_parser.add_argument('--reference-maps',
        help="""Path to a reference maps or template for such a file""")

    study_parser.add_argument('--population-map',
        help="""Path to a population map file or template for such a
        file""")

    study_parser.add_argument('--fit',
        help="""Path to a result file or template for such a file""")

    study_parser.add_argument('--strftime',
        help="""Format string of date and time. Convert time to string
        according to this format specification.""")

    study_parser.add_argument('--cohort',
        help = """Only use entries which belong to this cohort. If no
        study or protocol file is given, then this is the name of the
        cohort of the subject.""")

    study_parser.add_argument('--id',
        type=int,
        nargs='+',
        help = """Only process entries which belong to this subject ID.
        If no study or protocol file is given, then this is the id of
        the subject in its cohort.""")

    study_parser.add_argument('--datetime',
        help = """Only use entries which are stamped with this date and
        time. If no study or protocol file is given, then this is the
        date and time of the FMRI session of the subject.""")

    study_parser.add_argument('--paradigm',
        help = """Only use entries which belong to this stimulus design.
        The name of the stimulus will appear as part of the default file
        names produced by various fmristats command line tools. Unless
        you are planing to provide custom file names to all of
        fmristats' command line tools, try not to use any special
        characters in this name.""")

    study_parser.add_argument('--vb-image',
        help="""A template in standard space""")

    study_parser.add_argument('--vb-nii',
        help="""A template in standard space""")

    study_parser.add_argument('--vb-name',
        help="""The name of the standard space that shall be used in the
        analysis.""")

    study_parser.add_argument('--vb-background-image',
        help="""Image file with a background image in VB""")

    study_parser.add_argument('--vb-background-nii',
        help="""Image file with a background image in VB""")

    study_parser.add_argument('--vb-background-name',
        help="""Name of the background image in VB""")

    study_parser.add_argument('--vb-ati-image',
        help="""A template for the ATI reference field in standard
        space""")

    study_parser.add_argument('--vb-ati-nii',
        help="""A template for the ATI reference field in standard
        space""")

    study_parser.add_argument('--vb-ati-name',
        help="""A name for the template for the ATI reference field in
        standard space""")

    study_parser.add_argument('--diffeomorphism',
        help="""Name of the diffeomorphism that maps between standard
        space and subject reference space.""")

    study_parser.add_argument('--scale-type',
        choices=['diagonal','max','min'],
        help="""May be any string but special keywords are reserved:
        diagonal, max, and min. If SCALE_TYPE is `diagonal`, then SCALE
        (in fmrifit) is set to one half of the length of the diagonal of
        the orthorhombic measure lattice.  If SCALE_TYPE is `min` (or
        `max`), then SCALE is one half of the minimal (or maximal) edge
        length of the orthorhombic measure lattice.""")

    study_parser.add_argument('-o', '--out',
        help="""Save possibly modified study instance to OUT.""")

from ...epilog import epilog

def define_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=epilog)

    add_study_arguments(parser)

    ####################################################################
    # Verbosity
    ####################################################################

    control_verbosity  = parser.add_argument_group(
        """Control the level of verbosity""")

    control_verbosity.add_argument('-v', '--verbose',
        action='count',
        default=0,
        help="""Increase output verbosity""")

    return parser

def cmd():
    parser = define_parser()
    args = parser.parse_args()
    call(args)

cmd.__doc__ = __doc__

########################################################################
#
# Load libraries
#
########################################################################

from ... import load

from ...name import Identifier

from ...study import Study

from ...stimulus import Block

import pandas as pd

import datetime

import os

from os.path import isfile, isdir, join

from pathlib import Path

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

########################################################################

def get_study(args):

    # Parse template image in VB
    if args.vb_image:
        try:
            if args.verbose:
                print('Read (vb): {}'.format(args.vb_image))
            vb = load(args.vb_image)
            if args.vb_name is not None:
                vb.name = args.vb_name
        except Exception as e:
            vb = None
            print('Unable to read --vb-image: {}'.format(args.vb_image))
            print('Using fallback --vb-nii: {}'.format(e))
    elif args.vb_nii:
        try:
            if args.verbose:
                print('Read: {}'.format(args.vb_nii))
            vb = nii2image(ni.load(args.vb_nii), name=args.vb_name)
        except Exception as e:
            print('Unable to read --vb-nii: {}, {}'.format(args.vb_image, e))
            return
    else:
        vb = None

    # Parse background template image in VB
    if args.vb_background_image:
        try:
            if args.verbose:
                print('Read (vb-background): {}'.format(args.vb_background_image))
            vb_background = load(args.vb_background_image)
            if args.vb_background_name is not None:
                vb.name = args.vb_background_name
        except Exception as e:
            vb_background = None
            print('Unable to read --vb-background-image: {}'.format(args.vb_background_image))
            print('Using fallback --vb-background-nii: {}'.format(e))
    elif args.vb_background_nii:
        try:
            if args.verbose:
                print('Read: {}'.format(args.vb_background_nii))
            vb_background = nii2image(ni.load(args.vb_background_nii), name=args.vb_background_name)
        except Exception as e:
            print('Unable to read --vb-background-nii: {}, {}'.format(args.vb_background_image, e))
            return
    else:
        vb_background = None

    # Parse ATI reference template image in VB
    if args.vb_ati_image:
        try:
            if args.verbose:
                print('Read (vb-ati): {}'.format(args.vb_ati_image))
            vb_ati = load(args.vb_ati_image)
            if args.vb_ati_name is not None:
                vb.name = args.vb_ati_name
        except Exception as e:
            vb_ati = None
            print('Unable to read --vb-ati-image: {}'.format(args.vb_ati_image))
            print('Using fallback --vb-ati-nii: {}'.format(e))
    elif args.vb_ati_nii:
        try:
            if args.verbose:
                print('Read: {}'.format(args.vb_ati_nii))
            vb_ati = nii2image(ni.load(args.vb_ati_nii), name=args.vb_ati_name)
        except Exception as e:
            print('Unable to read --vb-ati-nii: {}, {}'.format(args.vb_ati_image, e))
            return
    else:
        vb_ati = None

    # Parse or create the protocol file
    if args.protocol:
        try:
            protocol = pd.read_pickle(args.protocol)
            if args.verbose:
                print('Read protocol: {}'.format(args.protocol))
        except Exception as e:
            print('Unable to read protocol file {}'.format(args.protocol))
            print('Exception: {}'.format(e))
            return
    else:
        protocol = None

    if protocol is None and args.id is not None:
        assert len(args.id) == 1, \
                """When no protocol is specified, only a single id is
                allowed"""

        if hasattr(args, 'epi_code'):
            epi_code = args.epi_code
        else:
            epi_code = None

        try:
            date = pd.to_datetime(args.datetime, format=args.strftime)
        except Exception as e:
            print('Unable to parse datetime string: {}'.format(args.datetime))
            print(e)
            return
        name = Identifier(cohort=args.cohort, j=args.id[0],
                datetime=date, paradigm=args.paradigm)
        protocol = name.to_data_frame(epi_code)

    force_single_subject = False
    if protocol is None:
        if hasattr(args, 'epi_code'):
            epi_code = args.epi_code
        else:
            epi_code = None

        for f in [args.session, args.reference_maps,
                args.population_map, args.fit, args.stimulus]:
            if f is not None:
                #print('Read: {}'.format(f))
                try:
                    instance = load(f)
                    protocol = instance.name.to_data_frame(epi_code)
                    print('Processing subject: {}'.format(
                        instance.name.name()))
                    force_single_subject = True
                    break
                except Exception as e:
                    pass
                    #print('…failed with {}'.format(e))

    # Parse or create the covariate file
    if args.covariates:
        try:
            covariates = pd.read_pickle(args.covariates)
            if args.verbose:
                print('Read covariates: {}'.format(args.covariates))
        except Exception as e:
            print('Unable to read covariates file {}'.format(args.covariates))
            print('Exception: {}'.format(e))
            return
    else:
        covariates = None

    if args.single_subject:
        if args.single_subject == 'yes':
            single_subject = True
        else:
            single_subject = args.single_subject
    elif force_single_subject:
        single_subject = True
    else:
        single_subject = False

    # Find study file
    if (args.study is None) and not single_subject:
        dir_home = str(Path.home())
        dir_curr = '.'
        dir_n = args.traverse_upwards

        for n in range(dir_n):
            if args.verbose > 1:
                print('Search for study file in: {}'.format(dir_curr))
            x = [f for f in os.listdir(dir_curr) if f.endswith(".study")]
            if x != []:
                break
            else:
                if os.path.abspath(dir_curr) == dir_home:
                    break
                else:
                    dir_curr = join('..', dir_curr)

        if len(x) > 1:
            print('Found more than one study: {}'.format(x))
            return
        elif len(x) == 1:
            study_dir = dir_curr
            study_file = x[0]
            print('Found study: {}'.format(os.path.join(study_dir, study_file)))

            os.chdir(study_dir)
            args.study = study_file
            print('Working directory has been changed to: {}'.format(os.getcwd()))

    # Parse or create the study file
    if args.study:
        try:
            if args.verbose:
                print('Read study: {}'.format(args.study))
            study = pd.read_pickle(args.study)
        except Exception as e:
            print('Unable to read study file {}'.format(args.study))
            print('Exception: {}'.format(e))
            return
    else:
        study = None

    # Exit if nothing to do
    if (study is None) and (protocol is None):
        print('Unable to find a valid study or protocol.')
        return

    file_layout = {
        'stimulus':args.stimulus,
        'session':args.session,
        'reference_maps':args.reference_maps,
        'result':args.fit,
        'population_map':args.population_map,
        }

    if (study is None) or single_subject:
        study = Study(
                protocol=protocol,
                covariates=covariates,
                vb = vb,
                vb_background = vb_background,
                vb_ati = vb_ati,
                diffeomorphism = args.diffeomorphism,
                scale_type = args.scale_type,
                file_layout=file_layout,
                root_dir = args.root_directory,
                strftime = args.strftime,
                single_subject=single_subject)

    else:
        study.update_layout(file_layout)
        if protocol is not None:
            if args.verbose > 1:
                print('Set (or reset) protocol file in study')
            study.protocol = protocol
        if covariates is not None:
            if args.verbose > 1:
                print('Set (or reset) covariates file in study')
            study.covariates = covariates
        if vb is not None:
            if args.verbose > 1:
                print('Set (or reset) template in standard space')
            study.vb = vb
        if vb_background is not None:
            if args.verbose > 1:
                print('Set (or reset) background in standard space')
            study.vb_background = vb_background
        if vb_ati is not None:
            if args.verbose > 1:
                print('Set (or reset) ATI reference field in standard space')
            study.vb_ati = vb_ati
        if args.scale_type is not None:
            if args.verbose > 1:
                print('Set (or reset) scale_type to: {}'.format(
                    args.scale_type))
            study.scale_type = args.scale_type
        if args.diffeomorphism is not None:
            if args.verbose > 1:
                print('Set (or reset) diffeomorphism to: {}'.format(
                    args.diffeomorphism))
            study.diffeomorphism = args.diffeomorphism
        if args.strftime is not None:
            if args.verbose > 1:
                print('Set (or reset) strftime to: {}'.format(
                    args.strftime))

    if args.protocol_update:
        for upfile in args.protocol_update:
            try:
                if args.verbose:
                    print('Update protocol with: {}'.format(upfile))
                df = pd.read_pickle(upfile)
            except Exception as e:
                print('Unable to read file {}, {}'.format(args.study, e))
                return
            study.update_protocol(df)

    if args.covariates_update:
        for upfile in args.covariates_update:
            try:
                if args.verbose:
                    print('Update covariates with: {}'.format(upfile))
                df = pd.read_pickle(upfile)
            except Exception as e:
                print('Unable to read file {}, {}'.format(args.study, e))
                return
            study.update_covariates(df)

    study.filter(cohort=args.cohort, j=args.id, paradigm=args.paradigm, inplace=True)

    if args.protocol_query:
        study.query(args.protocol_query, inplace=True)

    if args.covariates_query:
        study.query(args.covariates_query, inplace=True)

    return study

def call(args):

    study = get_study(args)

    if study is None:
        return

    if args.out is not None:
        if args.verbose:
            print('Save: {}'.format(args.out))

        dfile = os.path.dirname(args.out)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        study.save(args.out)
