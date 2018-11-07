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

Handling studies

"""

import fmristats.cmd.hp as hp

import argparse

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
        help="""Will overwrite the default and set the root directory of
        the study to ROOT_DIRECTORY.""")

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

    study_parser.add_argument('--protocol-query',
            help="""query the protocol.""")

    study_parser.add_argument('--covariates-query',
            help="""query the covariates.""")

    study_parser.add_argument('--stimulus',
        help="""Path to a stimulus file or template for such a
        file""")

    study_parser.add_argument('--session',
        help="""Path to a session file or template for such a file""")

    study_parser.add_argument('--reference-maps',
        help="""Path to a reference maps or template for such a file""")

    study_parser.add_argument('--population-map',
        help="""Path to a population map file or template for such a
        file""")

    study_parser.add_argument('--fit',
        help="""Path to a result file or template for such a file""")

    study_parser.add_argument('--strftime',
        help="""Format of date and time""")

    study_parser.add_argument('--cohort',
        help="""Cohort""")

    study_parser.add_argument('--id',
        type=int,
        nargs='+',
        help="""id""")

    study_parser.add_argument('--datetime',
        help="""Datetime""")

    study_parser.add_argument('--paradigm',
        help="""Paradigm""")

    study_parser.add_argument('--vb-name',
        default='self',
        help="""Name of the population space""")

    study_parser.add_argument('--diffeomorphism-name',
        default='identity',
        help="""Name of the diffeomorphism between population space and
        subject space""")

    study_parser.add_argument('--scale-type',
        default='max',
        choices=['diagonal','max','min'],
        help="""Scale type""")

    study_parser.add_argument('-o', '--out',
            help="""Save possibly modified study instance to OUT""")

def define_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=hp.epilog)

    add_study_arguments(parser)

    control_verbosity  = parser.add_argument_group(
            """Control the level of verbosity""")

    control_verbosity.add_argument('-v', '--verbose',
            action='count',
            default=0,
            help=hp.verbose)

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

import sys

import os

from os.path import isfile, isdir, join

from pathlib import Path

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

########################################################################

def get_study(args):

    # Parse or create the protocol file
    if args.protocol:
        try:
            protocol = pd.read_pickle(args.protocol)
            if args.verbose:
                print('Read: {}'.format(args.protocol))
        except Exception as e:
            print('Unable to read protocol file {}'.format(args.protocol))
            print('Exception: {}'.format(e))
            return
    else:
        protocol = None

    if protocol is None and args.id is not None:
        if hasattr(args, 'epi_code'):
            epi_code = args.epi_code
        else:
            epi_code = None

        assert len(args.id) == 1, \
                """When no protocol is specified, only a single id is
                allowed"""
        try:
            date = pd.to_datetime(args.datetime, format=args.strftime)
        except Exception as e:
            print('Unable to parse datetime string: {}'.format(args.datetime))
            print(e)
            return
        name = Identifier(cohort=args.cohort, j=args.id[0],
                datetime=date, paradigm=args.paradigm)
        protocol = name.to_data_frame(epi_code)

    if protocol is None:
        if hasattr(args, 'epi_code'):
            epi_code = args.epi_code
        else:
            epi_code = None

        for f in [args.session, args.reference_maps,
                args.population_map, args.fit, args.stimulus]:
            if f is not None:
                print('Read: {}'.format(f))
                try:
                    instance = load(f)
                    protocol = instance.name.to_data_frame(epi_code)
                except Exception as e:
                    print('…failed with {}'.format(e))

    # Parse or create the covariate file
    if args.covariates:
        try:
            covariates = pd.read_pickle(args.covariates)
            if args.verbose:
                print('Read: {}'.format(args.covariates))
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
            study = pd.read_pickle(args.study)
            if args.verbose:
                print('Read: {}'.format(args.study))
        except Exception as e:
            print('Unable to read study file {}'.format(args.study))
            print('Exception: {}'.format(e))
            return
    else:
        study = None

    # Exit if nothing to do
    if (study is None) and (protocol is None):
        print('Unable to find a valid protocol.')
        return

    file_layout = {
        'stimulus':args.stimulus,
        'session':args.session,
        'reference_maps':args.reference_maps,
        'result':args.fit,
        'population_map':args.population_map,
        'diffeomorphism_name':args.diffeomorphism_name,
        'scale_type':args.scale_type,
        }

    if (study is None) or single_subject:
        study = Study(protocol=protocol, covariates=covariates,
                file_layout=file_layout, strftime=args.strftime,
                single_subject=single_subject)

    else:
        study.update_layout(file_layout)
        if protocol is not None:
            study.protocol = protocol
        if covariates is not None:
            study.covariates = covariates

    study.filter(cohort=args.cohort, j=args.id, paradigm=args.paradigm, inplace=True)

    if args.protocol_query:
        study.query(args.protocol_query, inplace=True)

    if args.covariates_query:
        study.query(args.covariates_query, inplace=True)

    return study

def call(args):

    study = get_study(args)

    if study is None:
        sys.exit()

    if args.verbose > 1:
        #print(study.protocol.head())

        #if study.covariates is not None:
        #    print(study.covariates.head())

        for k, v in study.file_layout.items():
            print('{:<24}: {}'.format(k,v))

    if args.out is not None:
        if args.verbose:
            print('Save: {}'.format(args.out))

        dfile = os.path.dirname(args.out)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        study.save(args.out)
