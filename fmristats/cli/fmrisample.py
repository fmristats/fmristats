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

Create a sample of fitted statistic fields of FMRI signal model
for inference

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

    specific = parser.add_argument_group(
        """Create a sample of statistic fields""")

    specific.add_argument('sample',
            help="""Name of the output file""")

    specific.add_argument('-f', '--force',
            action='store_true',
            help="""Force re-writing the sample file.""")

    specific.add_argument('--mask',
        help="""Set the mask to use. If yes, true, apply or not given,
        both vb and vb_mask will apply. If no, false or ignore, neither
        mask will be applied. If vb, vb_background, vb_estimate or
        vb_mask, the respective mask will be applied.""")

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

from ..load import load

from ..lock import Lock

from ..study import Study

from ..smodel import Result

from ..sample import Sample

########################################################################

def call(args):

    if isfile(args.sample) and not args.force:
        print('Sample file already exists, use -f/--force to overwrite.')
        if args.verbose:
            print('Parse: {}'.format(args.sample))
        sample = load(args.sample)
        print('Description of standard space:')
        print(sample.vb.describe())
        print('Description of sample:')
        print(sample.describe())
        return
    else:
        sample_file = args.sample

    ###################################################################

    study = get_study(args)

    if study is None:
        print('Nothing to do.')
        return

    if study.vb is None:
        print("""
        You need to provide a template in standard space (for example by
        either setting a template in the study or by providing a
        template using --vb-image or --vb-nii).""")

    if study.vb_ati is None:
        print("""
        You need to provide an ATI reference field in standard space
        (for example by either setting the field in the study or by
        providing a field using --vb-ati-image or --vb-ati-nii).""")

    if (study.vb is None) or (study.vb_ati is None):
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
    # Options
    ####################################################################

    force   = args.force
    verbose = args.verbose

    ####################################################################
    # Create the iterator
    ####################################################################

    study_iterator = study.iterate('result',
            verbose=args.verbose,
            integer_index=True)

    df = study_iterator.df.copy()

    ####################################################################
    ####################################################################

    if args.verbose:
        print('Create population sample: {}'.format(sample_file))

    statistics = np.empty(study.vb.shape + (3,len(df),))
    statistics [ ... ] = np.nan

    for index, name, instances in study_iterator:
        result = instances['result']
        if result is not None:
            if type(result) is Lock:
                print('{}: Result file is locked. Still fitting?'.format(
                    name.name()))
            elif type(result) is not Result:
                print('{}: File does not contain Result. Skipping.'.format(
                    name.name()))
            else:
                if verbose > 2:
                    print("""{}: Description of the fit:
                        {}
                        {}""".format(name.name(), result.describe(),
                            result.population_map.describe()))

                result.mask(mask=mask, verbose=verbose)
                intercept = result.get_field('intercept','point')
                c = study.vb_ati.data / intercept.data

                beta = result.get_field('task','point')
                beta_stderr = result.get_field('task','stderr')

                statistics[...,0,index] = c*beta.data
                statistics[...,1,index] = c*beta_stderr.data
                statistics[...,2,index] = c
        else:
            study_iterator.df.ix[index,'valid'] = False

    df = study_iterator.df.copy()
    del df['result']

    sample = Sample(
            covariates = df,
            statistics = statistics,
            study      = study)

    sample = sample.filter()

    if args.verbose > 1:
        print('Description of the population space:')
        print(sample.vb.describe())
        print('Description of the sample:')
        print(sample.describe())

    try:
        if verbose:
            print('Save: {}'.format(sample_file))

        dfile = os.path.dirname(sample_file)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        sample.save(sample_file)
    except Exception as e:
        print('{}: Unable to create: {}, {}'.format(name.name(),
            sample_file, e))

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
