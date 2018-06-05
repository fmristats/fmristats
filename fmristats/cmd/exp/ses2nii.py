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

Export session file to Nifti1

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

    parser.add_argument('--session',
            default='../data/ses/{2}/{0}-{1:04d}-{2}-{3}.ses',
            help='input file;' + hp.session)

########################################################################
# Output arguments
########################################################################

    parser.add_argument('--nii',
            default='../results/nii/{2}/{0}-{1:04d}-{2}-{3}-{{:d}}.nii.gz',
            help="""output file""")

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-ses2nii.pkl',
            help=hp.protocol_log)

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
# Arguments which field to extract
########################################################################

    parser.add_argument('--cycle',
            type=int,
            nargs='+',
            help="""which cycle or cycles to extract.""")

########################################################################
# Miscellaneous
########################################################################

    parser.add_argument('-f', '--force',
            action='store_true',
            help=hp.force.format('result'))

    parser.add_argument('-v', '--verbose',
            action='store_true',
            help=hp.verbose)

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

from ...load import load_session

from ...name import Identifier

from ...protocol import layout_dummy, layout_sdummy

from ...diffeomorphisms import Image

from ...nifti import image2nii

import datetime

import sys

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

import nibabel as ni

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

    layout_dummy(df_layout, 'nii',
            template=args.nii,
            strftime=args.strftime
            )

    layout_dummy(df_layout, 'filename',
            template=args.session,
            strftime=args.strftime
            )

    ####################################################################
    # Apply wrapper
    ####################################################################

    def wm(r):
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)

        wrapper(name                  = name,
                df                    = df,
                index                 = r.Index,
                filename              = r.filename,
                template              = r.nii,
                verbose               = args.verbose,
                cycle                 = args.cycle,
                )

    it =  df_layout.itertuples()

    for r in it:
        wm(r)

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

def wrapper(name, df, index, filename, template, verbose,
        cycle):

    ####################################################################
    # Load fit from disk
    ####################################################################

    session = load_session(filename, name, df, index, verbose)
    if df.ix[index,'valid'] == False:
        return

    if verbose > 1:
        print('{}: Description of the session:'.format(session.name.name()))
        print(session.describe())

    for c in cycle:
        fname = template.format(c)
        dfile = os.path.dirname(fname)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        field = Image(
                reference=session.reference,
                data=session.raw[c],
                name=session.name.name()+'-{:d}'.format(c))

        if verbose:
            print('{}: Save: {} to: {}'.format(
                name.name(), field.name, fname))
        ni.save(image2nii(field), fname)
