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

Query covariate files

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

    parser.add_argument('protocol',
            help=hp.protocol)

########################################################################
# Output arguments
########################################################################

    parser.add_argument('-o', '--output',
            help="""output file""")

########################################################################
# Arguments specific for using the protocol API
########################################################################

    #to_process = parser.add_argument_group(
    #        """specifying the covariate entries to process""",
    #        """Arguments which give control which entries to
    #        process.""")

    #to_process.add_argument('--cohort',
    #        help=hp.cohort)

    #to_process.add_argument('--id',
    #        type=int,
    #        nargs='+',
    #        help=hp.j)

########################################################################
# Arguments specific for using fmriprotocol
########################################################################

    parser.add_argument('--covariate',
            help="""update field 'valid' in the protocol file such that
            only entries which are valid in the protocol *and* in
            COVARIATE file are marked as valid entries.""")

    parser.add_argument('--update',
            help="""update field 'valid' in the protocol file with the
            entries in the file UPDATE.""")

    parser.add_argument('--merge',
            help="""merge with covariate file MERGE.""")

    parser.add_argument('--query',
            help="""apply a QUERY to the data frame.""")

    parser.add_argument('--head',
            type=int,
            help="""print the first HEAD number of lines""")

    parser.add_argument('--to-csv',
            help="""export protocol file to CSV""")

########################################################################
# Miscellaneous
########################################################################

    misc = parser.add_argument_group(
            """miscellaneous""")

    misc.add_argument('-f', '--force',
            action='store_true',
            help=hp.force.format('protocol'))

    misc.add_argument('-v', '--verbose',
            action='count',
            default=0,
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

from ..df import filter_df

from ...name import Identifier

import pandas as pd

import sys

import os

from os.path import isfile, isdir, join

########################################################################

def call(args):
    try:
        df = pd.read_pickle(args.protocol)
        if args.verbose:
            print('Read: {}'.format(args.protocol))
    except Exception as e:
        print('Unable to read protocol file {}'.format(args.protocol))
        print('Exception: {}'.format(e))
        sys.exit()

    if args.verbose:
        valid = df.groupby(['cohort','valid']).id.agg(['count'])
        print('Number of entries in the list: {}'.format(len(df)))
        print("""\nOverview""")
        print(valid)
        print()

    if args.query:
        df = df.query(args.query).copy()

        if args.verbose:
            valid = df.groupby(['cohort','valid']).id.agg(['count'])
            print('Number of entries in the list after --query: {}'.format(len(df)))
            print("""\nOverview""")
            print(valid)
            print()

    if args.head:
        if args.verbose:
            print()
            print("""-------------------------------------------------------------""".format(args.head))
        print(df.head(args.head))
        print()

    if args.output:
        if not isfile(args.output) or args.force:
            if args.verbose:
                print('Save: {}'.format(args.output))
            df.to_pickle(args.output)
        else:
            print('File already exists: {}, use --force to overwrite'.format(args.output))

    if args.to_csv:
        if not isfile(args.to_csv) or args.force:
            if args.verbose:
                print('Save: {}'.format(args.to_csv))
            df.to_csv(args.to_csv, index=True)
        else:
            print('File already exists: {}, use --force to overwrite'.format(args.to_csv))
