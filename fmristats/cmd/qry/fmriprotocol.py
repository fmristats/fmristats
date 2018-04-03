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

Query protocol files

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

    to_process = parser.add_argument_group(
            """specifying the protocol entries to process""",
            """Arguments which give control which protocol entries to
            process.""")

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
# Arguments specific for using fmriprotocol
########################################################################

    parser.add_argument('--update',
            nargs='+',
            help="""update field 'valid' in the protocol file with the
            entries in the file UPDATE. This step will
            be applied *after* applying your filters.""")

    parser.add_argument('--query',
            help="""apply a QUERY to the data frame.""")

    parser.add_argument('--head',
            type=int,
            help="""print the first HEAD number of lines to standard
            output.""")

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
    if args.strftime == 'short':
        args.strftime = '%Y-%m-%d'

    try:
        df = pd.read_pickle(args.protocol)
        if args.verbose:
            print('Read: {}'.format(args.protocol))
    except Exception as e:
        print('Unable to read protocol file {}'.format(args.protocol))
        print('Exception: {}'.format(e))
        sys.exit()

    if args.verbose:
        valid = df.groupby(['cohort','paradigm','valid']).epi.agg(['count'])
        print('Number of entries in the list: {}'.format(len(df)))
        if args.verbose > 1:
            print()
            print(valid)
            print()

    df = filter_df(df, args.cohort, args.id, args.paradigm)

    if args.verbose:
        valid = df.groupby(['cohort','paradigm','valid']).epi.agg(['count'])
        print('Number of relevant entries in the list: {}'.format(len(df)))
        if args.verbose > 1:
            print()
            print(valid)
            print()

    if args.update:
        for upfile in args.update:
            up = pd.read_pickle(upfile)

            if 'date' in up.index.names:
                on = ['cohort','id','paradigm','date']
            elif 'paradigm' in up.index.names:
                on = ['cohort','id','paradigm']
            else:
                on = ['cohort','id']

            df = (df.join(up, on=on, lsuffix='_')
                    .assign(valid=lambda x: x.valid.fillna(True))
                    .assign(valid=lambda x: x.valid & x.valid_)
                    .drop('valid_', axis=1))

        if args.verbose:
            valid = df.groupby(['cohort','paradigm','valid']).epi.agg(['count'])
            print('Number of entries in the list after --update: {}'.format(len(df)))
            if args.verbose > 1:
                print()
                print(valid)
                print()
            if args.verbose > 2:
                print(df.head(args.head))
                print()

    if args.query:
        df = df.query(args.query).copy()

        if args.verbose:
            valid = df.groupby(['cohort','paradigm','valid']).epi.agg(['count'])
            print('Number of entries in the list after --query: {}'.format(len(df)))
            if args.verbose > 1:
                print()
                print(valid)
                print()

    df = filter_df(df, args.cohort, args.id, args.paradigm)

    if df is None:
        print('Number of relevant entries in the list: {}'.format(len(df)))
    else:
        if args.verbose:
            valid = df.groupby(['cohort','paradigm','valid']).epi.agg(['count'])
            print('Number of relevant entries in the list: {}'.format(len(df)))
            if args.verbose > 1:
                print()
                print(valid)
                print()

    if df is None:
        sys.exit()

    df = filter_df(df, args.cohort, args.id, args.paradigm)

    if args.head:
        if args.verbose:
            print()
            print("""First {:d} valid entries in the final protocol:
---------------------------------------------""".format(args.head))
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
            df.to_csv(args.to_csv, index=False)
        else:
            print('File already exists: {}, use --force to overwrite'.format(args.to_csv))
