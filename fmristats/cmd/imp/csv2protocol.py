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

Converts a csv to protocol

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

    parser.add_argument('csv', help="""csv file""")

    parser.add_argument('output', help="""output file""")

    parser.add_argument('--strftime',
            default='%Y-%m-%d-%H%M',
            help=hp.strftime)

    return parser

def cmd_pkl():
    parser = create_argument_parser()
    args = parser.parse_args()
    call(args, False)

def cmd_cov():
    parser = create_argument_parser()
    args = parser.parse_args()
    call(args, True)

cmd_pkl.__doc__ = __doc__
cmd_cov.__doc__ = __doc__

########################################################################
#
# Load libraries
#
########################################################################

import pandas as pd

########################################################################

def call(args, drop):
    df = pd.read_csv(args.csv,
            dtype={
                'cohort':str,
                'id':int,
                'paradigm':str,
                'date':str,
                'epi':int,
                'valid':bool})

    if args.strftime == 'short':
        args.strftime = '%Y-%m-%d'

    if not hasattr(df, 'valid'):
        df['valid'] = True

    if 'date' in df.columns:
        df['date'] = pd.to_datetime(df.date, format=args.strftime)
        df.set_index(
                keys=['cohort', 'id', 'paradigm', 'date'],
                drop=drop,
                inplace=True,
                verify_integrity=True)
        df.sort_index(inplace=True)
        if not drop:
            df.cohort = df.cohort.astype('category')
            df.paradigm = df.paradigm.astype('category')

    elif 'paradigm' in df.columns:
        df.set_index(
                keys=['cohort', 'id', 'paradigm'],
                drop=drop,
                inplace=True,
                verify_integrity=True)
        df.sort_index(inplace=True)
        if not drop:
            df.cohort = df.cohort.astype('category')
            df.paradigm = df.paradigm.astype('category')

    else:
        df.set_index(
                keys=['cohort', 'id'],
                drop=drop,
                inplace=True,
                verify_integrity=True)
        df.sort_index(inplace=True)
        if not drop:
            df.cohort = df.cohort.astype('category')

    df.to_pickle(args.output)
