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

Converts a csv to a covariate file

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

import pandas as pd

########################################################################

def call(args):
    df = pd.read_csv(args.csv,
            dtype={
                'cohort':str,
                'id':int,
                'valid':bool})

    df.set_index(keys=['cohort', 'id'],
            inplace=True,
            verify_integrity=True,
            drop=False)

    df['cohort'] = df.cohort.astype('category')

    if not hasattr(df, 'valid'):
        df['valid'] = True

    df.to_pickle(args.output)
