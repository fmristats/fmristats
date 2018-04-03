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

Quality assessment statistics for image quality, ability to track head
movements, and descriptive statistics of fitted effect fields.

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

    parser.add_argument('--reference-maps',
            default='../data/ref/{2}/{0}-{1:04d}-{2}-{3}.ref',
            help='input file;' + hp.reference_maps)

    parser.add_argument('--fit',
            default='../data/fit/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}-{5}.fit',
            help='input file;' + hp.sfit)

    parser.add_argument('--population-space',
            default='resolution',
            help=hp.population_space)

    parser.add_argument('--scale-type',
            default='max',
            choices=['diagonal','max','min'],
            help=hp.scale_type)

########################################################################
# Arguments specific for using the protocol API
########################################################################

    parser.add_argument('--protocol',
            help=hp.protocol)

    parser.add_argument('--cohort',
            help=hp.cohort)

    parser.add_argument('--id',
            type=int,
            nargs='+',
            help=hp.j)

    parser.add_argument('--datetime',
            help=hp.datetime)

    parser.add_argument('--paradigm',
            help=hp.paradigm)

    parser.add_argument('--strftime',
            default='%Y-%m-%d-%H%M',
            help=hp.strftime)

########################################################################
# Output arguments
########################################################################

    parser.add_argument('-o', '--protocol-log',
            default='logs/{}-fitassessment.pkl',
            help=hp.protocol_log)

    parser.add_argument('--assess-sessions',
            action='store_true',
            help="""whether to assess descriptive statistics of the
            sessions""")

    parser.add_argument('--assess-reference-maps',
            action='store_true',
            help="""whether to assess descriptive statistics of the
            reference maps""")

    parser.add_argument('--assess-fits',
            action='store_true',
            help="""whether to assess descriptive statistics of the
            fits""")

########################################################################
# Miscellaneous
########################################################################

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

from ...load import load_result, load_session, load_refmaps

from ...name import Identifier

from ...protocol import layout_dummy, layout_sdummy

#from ... import Session, ReferenceMaps, PopulationMap, SignalModel, Result

import pandas as pd

from pandas import Series

import datetime

import sys

import os

from os.path import isfile, isdir, join

import numpy as np

########################################################################

def call(args):
    output = args.protocol_log.format(
            datetime.datetime.now().strftime('%Y-%m-%d-%H%M'))

    ####################################################################
    # Parse protocol
    ####################################################################

    df = get_df(args)

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

    layout_dummy(df_layout, 'ref',
            template=args.reference_maps,
            strftime=args.strftime
            )

    layout_sdummy(df_layout, 'fit',
            template=args.fit,
            urname=args.population_space,
            scale_type=args.scale_type,
            strftime=args.strftime
            )

    ####################################################################
    # Add fields for summary statistics
    ####################################################################

    if args.assess_sessions:
        df['muf'] = Series(np.nan, index=df.index, dtype=float)
        df['mub'] = Series(np.nan, index=df.index, dtype=float)
        df['mfb'] = Series(np.nan, index=df.index, dtype=float)
        df['sdf'] = Series(np.nan, index=df.index, dtype=float)
        df['sdb'] = Series(np.nan, index=df.index, dtype=float)
        df['cvf'] = Series(np.nan, index=df.index, dtype=float)
        df['cvb'] = Series(np.nan, index=df.index, dtype=float)

    if args.assess_reference_maps:
        df['cyc']= Series(0, index=df.index, dtype=int)
        df['scn']= Series(0, index=df.index, dtype=int)
        df['cyo']= Series(0, index=df, dtype=int)
        df['sco']= Series(0, index=df, dtype=int)

    if args.assess_fits:
        df['vol'] = Series(np.nan, index=df, dtype=float)
        df['occ'] = Series(0, index=df, dtype=int)

    #################################################################
    # Apply wrappers
    #################################################################

    if args.assess_sessions:
        for r in df_layout.itertuples():
            name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)
            wrapper_session(
                    name     = name,
                    df       = df,
                    index    = r.Index,
                    filename = r.ses,
                    verbose  = args.verbose
                    )

    if args.assess_reference_maps:
        for r in df_layout.itertuples():
            name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)
            wrapper_refmaps(
                    name     = name,
                    df       = df,
                    index    = r.Index,
                    filename = r.ref,
                    verbose  = args.verbose
                    )

    if args.assess_fits:
        for r in df_layout.itertuples():
            name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)
            wrapper_fits(
                    name     = name,
                    df       = df,
                    index    = r.Index,
                    filename = r.fit,
                    urname   = args.population_space,
                    verbose  = args.verbose
                    )

    ####################################################################
    # Write protocol
    ####################################################################

    if args.verbose:
        print('Save: {}'.format(output))

    dfile = os.path.dirname(output)
    if dfile and not isdir(dfile):
       os.makedirs(dfile)

    df.to_pickle(output)

#######################################################################

def wrapper_session(name, df, index, filename, verbose):
    session = load_session(filename, name, df, index, verbose)
    if df.ix[index,'valid'] == False:
        return

    dstats = session.descriptive_statistics()
    dstats_names = ['muf', 'mub', 'mfb', 'sdf', 'sdb', 'cvf', 'cvb']
    dstats = Series(dstats, index=dstats_names)
    for s in dstats_names:
        df.ix[index,s] = dstats[s]

def wrapper_refmaps(name, df, index, filename, verbose):
    reference_maps = load_refmaps(filename, name, df, index, verbose)
    if df.ix[index,'valid'] == False:
        return

    dstats = reference_maps.descriptive_statistics()
    dstats_names = ['cyc', 'scn', 'cyo', 'sco']
    dstats = Series(dstats, index=dstats_names)
    for s in dstats_names:
        df.ix[index,s] = dstats[s]

def wrapper_fits(name, df, index, filename, urname, verbose):
    result = load_result(filename, name, df, index, urname, verbose)
    if df.ix[index,'valid'] == False:
        return

    dstats = result.descriptive_statistics()
    dstats_names = ['vol', 'occ']
    dstats = Series(dstats, index=dstats_names)
    for s in dstats_names:
        df.ix[index,s] = dstats[s]
