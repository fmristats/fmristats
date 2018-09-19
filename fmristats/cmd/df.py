# Copyright 2018 Thomas W. D. Möbius
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

from ..load import load

from ..protocol import filter_protocol

import pandas as pd

import numpy as np

def filter_df(df, cohort, j, paradigm):

    if cohort is None:
        cohort = slice(None)

    if paradigm is None:
        paradigm = slice(None)

    if j is None:
        j = slice(None)
    elif len(j) == 1:
        j = j[0]
    elif len(j) == 2:
        j = slice(j[0], j[1])

    df = df.loc(axis=0)[cohort, j, paradigm]

    if len(df) < 1:
        return

    df = df[ df.valid == True].copy()
    df['epi'] = df.epi.astype(int)
    df.sort_index(inplace=True)
    return df

def get_df(args, fall_back=None):
    if args.protocol:
        try:
            df = pd.read_pickle(args.protocol)
            if args.verbose:
                print('Read: {}'.format(args.protocol))
        except Exception as e:
            print('Unable to read protocol file {}'.format(args.protocol))
            print('Exception: {}'.format(e))
            return

        return filter_df(df, args.cohort, args.id, args.paradigm)

    elif args.id is not None:
        if len(args.id) > 1:
            print('When no protocol file is specified, then only one id is allowed')
            return

        try:
            date = pd.to_datetime(args.datetime, format=args.strftime)
        except Exception as e:
            print('Unable to parse date time string: {}'.format(args.datetime))
            print(e)
            return

        if hasattr(args, 'epi_code'):
            assert type(args.epi_code) is int, \
                    """

You need to provide an EPI code via --epi-code, it must be integer,
within [-3,3], and not null.
"""
            df = pd.DataFrame({
                'cohort'   : args.cohort,
                'id'       : args.id[0],
                'date'     : date,
                'paradigm' : args.paradigm,
                'epi'      : args.epi_code,
                'valid'    : True,
                }, index=args.id)
        else:
            df = pd.DataFrame({
                'cohort'   : args.cohort,
                'id'       : args.id[0],
                'date'     : date,
                'paradigm' : args.paradigm,
                'valid'    : True,
                }, index=args.id)

        df.set_index(
                keys=['cohort', 'id', 'paradigm', 'date'],
                drop=False,
                inplace=True,
                verify_integrity=True)
        df.sort_index(inplace=True)
        df.cohort = df.cohort.astype('category')
        df.paradigm = df.paradigm.astype('category')

        return df

    elif fall_back is not None:
        try:
            instance = load(fall_back)
        except:
            try:
                instance = load(fall_back[0])
            except:
                instance = load(fall_back[1])

        if hasattr(args, 'epi_code'):
            assert type(args.epi_code) is int, \
                    """

You need to provide an EPI code via --epi-code, and it must be integer,
within [-3,3], and not null.
"""
            df = pd.DataFrame({
                'cohort'   : instance.name.cohort,
                'id'       : instance.name.j,
                'date'     : instance.name.datetime,
                'paradigm' : instance.name.paradigm,
                'epi'      : args.epi_code,
                'valid'    : True,
                }, index=[instance.name.j])
        else:
            df = pd.DataFrame({
                'cohort'   : instance.name.cohort,
                'id'       : instance.name.j,
                'date'     : instance.name.datetime,
                'paradigm' : instance.name.paradigm,
                'valid'    : True,
                }, index=[instance.name.j])

        df.set_index(
                keys=['cohort', 'id', 'paradigm', 'date'],
                drop=False,
                inplace=True,
                verify_integrity=True)
        df.sort_index(inplace=True)
        df.cohort = df.cohort.astype('category')
        df.paradigm = df.paradigm.astype('category')

        return df
