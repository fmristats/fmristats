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

from .lock import Lock

from .name import Identifier

from .affines import isclose

from .irritation import Block

from .session import Session

from .reference import ReferenceMaps

from .pmap import PopulationMap

from .smodel import SignalModel, Result

from pandas import DataFrame

import numpy as np

import pickle

def load(file):
    """
    Load instances from disk

    This will load a call instances from disk.

    Parameters
    ----------
    file : str
        File name.
    """
    with open(file, 'rb') as input:
        self = pickle.load(input)

    return self

def load_block_irritation(file, name, df, index, verbose=True):
    """
    Load a irritation file

    This will performe some tests whether the correct file has been
    loaded from disk.  This function is frequently used in the command
    line interface for fmristats.

    Parameters
    ----------
    file : str
    name : Identifier
    df : DataFrame
    verbose : bool

    Returns
    -------
    Irritation
    """
    assert type(name) is Identifier, 'name must be Identifier'
    assert type(df) is DataFrame, 'df must be DataFrame'

    try:
        irritation = load(file)
        if verbose:
            print('{}: Read: {}'.format(name.name(), file))
    except Exception as e:
        print('{}: Unable to read: {}'.format(name.name(), file))
        #print('{}: Exception: {}'.format(name.name(), e))
        df.ix[index,'valid'] = False
        return

    if type(irritation) is Block and name.is_equal(irritation.name):
        return irritation
    elif type(irritation) is Lock:
        df.ix[index,'valid'] = False
        print('{}: Irritation is currently locked'.format(name.name(), file))
        return irritation
    else:
        df.ix[index,'valid'] = False
        print("""{}: Irritation does not match:
        Expected: {}
        Found:    {}""".format(name.name(), name.name(True), irritation.name.name(True)))
        return irritation

def load_session(file, name, df, index, verbose=True):
    """
    Load a session file

    This will performe some tests whether the correct file has been
    loaded from disk.  This function is frequently used in the command
    line interface for fmristats.

    Parameters
    ----------
    file : str
    name : Identifier
    df : DataFrame
    verbose : bool

    Returns
    -------
    Session
    """
    assert type(name) is Identifier, 'name must be Identifier'
    assert type(df) is DataFrame, 'df must be DataFrame'

    try:
        session = load(file)
        if verbose:
            print('{}: Read: {}'.format(name.name(), file))
    except Exception as e:
        print('{}: Unable to read: {}'.format(name.name(), file))
        #print('{}: Exception: {}'.format(name.name(), e))
        df.ix[index,'valid'] = False
        return

    if type(session) is Session and name.is_equal(session.name):
        return session
    elif type(session) is Lock:
        df.ix[index,'valid'] = False
        print('{}: Session is currently locked'.format(name.name(), file))
        return session
    else:
        df.ix[index,'valid'] = False
        print("""{}: Session does not match:
        Expected: {}
        Found:    {}""".format(name.name(), name.name(True), session.name.name(True)))
        return session

def load_refmaps(file, name, df, index, verbose=True):
    """
    Load a file

    This will performe some tests whether the correct file has been
    loaded from disk.  This function is frequently used in the command
    line interface for fmristats.

    Parameters
    ----------
    file : str
    name : Identifier
    df : DataFrame
    index : int
    verbose : bool

    Returns
    -------
    ReferenceMaps
    """
    assert type(name) is Identifier, 'name must be Identifier'
    assert type(df) is DataFrame, 'df must be DataFrame'

    try:
        reference_maps = load(file)
        if verbose:
            print('{}: Read: {}'.format(name.name(), file))
    except Exception as e:
        print('{}: Unable to read: {}'.format(name.name(), file))
        #print('{}: Exception: {}'.format(name.name(), e))
        df.ix[index,'valid'] = False
        return

    if type(reference_maps) is ReferenceMaps and name.is_equal(reference_maps.name):
        return reference_maps
    elif type(reference_maps) is Lock:
        df.ix[index,'valid'] = False
        print('{}: ReferenceMaps is currently locked'.format(name.name(), file))
        return reference_maps
    else:
        df.ix[index,'valid'] = False
        print("""{}: ReferenceMaps does not match:
        Expected: {}
        Found:    {}""".format(name.name(), name.name(True), reference_maps.name.name(True)))
        return reference_maps

def load_result(file, name, df, index, vb, verbose=True):
    """
    Load a file

    This will performe some tests whether the correct file has been
    loaded from disk.  This function is frequently used in the command
    line interface for fmristats.

    Parameters
    ----------
    file : str
    name : Identifier
    df : DataFrame
    index : int
    verbose : bool

    Returns
    -------
    Result
    """
    assert type(name) is Identifier, 'name must be Identifier'
    assert type(df) is DataFrame, 'df must be DataFrame'

    try:
        result = load(file)
        if verbose:
            print('{}: Read: {}'.format(name.name(), file))
    except Exception as e:
        print('{}: Unable to read: {}'.format(name.name(), file))
        #print('{}: Exception: {}'.format(name.name(), e))
        df.ix[index,'valid'] = False
        return

    if type(result) is Result and name.is_equal(result.population_map.name):
        if vb is None:
            if np.isnan(result.statistics).all():
                df.ix[index,'valid'] = False
                print('{}: Result is completely empty'.format(name.name()))
                return result
            else:
                return result
        else:
            if vb == result.population_map.diffeomorphism.vb:
                if np.isnan(result.statistics).all():
                    df.ix[index,'valid'] = False
                    print('{}: Result is completely empty'.format(name.name()))
                    return result
                else:
                    return result
            else:
                df.ix[index,'valid'] = False
                print("""{}: Result does not match:
                Expected: {}
                Found:    {}""".format(name.name(), vb, population_map.diffeomorphism.vb))
                return result
    elif type(result) is Lock:
        df.ix[index,'valid'] = False
        print('{}: Result is currently locked'.format(name.name()))
        return result
    else:
        df.ix[index,'valid'] = False
        print("""{}: Result does not match:
        Expected: {}
        Found:    {}""".format(name.name(), name.name(True), result.name.name(True)))
        return result

def load_population_map(file, name, df, index, vb, verbose=True):
    """
    Load a file

    This will perform some tests whether the correct file has been
    loaded from disk.  This function is frequently used in the command
    line interface for fmristats.

    Parameters
    ----------
    file : str
    name : Identifier
    df : DataFrame
    index : int
    verbose : bool

    Returns
    -------
    PopulationMap
    """
    assert type(name) is Identifier, 'name must be Identifier'
    assert type(df) is DataFrame, 'df must be DataFrame'

    try:
        population_map = load(file)
        if verbose:
            print('{}: Read: {}'.format(name.name(), file))
    except Exception as e:
        print('{}: Unable to read: {}'.format(name.name(), file))
        #print('{}: Exception: {}'.format(name.name(), e))
        df.ix[index,'valid'] = False
        return

    if type(population_map) is PopulationMap and name.is_equal(population_map.diffeomorphism.nb):
        if vb is None:
            return population_map
        else:
            if vb == population_map.diffeomorphism.vb:
                return population_map
            else:
                df.ix[index,'valid'] = False
                print("""{}: PopulationMap does not match:
                Expected: {}
                Found:    {}""".format(name.name(), vb, population_map.diffeomorphism.vb))
                return population_map
    elif type(population_map) is Lock:
        df.ix[index,'valid'] = False
        print('{}: PopulationMap is currently locked'.format(name.name(), file))
        return population_map
    else:
        df.ix[index,'valid'] = False
        print("""{}: PopulationMap does not match:
        Expected: {}
        Found:    {}""".format(name.name(), name.name(True), population_map.name.name(True)))
        return population_map
