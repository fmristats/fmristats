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

from .name import Identifier

import os

import datetime

import pickle

class Lock:
    """
    An instance of this class may by used to temporary look
    the file name while processing or fitting data.

    Parameters
    ----------
    name : Identifier
    who : str
        name of the program that currently locks this file.
    fname :
        file name that is being locked.
    timestamp : datetime
        time of lock.
    """
    def __init__(self, name, who, fname, timestamp=None):
        assert type(name) is Identifier, 'name must be an Identifier'

        if timestamp is None:
            timestamp = datetime.datetime.now()

        self.name = name
        self.fname = fname
        self.who = who
        self.timestamp = timestamp

    def describe(self, strftime='%Y-%m-%d-%H%M'):
        description = """
        Cohort:     {}
        Subject:    {:d}
        Paradigm:   {}
        File name:  {}
        Created on: {}
        Looked by:  {}
        """
        return description.format(
                self.name.cohort,
                self.name.j,
                self.name.paradigm,
                self.fname,
                self.timestamp.strftime(strftime),
                self.who,
                )

    def unlock(self):
        os.remove(self.fname)

    def conditional_unlock(self, df, index, verbose, force=False):
        if (not df.ix[index,'valid']) or force:
            if verbose:
                print('{}: Unlock: {}'.format(self.name.name(), self.fname))
            self.unlock()
            df.ix[index,'locked'] = False
            return True
        else:
            return False

    def save(self, file=None, **kwargs):
        """
        Save instance to disk

        This will save the instance to disk.

        Parameters
        ----------
        file : str
            File name.
        """
        if file is None:
            file = self.fname
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)
