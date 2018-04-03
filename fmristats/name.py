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

import numpy as np

from datetime import datetime

class Identifier:
    """
    Identifier

    Parameters
    ----------
    cohort : str
    j : int
    datetime : datetime
    paradigm : str
    """

    template = '{}-{:04d}-{}'
    template_long = '{}-{:04d}-{}-{}'
    template_strftime = '%Y-%m-%d-%H%M'

    def __init__(self, cohort, j, datetime, paradigm):
        assert type(cohort) == str, 'cohort must be str'
        assert (type(j) == int) or (type(j) == np.int64), 'j must be integer'
        assert type(paradigm) == str, 'paradigm must be str'
        self.cohort = cohort
        self.j = j
        self.datetime = datetime
        self.paradigm = paradigm

    def is_equal(self, x):
        assert type(x) is Identifier

        if self.cohort   != x.cohort:
            return False
        if self.j        != x.j:
            return False
        if self.datetime != x.datetime:
            return False
        if self.paradigm != x.paradigm:
            return False

        return True

    def name(self, long=False):
        if long:
            return self.template_long.format(
                    self.cohort,
                    self.j,
                    datetime.strftime(self.datetime, self.template_strftime),
                    self.paradigm)
        else:
            return self.template.format(self.cohort, self.j, self.paradigm)

    def describe(self):
        description = """
        Cohort:   {}
        Subject:  {}
        Date:     {}
        Paradigm: {}"""
        return description.format(
                self.cohort,
                self.j,
                datetime.strftime(self.datetime, self.template_strftime),
                self.paradigm,
                )
