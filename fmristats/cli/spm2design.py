# Copyright 2019 Thomas W. D. Möbius
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

Takes the SPM.mat file from SPM as input and outputs the design as a
(pandas and pickled) data frame.

"""

from ..epilog import epilog

import argparse

def define_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=epilog)

    parser.add_argument('mat', help="""Path to the SPM.mat file""")

    parser.add_argument('design', help="""Name of the output file""")

    parser.add_argument('-v', '--verbose',
            action='count',
            default=0,
            help="""Increase output verbosity""")

    return parser

def cmd():
    parser = define_parser()
    args = parser.parse_args()
    call(args, args.verbose)

cmd.__doc__ = __doc__

########################################################################
# Load libraries
########################################################################

import pandas as pd
from pandas import DataFrame

import numpy as np

from scipy.io import loadmat

########################################################################

def call(args, verbose):

    mat = loadmat(args.mat)

    design = mat['SPM']['xX'][0,0]['X'][0,0]
    design = DataFrame(design)
    design.to_pickle(args.design)

    if verbose:
        print("""The design matrix has shape: {}\n{}""".format(
            design.shape, design.round(3).head(12)))
