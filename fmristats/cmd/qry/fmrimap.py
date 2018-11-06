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

Functions that return the image of a coordinate given in one space in
the coordinates of the image space of a given map.

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

    parser.add_argument('input',
        help="""a population space, a population map, any
        diffeomorphism, or an image file.""")

    parser.add_argument('first',  help='first coordinate',  type=float)
    parser.add_argument('second', help='second coordinate', type=float)
    parser.add_argument('third',  help='third coordinate',  type=float)

    parser.add_argument('-i', '--index',
        action='store_true',
        help="""by default coordinates are assumed to live in Euclidean
        space. If --index is set, coordinate are interpreted to live in
        index space instead.""")

    parser.add_argument('--inverse',
        action='store_true',
        help="""inverse""")

    parser.add_argument('-v', '--verbose',
            action='store_true',
            help='verbose output')

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

from ...load import load

from ...pmap import PopulationMap

from ...smodel import Result

from ...diffeomorphisms import Warp, Displacement

import sys

import numpy as np

########################################################################

def call(args):
    """
    Functions that return the image of a coordinate given in one space in
    the coordinates of the image space of a given map.

    Parameters
    ----------
    args : Arguments
    """
    try:
        x = load(args.input)
    except Exception as e:
        print('Cannot read: {}'.format(args.input))
        print('Failed with: {}'.format(e))
        sys.exit()

    if type(x) is PopulationMap:
        x = x.diffeomorphism

    if type(x) is Result:
        x = x.population_map.diffeomorphism

    if args.index:
        index = tuple(np.array((
            args.first,
            args.second,
            args.third), dtype=int))

        if args.verbose:
            print(
"""index (in RAS+): {}
has coordinates in the domain: {}
and coordinates in the image:  {}""".format(
                index,
                x.reference.apply_to_index(index),
                x.apply_to_index(index)))
        else:
            print(x.apply_to_index(index))

    else:
        coordinate = np.array((
            args.first,
            args.second,
            args.third), dtype=float)

        if args.inverse:
            if type(x) is Warp or type(x) is Displacement:
                coordinates = x.coordinates()
                distances   = ((coordinates - coordinate)**2).sum(axis=-1)
                domain_index = np.unravel_index(np.argmin(distances), distances.shape)
                domain_coordinate = x.reference.apply_to_index(domain_index)
                if args.verbose:
                    print(
    """index (in RAS+): {}
    has coordinates in the domain: {}
    and coordinates in the image:  {}""".format(
                        domain_index,
                        domain_coordinate,
                        coordinate))
                else:
                    print(domain_coordinate)
            else:
                print('fmrimap --inverse is currently only implemented for Warp and Displacement')

        else:
            if args.verbose:
                print(
"""index (in RAS+): {}
has coordinates in the domain: {}
and coordinates in the image:  {}""".format(
                    x.reference.inv().apply(coordinate),
                    coordinate,
                    x.apply(coordinate)))
            else:
                print(x.apply(coordinate))
