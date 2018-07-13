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

"""

Creates a scalar unit field for ATI reference by smoothing the input
image to given curvature.

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

    parser.add_argument('scale',
            type=float,
            help="""scale (standard deviation) of the kernel""")

    parser.add_argument('input',
            help="""input image file""")

    parser.add_argument('output',
            help="""output image file""")

    parser.add_argument('--ignore-mask',
            action='store_true',
            help="""ignore image mask""")

    parser.add_argument('--name',
            help="""if not given, the name will be set equal to the
            input image.""")

    parser.add_argument('--factor',
            type=float,
            default=3,
            help="""factor""")

    parser.add_argument('-v', '--verbose',
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

from ... import load

from ...diffeomorphisms import Image

from ...ati import fit_field, extract_field

import pandas as pd

import numpy as np

import time

########################################################################

def call(args):
    img = load(args.input)
    assert type(img) is Image

    if args.verbose > 1:
        print(img)

    if args.name:
        name = args.name
    else:
        name = img.name

    if not args.ignore_mask:
        img.mask()

    coordinates = img.coordinates()
    coordinates.shape

    if args.ignore_mask:
        mask = None
    else:
        mask          = img.get_mask()

    endog         = img.data.ravel()
    exog          = np.ones_like(endog).reshape((-1,1))
    agc           = coordinates.reshape((-1,3))

    old_settings = np.seterr(divide='raise', invalid='raise')
    time0 = time.time()

    statistics, parameter_dict, value_dict = fit_field(
            coordinates   = coordinates,
            mask          = mask,
            endog         = endog,
            exog          = exog,
            agc           = agc,
            scale         = args.scale,
            radius        = args.factor*args.scale,
            verbose       = args.verbose,
            name          = img.name)

    time1 = time.time()
    np.seterr(**old_settings)

    field = extract_field(
            field=statistics,
            param='intercept',
            value='point',
            parameter_dict={'intercept':0},
            value_dict=value_dict)

    result = Image(reference=img.reference, data=field, name=name)
    result.save(args.output)

    if args.verbose:
        time_spend = time1 - time0
        print('{}: Time needed for the fit: {:.2f} min'.format(
            result.name, time_spend / 60))
        print('{}: Time needed for the fit: {:.2f} h'  .format(
            result.name, time_spend / 60**2))
