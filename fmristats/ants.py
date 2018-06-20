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

Wrapper for ANTS

"""

from .pmap import PopulationMap

from .diffeomorphisms import Warp

from .nifti import nii2image, image2nii

import nibabel as ni

import numpy as np

import os

from os.path import isfile, isdir, join

from nipype.interfaces.ants import RegistrationSynQuick, ApplyTransformsToPoints

from contextlib import contextmanager

import sys

# There seems to be a bug either in nipype.interfaces or ANTS. Expected
# by nipype is
#     RegistrationSynQuick
# provided by ANTS is
#     RegistrationSyNQuick
# whether this is due to the installation by NeuroDebian?

def fit_population_map(vb_image, nb_image, output_prefix, threads=4,
        verbose=True, vb=None, name='ants'):
    """
    Fits a diffeomorphism ψ from `vb` (the domain of ψ) to `nb` (the
    image of ψ) using the images `vb_image` and `nb_image` as references in
    `vb` and `nb` respectively.

    Parameters
    ----------
    vb_image : Image
        Domain / Vorbereich
    nb_image : Image
        Image / Nachbereich
    """

    dfile = os.path.dirname(output_prefix)
    if dfile and not isdir(dfile):
       os.makedirs(dfile)

    vb_file = output_prefix + 'vb.nii.gz'
    nb_file = output_prefix + 'nb.nii.gz'

    ni.save(image2nii(vb_image), vb_file)
    ni.save(image2nii(nb_image), nb_file)

    reg = RegistrationSynQuick()
    reg.inputs.fixed_image  = vb_file
    reg.inputs.moving_image = nb_file
    reg.inputs.num_threads  = threads
    reg.inputs.output_prefix = output_prefix

    if verbose:
        print()
        print(reg.cmdline)

    reg.run()

    vb_coord = output_prefix + 'vb-coordinates.csv'
    nb_coord = output_prefix + 'nb-coordinates.csv'

    vb_grid = vb_image.coordinates()
    vb_grid[...,:2] = -vb_grid[...,:2]
    vb_grid = vb_grid.reshape(-1,3)

    np.savetxt(vb_coord, X=vb_grid, delimiter=',', fmt='%.2f',
        header='x,y,z', comments='')

    transM  = output_prefix + '0GenericAffine.mat'
    transW  = output_prefix + '1Warp.nii.gz'
    transWI = output_prefix + '1InverseWarp.nii.gz'

    at = ApplyTransformsToPoints()
    at.inputs.dimension  = 3
    at.inputs.input_file = vb_coord
    at.inputs.transforms = [transW, transM]
    at.inputs.invert_transform_flags = [False, False]
    at.inputs.output_file = nb_coord

    if verbose:
        print()
        print(at.cmdline)

    at.run()

    coordinates = np.loadtxt(nb_coord, delimiter=',', skiprows=1)
    coordinates = coordinates.reshape(vb_image.shape+(3,))
    coordinates[...,:2] = -coordinates[...,:2]

    diffeomorphism = Warp(
            reference=vb_image.reference,
            warp=coordinates,
            vb=vb_image.name,
            nb=nb_image.name,
            name=name,
            metadata={
                'vb_file': vb_file,
                'nb_file': nb_file,
                'transW' : transW,
                'transM' : transM,
                'RegistrationSyNQuick' : reg.cmdline,
                'ApplyTransformsToPoints' : at.cmdline
                }
            )

    try:
        vb_estimate_file = output_prefix + 'Warped.nii.gz'
        vb_estimate = nii2image(ni.load(vb_estimate_file), name='vb_estimate')
    except Exception as e:
        print('{}: Unable to read: {}'.format(nb_image.name.name(), vb_estimate_file))
        print('{}: Exception: {}'.format(nb_image.name.name(), e))

    try:
        nb_estimate_file = output_prefix + 'InverseWarped.nii.gz'
        nb_estimate = nii2image(ni.load(nb_estimate_file), name='nb_estimate')
    except Exception as e:
        print('{}: Unable to read: {}'.format(nb_image.name.name(), nb_estimate_file))
        print('{}: Exception: {}'.format(nb_image.name.name(), e))

    return PopulationMap(diffeomorphism,
            vb=vb_image,
            nb=nb_image,
            vb_estimate=vb_estimate,
            nb_estimate=nb_estimate,
            name=vb,
            )

# Some remarks:
#
# From:
#     https://sourceforge.net/p/advants/discussion/840261/thread/2a1e9307/
#
# The input and output point coordinates are in "physical-LPS"
# coordinates, that is, they are the coordinates encoded in the nifti
# header (which are in RAS orientation), but with X and Y's sign flipped
# (so it becomes LPS orientation). To use the function, I first do -X
# and -Y for the points I want to use, feed them to the function, and
# again take -X and -Y of what comes out to recover the physical
# coordinates I expected.
#
# With one word: strange.
