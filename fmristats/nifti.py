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

from . import Image, Session

import numpy as np

import nibabel as ni

#######################################################################
# Image
#######################################################################

def nii2image(nii, name=None):
    """
    Converts Nifti1Image to Image
    """
    data = nii.get_data()
    assert len(data.shape) == 3, 'not a 3D-image'
    return Image(reference=nii.affine, data=data, name=name)

def image2nii(image, nan_to_zero=False):
    """
    Converts Image to Nifti1Image
    """
    data = image.data.copy()
    if nan_to_zero:
        data [ ~image.get_mask() ] = 0
    return ni.Nifti1Image(data, image.reference.affine)

#######################################################################
# Session
#######################################################################

def nii2session(name, nii, epi_code):
    """
    Create a session instance from a Nifti1Image

    Parameters
    ----------
    name : str
    nii : nibabel.nifti1.Nifti1Image
    plain : int
    """
    zooms = nii.header.get_zooms()
    return Session(name=name,
            data  = np.rollaxis(nii.get_data(), -1),
            epi_code = epi_code,
            spacial_resolution = zooms[:-1],
            temporal_resolution = zooms[-1],
            reference = nii.affine)
