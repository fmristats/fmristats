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

Wrapper for FSL

"""

from .diffeomorphisms import Warp

from .nifti import nii2image, image2nii

import nibabel as ni

import numpy as np

from subprocess import run, PIPE

import os

from os.path import isdir

def bet(intercept, intercept_file, mask_file, cmd='fsl5.0-bet', variante='R',
        verbose=0):

    dfile = os.path.dirname(intercept_file)
    if dfile and not isdir(dfile):
       os.makedirs(dfile)

    ni.save(image2nii(intercept), intercept_file)

    command = [cmd]
    command.append(intercept_file)
    command.append(mask_file)
    command.append('-{}'.format(variante))

    if verbose:
        print()
        print('\n  '.join(command))

    command = ' '.join(command)

    try:
        bout = run([command], shell=True, stdout=PIPE, check=True)
        bout = bout.stdout.decode('utf-8').strip()
        if verbose > 1 and bout != '':
            print('Output of bet: {}'.format(bout))
    except Exception as e:
        print('Unable to run: {}'.format(command))
        print('Failed with: {}'.format(e))
        return

    try:
        mask = ni.load(mask_file)
    except Exception as e:
        print('Unable to read: {}'.format(mask_file))
        print('Failed with: {}'.format(e))
        return

    template = nii2image(mask)
    return template

def fnirt(warpcoef_file, nb_nii, vb_estimate_nii=None,
        vb_nii=None, vb_mask=None, nb_mask=None,
        config='T1_2_MNI152_2mm', cmd='fsl5.0-fnirt', verbose=True):
    """
    Run FNIRT to fit a diffeomorphism :math:`ψ` given data :math:`M` in
    the domain (vb) and data :math:`R` in the image (nb) of the
    diffeomorphism.

    Parameters
    ----------
    warpcoef_file : str
        File name in which to save the warp coefficients.
    nb_nii : str
        File name where to find :math:`ψ[m]`.
    vb_estimate_nii : None or str
        File name where to save :math:`ψ^{-1}[r]`.
    vb_nii : None or str
        File name where to find the data of :math:`m` (the template) in
        standard space that is used as reference for the domain (nb) of
        the diffeomorphism.
    vb_mask : None or str
        Name of file with mask in standard space.
    nb_mask : None or str
        Name of file with mask in subject reference space.
    config : str
        Name of the configuration file for FNIRT.
    cmd : str
        Name of the FSL FNIRT command line program.
    verbose : bool
        Control verbosity.
    """
    command = [cmd]

    command.append('--in={}'.format(nb_nii))

    if nb_mask:
        command.append('--inmask={}'.format(nb_mask))

    if vb_nii:
        command.append('--ref={}'.format(vb_nii))

    if vb_mask:
        command.append('--refmask={}'.format(vb_mask))

    if config:
        command.append('--config={}'.format(config))

    if vb_estimate_nii:
        command.append('--iout={}'.format(vb_estimate_nii))

    command.append('--cout={}'.format(warpcoef_file))

    if verbose:
        print()
        print('\n  '.join(command))

    command = ' '.join(command)

    try:
        sout = run([command], shell=True, stdout=PIPE, check=True)
        if verbose > 1:
            print(sout.stdout.decode('utf-8').strip())
        return True
    except Exception as e:
        print('Unable to run: {}'.format(command))
        print('Failed with: {}'.format(e))
        return False

def splines2warp(warpcoef_file, vb, vb_nii, nb_nii, name,
        coefficients_vb, coefficients_nb, new_diffeomorphism='fnirt',
        cmd='fsl5.0-std2imgcoord'):
    """
    This will run FSL img2stdcoord to turn the spline coefficient file
    produced by FNIRT into a diffeomorphism instance.

    Parameters
    ----------
    warpcoef_file : str
        File name of the spline coefficient file.
    vb : Image
        Standard space (VB).
    vb_nii : Image
        Path to standard space (VB) as Nifti1.
    nb_name : str or Identifier
        name or identifier for the image (nb).
    cmd : str
        Path to FSL std2imgcoord.
    coefficients_vb : str
        Coordinates of the template's image grid in standard space.
    coefficients_nb : str
        Coordinates of the template's image grid in subject space.

    Returns
    -------
    Warp : The diffeomorphism as a warp field.
    """
    vb_grid = vb.coordinates()
    vb_grid = vb_grid.reshape(-1,3)
    np.savetxt(coefficients_vb, X=vb_grid, delimiter=' ', fmt='%.2f')

    command = '{} -std {} -img {} -warp {} -mm {} > {}'.format(
            cmd, vb_nii, nb_nii, warpcoef_file,
            coefficients_vb, coefficients_nb)

    try:
        soutput = run([command], shell=True, stdout=PIPE, check=True)
    except Exception as e:
        print('Unable to run: {}, {}'.format(command, e))
        return

    try:
        coordinates = np.loadtxt(coefficients_nb)
    except Exception as e:
        print('Unable to read: {}, {}'.format(coefficients_nb, e))
        return

    coordinates = coordinates.reshape(vb.shape+(3,))

    return Warp(
            reference=vb.reference,
            warp=coordinates,
            vb=vb.name,
            nb=name,
            name=new_diffeomorphism,
            metadata={
                'vb_nii': vb_nii,
                'nb_nii': nb_nii,
                'splines' : warpcoef_file,
                }
            )
