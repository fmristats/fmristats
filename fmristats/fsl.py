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

from .pmap import PopulationMap

from .diffeomorphisms import Warp

from .nifti import nii2image, image2nii

import nibabel as ni

import numpy as np

from subprocess import PIPE

from subprocess import run

import os

from os.path import isfile, isdir, join

def bet(intercept, intercept_file, mask_file, cmd, variante='R',
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

def fit_warpcoef(nb_file, warpcoef_file, preimage_file=None,
        vb_file=None, vb_mask=None, nb_mask=None,
        config='T1_2_MNI152_2mm', cmd='fsl5.0-fnirt', verbose=True):
    """
    Run FNIRT to fit a diffeomorphism :math:`ψ` given data :math:`M` in
    the domain (vb) and data :math:`R` in the image (nb) of the
    diffeomorphism.

    Parameters
    ----------
    nb_file : str
        File name where to find the data of :math:`R=ψ[M]` which shall
        be used as reference for the image space (nb) of the
        diffeomorphism.
    warpcoef_file : str
        File name in which to save the warp coefficients.
    target_file : None or str
        File name in which to save the target :math:`ψ^{-1}[R]`. If
        None, the target will not be saved. This is the default.
    config : str
        Name of the configuration file FNIRT shall use.
    vb_file : None or str
        File name where to find the data of :math:`M` (the template)
        which shall be used as reference for the domain space (nb) of
        the diffeomorphism.
    vb_mask : None or str
        Name of file with mask in reference space. May be provided if
        vb_file is given.
    cmd : str
        Name of the FSL's FNIRT command line program.
    verbose : bool
        Print concatinated command to stdout
    """
    command = [cmd]

    command.append('--in={}'.format(nb_file))

    if nb_mask:
        command.append('--inmask={}'.format(nb_mask))

    if vb_file:
        command.append('--ref={}'.format(vb_file))

    if vb_mask:
        command.append('--refmask={}'.format(vb_mask))

    if config:
        command.append('--config={}'.format(config))

    if preimage_file:
        command.append('--iout={}'.format(preimage_file))

    command.append('--cout={}'.format(warpcoef_file))

    if verbose:
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

def warpcoef2pmap(warpcoef_file, vb_file, vb_name, nb_file, nb_name,
        cpopulation_file, csubject_file, cmd='fsl5.0-std2imgcoord'):
    """
    Run FSL's img2stdcoord to turn the warp coefficient file produced by
    FNIRT into a fmristats' diffeomorphism intance.

    Parameters
    ----------
    vb_file : str
        file name where to find the data of :math:`M` (the template)
        which shall be used as reference for the domain space (nb) of
        the diffeomorphism, i.e. file name in which the template is
        stored.
    nb_file : str
        file name where to find the data of :math:`R=ψ[M]` which shall
        be used as reference for the image space (nb) of the
        diffeomorphism.
    vb_name :
        name or identifyer for the domain (vb).
    nb_name :
        name or identifyer for the image (nb).
    warpcoef_file : str
        File name of the warp coefficient file.
    cmd : str
        Name of FSL's std2imgcoord command line program.
    cpopulation_file : str
        Coordinates of the template's image grid in population space
    csubject_file : str
        Coordinates of the template's image grid in subject space

    Returns
    -------
    ndarray, shape template.shape
        The coordinates.
    """
    template = nii2image(ni.load(vb_file))
    template_grid = template.coordinates()
    template_grid = template_grid.reshape(-1,3)
    np.savetxt(cpopulation_file, X=template_grid, delimiter=' ', fmt='%i')

    command = '{} -std {} -img {} -warp {} -mm {} > {}'.format(
            cmd, vb_file, nb_file, warpcoef_file,
            cpopulation_file, csubject_file)

    try:
        soutput = run([command], shell=True, stdout=PIPE, check=True)
    except Exception as e:
        print('Unable to run: {}'.format(command))
        print('Failed with: {}'.format(e))
        return

    try:
        coordinates = np.loadtxt(csubject_file)
    except Exception as e:
        print('Unable to read: {}'.format(csubject_file))
        print('Failed with: {}'.format(e))
        return

    coordinates = coordinates.reshape(template.shape+(3,))

    diffeomorphism = Warp(
            reference=template.reference,
            warp=coordinates,
            vb=vb_name,
            nb=nb_name,
            metadata={
                'vb_file':vb_file,
                'nb_file':nb_file,
                'warpcoef':warpcoef_file,
                }
            )

    return PopulationMap(diffeomorphism, template=template)
