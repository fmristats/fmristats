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

Code for tracking a rigid body which has been sampled at a fixed lattice
in $ℝ^3$ by the principle axis method.

"""

import numpy as np

from fmristats.affines import Affine, Affines, cartesian2homogeneous

def fit_by_pcm(data, reference):
    """
    Tracking a rigid body

    This will track a rigid body, which has been sampled at a fixed
    lattice in :math:`ℝ^3`, by the so-called principle axis method.

    Parameters
    ----------
    data : ndarray, shape (n,x,y,z)
        Matrix of observations
    reference : Affine
        Affine transformation that maps an index to (physical) position

    Returns
    -------
    scan_inverse_references : Affines
        Rigid body transformations that map from the given coordinate
        system (specified by the affine `reference`) to the coordinate
        system of the body, in which zero corresponds to the bary centre
        of the body and the base axis correspond to its principle axis.

        Note that this is the other way around: In the terminology of
        the MB estimator, this corresponds to :math:`ρ_t^{-1}`.
    w : ndarray, shape (n,3)
        Principle semi axis lengths of the rigid body.

    Notes
    -----
    The function implements a (or the) principle axis method for rigid
    body tracking.  The method fits a orthonormal base (i.e. a
    coordinate system) onto data, such that the base vectors correspond
    to the principle axis (v) of the body.  The returned scan_inverse_references
    are affine transformations that map from the specified coordinate
    system to the coordinates of these points in the body.
    """
    numob, i, j, k = data.shape
    data = data.copy()
    data [ np.isnan(data) ] = 0

    n,x,y,z = data.shape
    indices = ((slice(0,x), slice(0,y), slice(0,z)))
    lattice = reference.apply_to_indices(indices)
    lattice = np.moveaxis(lattice, -1 ,0)

    # centre of masses
    com = np.empty((3,numob))
    for t in range(numob):
        com[:,t] = (data[t] * lattice).sum(axis=(1,2,3)) / data[t].sum()

    # demeaned lattice
    dm = np.repeat(lattice[...,None], numob, axis=-1) - com[:,None,None,None,:]

    # unweighted variance-covariance for each time point
    cov = dm[:,None] * dm

    # variance-covariance
    var = np.empty((3,3,numob))
    for t in range(numob):
        var[:,:,t] = (data[t] * cov[...,t]).sum(axis=(2,3,4)) / data[t].sum()

    # eigen decomposition: eigenvalues in w, eigenvectors in v
    var = np.moveaxis(var, -1, 0)
    w, v = np.linalg.eig(var)

    # scan inverse references
    scan_inverse_references = np.empty((numob,4,4))
    for t in range(numob):
        scan_inverse_references[t] = cartesian2homogeneous(v[t], com[:,t])

    return Affines(scan_inverse_references), w
