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

Defines classes for affine transformations, rigid body transformations,
and orthonormal transformations, and functions for manipulating and
calculating inverses of affine transformations.

"""

from .euler import rotation_matrix_to_euler_angles

# TODO: make importing nibabel obsolete, here!
import nibabel as ni

import pickle

import numpy as np

from numpy.linalg import inv, det, norm, svd

########################################################################
#
# Basic functions
#
########################################################################

def cartesian2homogeneous(mat, vec):
    return ni.affines.from_matvec(mat, vec)

def from_cartesian(mat, vec):
    return Affine(cartesian2homogeneous(mat, vec))

def is_rotation_matrix(a, **kwargs):
    """
    Test if a 3x3-matrix is a rotation matrix

    Parameters
    ----------
    a : ndarray, shape (n,n)
        an darray
    **kwargs :
        parameters passed on to `np.close`.

    Returns
    -------
    bool :
        True or False.

    Notes
    -----
    The function assumes that `a` is square.
    """
    return np.isclose(
            np.dot(np.transpose(a), a),
            np.identity(3), **kwargs).all()

def isclose(x, y, **kwargs):
    """
    Test if two affines are close to each other

    Parameters
    ----------
    x : Affine
        an affine transformation.
    y : Affine
        an affine transformation.
    **kwargs :
        parameters passed on to `np.close`.

    Returns
    -------
    bool :
        True or False.
    """
    assert type(x) is Affine and type(y) is Affine, 'x and y must be of type Affine'
    return np.isclose(x.affine, y.affine, **kwargs).all()

def isinverse(x, y, **kwargs):
    """
    Test if two affines are inverses of each other

    Parameters
    ----------
    x : Affine
        an affine transformation.
    y : Affine
        an affine transformation.
    **kwargs :
        parameters passed on to `np.close`.

    Returns
    -------
    bool :
        True or False.
    """
    assert type(x) is Affine and type(y) is Affine, 'x and y must be of type Affine'
    return np.isclose(x.dot(y).affine, np.eye(4), **kwargs).all()


########################################################################
#
# Class for working with a single affine transformation
#
########################################################################

class Affine:
    """
    An affine transformation

    Affine transformations are stored in homologous coordinates.

    Parameters
    ----------
    affine : ndarray, shape (4,4)
        An affine transformation in homologous coordinates.
    """
    def __init__(self, affine):
        assert affine.shape == (4,4), 'wrong shape'
        assert (affine[3] == np.array([0,0,0,1])).all(), '3rd row must be [0,0,0,1]'
        self.is_rigid = is_rotation_matrix(affine[:3,:3])
        self.affine = affine

    def apply(self, x):
        """
        Apply affine transformation

        Will apply the affine transformation to the matrix or vector x.

        Parameters
        ----------
        x : ndarray, shape (3,) or (...,3)
            Vector or matrix of (row) vectors.

        Return
        ------
        ndarray:
            The result of the affine applied to the matrix (or vector)
            x, in particular, it will be of the same shape.
        """
        tmp = np.append(x,np.ones(x.shape[:-1] + (1,)), axis=-1)
        return np.einsum('...ij,...j', self.affine, tmp)[...,:3]

    def apply_to_index(self, index):
        """
        Apply affine transformation to index

        Parameters
        ----------
        index : 3-tuple
            3-tuple of int

        Returns
        -------
        ndarray:
            The result of the affine applied to the index.
        """
        return self.apply(np.array(index))

    def apply_to_indices(self, indices):
        """
        Apply affine transformation to indices

        Parameters
        ----------
        index : 3-tuple
            3-tuple of Slice or Ellipse objects.

        Returns
        -------
        ndarray:
            The result of the affine applied to the indices.
        """
        return self.apply(np.moveaxis(np.mgrid[indices], 0, -1))

    def inv(self):
        """
        Calculates the inverse of an affine transformation

        If the affine transformation is rigid, i.e., of the form A⋅x + b
        with A being a rotation matrix, then a more numerical stable
        algorithm for calculating the inverse is used than otherwise.

        Return
        ------
        The inverse of the affine transformation.
        """
        if self.is_rigid:
            u = self.affine[:3,:3]
            p = self.affine[:3, 3]

            inverse = np.empty_like(self.affine)
            inverse[:3,:3] =  u.T
            inverse[:3, 3] = -u.T.dot(p)
            inverse[3] = (0,0,0,1)

            return Affine(inverse)
        else:
            return Affine(inv(self.affine))

    def dot(self, x):
        """
        Matrix multiplication with x

        Parameters
        ----------
        x : Affine or Affines
            An affine transformation.

        Returns
        -------
        Affines or Affines:
            Output has the same type as x.
        """
        assert (type(x) is Affine) or (type(x) is Affines), 'x must be Affine or Affine'

        if type(x) is Affine:
            return Affine(self.affine.dot(x.affine))

        if type(x) is Affines:
            results = np.empty_like(x.affines)
            for t in range(x.n):
                results[t] = self.affine.dot(x.affines[t])
            return Affines(results)

    def euler(self):
        """
        Calculate Euler angles
        """
        assert self.is_rigid, 'does only work if affine is rigid'
        return rotation_matrix_to_euler_angles(self.affine[:3,:3])

    def resolution(self):
        """
        Resolution

        The resolution is the length of the base vectors (1,0,0),
        (0,1,0), and (0,0,1) after being mapped to by this affine
        transformations. The length of ``A⋅b``, for ``A`` being this
        affine transformation and ``b`` being one of the three basis
        vectors is the norm of ``A⋅b - A⋅0``.
        """
        return norm(self.affine, axis=0)[:3]

    def volume(self):
        """
        Volume of the parallelepiped spanned by the image of the
        canonical basis
        """
        return abs(det(self.affine))

    def diagonal(self):
        return norm(self.resolution())

    def aspect(self, x, y):
        """
        Parameters
        ----------
        x : int, either 0,1,2
            x
        y : int, either 0,1,2
            y

        Returns
        -------
        float :
            The aspect ratio of the respected resolution.
        """
        resolution = self.resolution()
        return resolution[x] / resolution[y]

    def describe(self, is_rigid=False):
        if is_rigid:
            description = """
        Resolution (left to right):         {:>5.2f} mm
        Resolution (posterior to anterior): {:>5.2f} mm
        Resolution (inferior to superior):  {:>5.2f} mm
        Diagonal of one voxel:              {:>5.2f} mm
        Volume of one voxel:                {:>5.2f} mm^3
        Aspect 0 on 1:                      {:>5.2f}
        Aspect 0 on 2:                      {:>5.2f}
        Aspect 1 on 2:                      {:>5.2f}
        Rigid transformation: {}"""
        else:
            description = """
        Resolution (left to right):         {:>5.2f} mm
        Resolution (posterior to anterior): {:>5.2f} mm
        Resolution (inferior to superior):  {:>5.2f} mm
        Diagonal of one voxel:              {:>5.2f} mm
        Volume of one voxel:                {:>5.2f} mm^3
        Aspect 0 on 1:                      {:>5.2f}
        Aspect 0 on 2:                      {:>5.2f}
        Aspect 1 on 2:                      {:>5.2f}"""

        resolution = self.resolution()
        return description.format(
                resolution[0],
                resolution[1],
                resolution[2],
                self.diagonal(),
                self.volume(),
                self.aspect(0,1),
                self.aspect(0,2),
                self.aspect(1,2),
                self.is_rigid
            )

    def save(self, file, **kwargs):
        """
        Save instance to disk

        Parameters
        ----------
        file : str
            A file name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)

########################################################################
#
# Class for working with an array of affine transformations
#
########################################################################

class Affines:
    """
    Affine transformations

    Affine transformations are stored in homologous coordinates.
    """
    def __init__(self, affines):
        """
        Parameters
        ----------
        affines : ndarray, shape (n,4,4)
            Affine transformations in homologous coordinates.
        """
        # TODO: allow higher dimensional indices: affines.shape = (m,n,4,4)

        if type(affines) is Affines:
            affines = affines.affines

        n, x, y = affines.shape
        assert (x,y) == (4,4), 'affines must be given in homologous coordinates'
        assert np.isclose(affines[:,3], np.array((0,0,0,1))).all(), 'last row must be [0,0,0,1]'
        rigid = np.repeat(False, n)
        for t in range(n):
            rigid[t] = is_rotation_matrix(affines[t,:3,:3])
        self.rigid = rigid.all()
        self.affines = affines
        self.n = n

    def index(self, j):
        return Affine(self.affines[j])

    def mean(self):
        """
        Calculates the mean affine transformation

        Returns
        -------
        Affine :
            The mean affine transformation.

        Notes
        -----
        Some basic non-Euclidean statistics is used here.
        """
        # https://en.wikipedia.org/wiki/Spherical_coordinate_system
        # https://en.wikipedia.org/wiki/Mean_of_circ
        m = self.affines.mean(axis=0)
        m[3] = (0,0,0,1)
        return Affine(m)

    def mean_rigid(self):
        """
        Calculates the mean rigid transformation of rigid
        transformations

        Returns
        -------
        Affine :
            The mean rigid affine transformation.

        Notes
        -----
        Some basic non-Euclidean statistics is used here.
        """
        assert self.rigid, 'affines must be rigid'
        m = self.mean().affine
        u, s, v = svd(m[:3,:3])
        m[:3,:3] = u.dot(v)
        return Affine(m)

    def mean_within_windows(self, r=0, skip=None):
        """
        Calculates the mean rigid transformation of the rigid
        transformation within a window of [-r,r] around each
        transformation.
        """
        assert self.rigid, 'affines must be rigid'
        assert r >= 0, 'r must be non-negative'

        mean_affines = self.affines.copy()

        if r > 0:
            if skip is None:
                for t in range(self.n):
                    xs = self.affines[max(0,t-r):t+r+1]
                    if xs.shape[0] > 0:
                        mean_affines[t] = Affines(xs).mean_rigid().affine
            else:
                for t in range(self.n):
                    xs = self.affines[max(0,t-r):t+r+1]
                    wh = np.where(~skip[max(0,t-r):t+r+1])
                    xs = xs[wh]
                    if xs.shape[0] > 0:
                        mean_affines[t] = Affines(xs).mean_rigid().affine

        mean_affines[:,3] = (0,0,0,1)
        return Affines(mean_affines)

    def euler(self):
        """
        Calculate Euler angles for all affine transformations

        Notes
        -----
        This only makes sense, if the affine transformations are rigid.
        """
        assert self.rigid, 'affines must be rigid'
        euler = np.zeros((self.n,3))
        for t in range(self.n):
            euler[t] = rotation_matrix_to_euler_angles(self.affines[t,:3,:3])
        return euler

    def apply(self, x):
        """
        Apply affine transformations to a vector

        Will apply the affine transformations to the vector x.

        Parameters
        ----------
        x : ndarray, shape (3,)
            Vector to which the affine transformations are applied

        Returns
        -------
        The result of the affine transformations applied to x.
        """
        assert x.shape == (3,), 'only makes sense if x is a 3-dim vector'
        tmp = np.append(x,np.ones(x.shape[:-1] + (1,)), axis=-1)
        return np.einsum('...ij,...j', self.affines, tmp)[...,:3]

    def apply_to_index(self, index):
        """
        Apply affine transformation to index

        Parameters
        ----------
        index : 3-tuple
            3-tuple of int

        Returns
        -------
        ndarray:
            The result of the affine applied to the index.
        """
        return self.apply(np.array(index))

    def inv(self):
        """
        Calculates the inverse of an affine transformation

        If the affine transformation is rigid, i.e., of the form A⋅x + b
        with A being a rotation matrix, then a more numerical stable
        algorithm for calculating the inverse is used than otherwise.

        Return
        ------
        The inverse of the affine transformation.
        """
        if self.rigid:
            u = self.affines[:,:3,:3]
            p = self.affines[:,:3, 3]
            uT = np.moveaxis(u, [1,2], [2,1])

            inverses = np.empty_like(self.affines)
            for t in range(self.n):
                inverses[t,:3,:3] =  uT[t]
                inverses[t,:3, 3] = -uT[t].dot(p[t])

            inverses[:,3] = (0,0,0,1)
            return Affines(inverses)
        else:
            inverses = np.empty_like(self.affines)
            for t in range(self.n):
                inverses[t] = inv(self.affines[t])
            return Affines(inverses)

    def dot(self, x):
        """
        Matrix multiplication (a.k.a. dot product) with the affines
        transformation x

        This will performe a matrix multiplication of each affine
        transformation with x.

        Parameters
        ----------
        x : Affine
            An affine transformation.

        Returns
        -------
        The result of the calculation is again of type Affines.
        """
        assert type(x) is Affine, 'x must be Affine'
        return Affines(self.affines.dot(x.affine))

    def save(self, file, **kwargs):
        """
        Save instance to disk

        Parameters
        ----------
        file : str
            A file name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)
