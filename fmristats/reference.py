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

The reference maps which map from a subject specific reference space to
the location and orientation of the subject's brain in the scanner
during the course of a FMRI session at a given time.

"""

from .tau import tau

from .affines import Affine, Affines

from .grubbs import grubbs

from .tracking import fit_by_pca

from .euler import rotation_matrix_to_euler_angles

import numpy as np

from numpy.linalg import inv, norm

import pickle

class ReferenceMaps:
    """
    The reference maps which map from a subject specific reference space
    to the location and orientation of the subject's brain in the
    scanner during the course of a FMRI session at a given time.

    Parameters
    ----------
    name : Identifier
    """

    def __init__(self, name):
        self.name = name

    def fit(self, session, use_raw=True):
        """
        Fit head movement

        This will fit rigid body transformations to the data that estimate
        position and bearing of the head.

        Parameters
        ----------
        session : Session
        use_raw : bool

        Notes
        -----
        The function implements the principle axis method for rigid body
        tracking.
        """
        self.shape = (session.numob, session.shape[session.ep])
        self.epi_code = session.epi_code
        self.ep = abs(self.epi_code)-1

        self.temporal_resolution = session.temporal_resolution
        self.slice_time = session.slice_time

        if use_raw:
            scan_references, w = fit_by_pca(
                    data=session.raw,
                    reference=session.reference)
        else:
            scan_references, w = fit_by_pca(
                    data=session.data,
                    reference=session.reference)

        self.semi_axis_norms = w

        self.set_scan_references(scan_references=scan_references)

    ####################################################################
    # Flights between space and time
    ####################################################################

    def set_scan_references(self, scan_references):
        """
        Sets the reference space of the fMRI experiment

        This will define the reference space of an fMRI experiment and
        calculates the corresponding observation grid of to this space.

        Parameters
        ----------
        scan_references : Affines or ndarray, shape (n,4,4), dtype: float
            The rigid body transformations.  Affine transformations that
            map coordinates, given with respect to the rigid body, to
            the coordinates of these points in (physical) scanner space
            at the time of the respected scan.

        Notes
        -----
        In particular, this function will set the scan references and
        their respected inverses.  In fact, the reference space is
        solely defined through these affine transformations.  A scan
        reference is a rigid body transformation that maps a point in
        the reference space `R` of a rigid body to the position of this
        point in the scanner `V` during the respected scan, in signs:

        :math:` R \to V`

        Note that the map goes from *reference* to *specific*, or from
        *the* scan to *some* scan, or from *fixed* to *moving*.

        The reference space is typically set by a function that fits
        head movements to observed intensities at a fixed (say, scanner
        specific) grid, e.g., `fit_head_movements`.
        """
        if type(scan_references) is Affines:
            self.scan_references = scan_references
        else:
            self.scan_references = Affines(scan_references)

    def reset_reference_space(self, x=None):
        """
        This will reset the reference space

        Parameters
        ----------
        x : None or Affine or ndarray, shape (4,4), dtype: float
            an affine transformation or None (default)

        Notes
        -----
        This will reset the reference space by moving origin and base
        vectors to the new position specified by the affine
        transformation x.  The affine transformation goes **from this**
        reference space **to the new** reference. If the affine
        transformation has the form Ax+b, then b is the new origin and A
        defines the new orienation.

        If x=None, then x will be set to the mean of the scan
        references.

        As scan references :math:`ρ_t` go from scanner space to
        reference space, the new scan reference are

        .. math::

            ρ_t ∘ x^{-1}

        A warning: this makes all population maps which have been
        defined for this session obsolete. You should run this
        function before setting any population maps!
        """
        if x is None:
            x = self.scan_references.mean_rigid()

        if type(x) is not Affine:
            x = Affine(x)

        assert type(x) is Affine, 'x must be Affine'

        self.reference_affine = x
        scan_references  = x.inv().dot(self.scan_references)
        self.set_scan_references(scan_references=scan_references)

    ####################################################################
    # Outlier detection
    ####################################################################

    def detect_outlying(self, sgnf):
        """
        Detect outlying scans

        Parameters
        ----------
        sgnf : float
            Significance level at which a test for the existence of a
            outlier is tested.

        Returns
        -------
        ndarray of bool
             True if the particular full scan cycle is considered an outlier.

        Notes
        -----
        Uses the eigenvalues of the pcm method to for outlier detection.
        (Currently only full scan cycles supported.)

        Remove scan cycles which have eigenvalues which differ
        significantly from all other volumes, as these likely are the
        result of severe head movement during this particular
        measurement period.

        Grubbs' test is used recursively for the outlier detection. The
        norm of all three semi axis length and each single semi axis
        length is tested for outlier separately.

        It should simply give you an idea on how severe head movements
        might have been, and allow you to remove the most obvious
        outliers.
        """
        # TODO:
        # 1. Implement 'less' in `grubbs_test` and change outlying
        #    cycle detection to less.
        # 2. Do this with respect to each slice not with respect to a
        #    full scan cycle...
        # 3. Use as the reference distribution for grubbs only scans
        #    within blocks of irritation

        # Look at the length of the semi axis norms
        w0 = self.semi_axis_norms[...,0]
        w1 = self.semi_axis_norms[...,1]
        w2 = self.semi_axis_norms[...,2]
        _, args0 = grubbs(w0, sgnf)
        _, args1 = grubbs(w1, sgnf)
        _, args2 = grubbs(w2, sgnf)

        # Look at the bary centres of the scan cycles
        x0 = self.scan_references.affines[:,0,3]
        x1 = self.scan_references.affines[:,1,3]
        x2 = self.scan_references.affines[:,2,3]
        _, args3 = grubbs(x0, sgnf)
        _, args4 = grubbs(x1, sgnf)
        _, args5 = grubbs(x2, sgnf)

        euler = self.scan_references.euler().T

        euler[ abs(euler) > tau / 5 ] = np.nan

        _, args6 = grubbs(euler[0], sgnf)
        _, args7 = grubbs(euler[1], sgnf)
        _, args8 = grubbs(euler[2], sgnf)

        values = np.vstack((
            w0, w1, w2,
            x0, x1, x2,
            euler[0], euler[1], euler[2]))

        outlying = np.vstack((
            args0, args1, args2,
            args3, args4, args5,
            args6, args7, args8))

        self.outlying = outlying
        self.values   = values

    def detect_outlying_scans(self, sgnf):
        self.detect_outlying(sgnf)

        self.outlying_cycles = self.outlying.any(axis=0)
        self.outlying_scans = np.repeat(
                self.outlying_cycles,
                self.shape[-1]).reshape(self.shape)

    ####################################################################
    # Descriptive statistics
    ####################################################################

    def descriptive_statistics(self, sgnf=.1):
        """
        Give descriptive statistics of the instance
        """
        self.detect_outlying_scans(sgnf=sgnf)
        n, m = self.shape
        o = self.outlying_cycles.sum()
        return n, m*n, o, m*o

    def describe(self):
        """
        Give a description of the instance
        """
        description = """
        Scans:        {:>4d}
        Valid:        {:>4d}
        Outlying:     {:>4d} ({:.2f}%)
        Scan cycles:  {:>4d}
        Valid:        {:>4d}
        Outlying:     {:>4d} ({:.2f}%)
        """
        n, m, on, om = self.descriptive_statistics()

        return description.format(
                m, m-om, om, 100*(om/m),
                n, n-on, on, 100*(on/n))

    ####################################################################
    # Save and load class instance to and from disk
    ####################################################################

    def save(self, file, **kwargs):
        """
        Save instance to disk

        This will save the instance to disk for later use.

        Parameters
        ----------
        file : str
            File name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)
