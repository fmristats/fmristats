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

Subjects may tend to move their head in the scanner ever so slightly,
which means that the tilt and position of the head during the FMRI
session may change over time. For this reason, we need to track the
subject head in the scanner, i.e. we need to find and estimate its tilt
and position at a given time point. Mathematically, head movements are
affine transformations (rigid body transformations) that map from a
space with a subject-specific coordinate system to the scanner space (or
in other words from a coordinate system that is fixed with respect to
the subject to a coordinate system that is fixed with respect to the
scanner). The *subject reference space* will always be denoted by
:math:`R`, and the scanner space by :math:`S`. We may then identify head
movements with the maps:

.. math::

    ρ_t : R \\to S.

These maps will be called reference maps: :math:`ρ_t` is called the
*reference map of scan* :math:`t`. Reference maps :math:`ρ_t` are
functions that map from subject reference space to the location and
orientation of the subject brain in the scanner at a given time point.

The inverse of :math:`ρ_t` is called the *acquisition map of scan*
:math:`t`:

.. math::

    ρ_t^{-1} : S \\to R.

Hence, the acquisition map of scan :math:`t` maps from scanner space to
subject reference space.

"""

from .tau import tau

from .affines import Affine, Affines

from .grubbs import grubbs

from .tracking import fit_by_pcm

from .euler import rotation_matrix_to_euler_angles

import numpy as np

from numpy.linalg import inv, norm

import pickle

class ReferenceMaps:
    """
    Reference maps are functions that map from the subject reference
    space to the location and orientation of the subject brain in the
    scanner at a given time point.

    Parameters
    ----------
    name : Identifier
        The identifier of the subject
    """

    def __init__(self, name):
        self.name = name
        self.outlying_scans = None
        self.outlying_cycles = None
        self.outlying = None

    def fit(self, session, use_raw=True):
        """
        Fit head movement

        This will fit rigid body transformations to the data, i.e. it
        will estimate the position and bearing of the head in each scan
        cycle.

        Parameters
        ----------
        session : Session
            The FMRI session data of the subject.
        use_raw : bool
            Use the raw data or only the data in the foreground.

        Notes
        -----
        The function implements a principle axis method for rigid body
        tracking.
        """
        self.shape = (session.numob, session.shape[session.ep])
        self.epi_code = session.epi_code
        self.ep = abs(self.epi_code)-1

        self.temporal_resolution = session.temporal_resolution
        self.slice_timing = session.slice_timing

        if use_raw:
            acquisition_maps, w = fit_by_pcm(
                    data=session.raw,
                    reference=session.reference)
        else:
            acquisition_maps, w = fit_by_pcm(
                    data=session.data,
                    reference=session.reference)

        self.reference       = session.reference
        self.semi_axis_norms = w

        self.set_acquisition_maps(maps = acquisition_maps)

    ####################################################################
    # Flights between space and time
    ####################################################################

    def set_acquisition_maps(self, maps):
        """
        The subject reference space of the fMRI experiment

        This will define the subject reference space of an FMRI
        experiment.

        Parameters
        ----------
        maps : Affines or ndarray, shape (n,4,4), dtype: float
            Acquisition maps. Rigid body transformations that
            map from scanner space to subject reference space.

        Notes
        -----
        The subject reference space is uniquely defined by these
        affine transformations.
        """
        if type(maps) is Affines:
            self.acquisition_maps = maps
        else:
            self.acquisition_maps = Affines(maps)

    def reset_reference_space(self, x=None, cycle=None):
        """
        This will reset the coordinates system of the reference space

        Parameters
        ----------
        x : None or Affine or ndarray, shape (4,4), dtype: float
            An affine transformation or None (default)
        cycle : None or int
            Index of a scan cycle or None

        Notes
        -----
        This will reset the coordinates system of the subject reference
        space by moving origin and base vectors to the new position
        specified by the transformation x. The affine transformation
        :math:`x` goes **from this** reference space **to the new**
        reference space.

        .. math::

            x : R \\to R'.

        If the affine transformation has the form
        x(:math:`x`)= :math:`Ax+b`, then :math:`b` is the new origin and
        :math:`A` defines the new orientation.

        If x=None, then x will be set to the inverse of the mean of the
        acquisition maps. This has the consequence that the new
        reference space is identical to the average position of the
        subject head in the scanner. If x=None and cycle is
        integer :math:`t_0`, then :math:`x` is set to:

        .. math::

            x = ρ_{t_0}.

        This has the consequence that the new subject reference space
        equals the position of the subject head during scan cycle
        :math:`t_0`.

        References maps are the maps :math:`ρ_t` which map from subject
        reference space to scanner:

        .. math::

            ρ_t : R \\to S.

        Acquisition maps are the inverses of the maps :math:`ρ_t` and
        map from scanner space to subject reference space:

        .. math::

            ρ_t^{-1} : S \\to R.

        The new subject reference space :math:`R'` is equal to
        :math:`R':=x[R]`, and the new acquisition maps are therefore:

        .. math::

            ρ_t^{'-1} = x ∘ ρ^{-1}_t

        Then:

        .. math::

            ρ'_t : R' \\to S

        **Warning**: Resetting the reference space makes all population
        maps which have been defined for this subject reference space
        obsolete. You should thus perform this operation prior to
        fitting any population maps.
        """
        if x is None:
            if cycle is None:
                x = self.acquisition_maps.mean_rigid().inv()
            else:
                x = Affine(self.acquisition_maps.affines[cycle]).inv()

        if type(x) is not Affine:
            x = Affine(x)

        assert type(x) is Affine, 'x must be Affine'
        assert x.is_rigid, 'x must be rigid'

        self.x = x
        self.cycle = cycle

        new_acquisition_maps = x.dot(self.acquisition_maps)
        self.set_acquisition_maps(maps = new_acquisition_maps)

    def mean(self, r, skip_outlying=True):
        assert r >= 0, 'r must be non-negative'

        if r > 0:
            if skip_outlying \
                    and (self.outlying_cycles is not None) \
                    and  self.outlying_cycles.any():
                self.acquisition_maps = self.acquisition_maps. \
                        mean_within_windows(r, skip=self.outlying_cycles)
            else:
                self.acquisition_maps = self.acquisition_maps. \
                        mean_within_windows(r, skip=None)

    ####################################################################
    # Outlier detection
    ####################################################################

    def detect_outlying_scans(self, sgnf):
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
        Uses the eigenvalues of the pcm method for outlier detection.
        (Currently only full scan cycles supported.)

        The algorithm will mark scan cycles as False which have
        eigenvalues which differ significantly from the eigenvalues all
        other scan cycles, as this likely is the result of severe head
        movement during the measurement of this cycle.

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
        #    within blocks of stimulus

        # Look at the length of the semi axis norms
        w0 = self.semi_axis_norms[...,0]
        w1 = self.semi_axis_norms[...,1]
        w2 = self.semi_axis_norms[...,2]
        _, args0 = grubbs(w0, sgnf)
        _, args1 = grubbs(w1, sgnf)
        _, args2 = grubbs(w2, sgnf)

        # Look at the bary centres of the scan cycles
        x0 = self.acquisition_maps.affines[:,0,3]
        x1 = self.acquisition_maps.affines[:,1,3]
        x2 = self.acquisition_maps.affines[:,2,3]
        _, args3 = grubbs(x0, sgnf)
        _, args4 = grubbs(x1, sgnf)
        _, args5 = grubbs(x2, sgnf)

        euler = self.acquisition_maps.euler().T

        euler[ abs(euler) > tau / 5 ] = np.nan

        _, args6 = grubbs(euler[0], sgnf)
        _, args7 = grubbs(euler[1], sgnf)
        _, args8 = grubbs(euler[2], sgnf)

        outlying = np.vstack((
            args0, args1, args2,
            args3, args4, args5,
            args6, args7, args8))

        self.outlying = outlying

        self.outlying_cycles = self.outlying.any(axis=0)

        self.outlying_scans = np.repeat(
                self.outlying_cycles,
                self.shape[-1]).reshape(self.shape)

    ####################################################################
    # Descriptive statistics
    ####################################################################

    def descriptive_statistics(self):
        """
        Give descriptive statistics of the instance
        """
        # TODO: check if outlying_cycles exists!
        n, m = self.shape
        o = self.outlying_cycles.sum()
        return n, m*n, o, m*o

    def describe(self):
        """
        Give a description of the instance
        """
        # TODO: check if outlying_cycles exists!
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
