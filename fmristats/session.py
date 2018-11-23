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

Defines a class for FMRI sessions, which allows to define the
experimental setting and design, and to set up necessary parameters
needed for the statistical analysis of single subject data.

"""

from .name import Identifier

from .affines import Affine, isinverse

from .stimulus import Stimulus

from .filters import fit_foreground

from .diffeomorphisms import Image

import numpy as np

import pickle

#######################################################################
#######################################################################
#
# Describe the experimental design of an FMRI session
#
#######################################################################
#######################################################################

def create_slice_timing(shape, epi_code, temporal_resolution,
        interleaved=None):
    """
    Creates a time vector from parameters

    Parameters
    ----------
    shape : tuple, type: int
    epi_code : int or None
    temporal_resolution : float
        The temporal resolution.
    interleaved : bool

    Returns
    -------
    ndarray of float, shape (M,) or (M,N,O,P)
        The time vector that holds the slices time information. The
        vector holds the start point of the time at which a particular
        slice has been measured.  If ep is None, slice_time has
        shape (M,) with M equal to the number of full scan cycles.  If
        ep is int, then slice_time has
        the same shape as the data.

    """
    # TODO: implement interleaved
    assert not interleaved, 'sorry, interleaved not implemented yet'

    n = shape[0]
    m = shape[abs(epi_code)]
    slice_time = np.arange(
            start = 0,
            stop  = n*temporal_resolution,
            step  = temporal_resolution/m).reshape((n,m))
    return slice_time

def fmrisetup(session, stimulus):
    """

    If a session instance was created using the basic information in a
    Nifti1 file, this will complete the set-up of the session by
    providing the rest of the necessary FMRI session information:
    stimulus design (a.k.a paradigm).

    Parameters
    ----------
    session : Session
    stimulus : Stimulus
    """
    session.set_slice_timing()
    session.set_stimulus(stimulus)

#######################################################################

class Session:
    """
    A FMRI session

    Collects attributes which reasonably define a FMRI
    session, e.g., experimental design, the on- and offsets of the
    blocks or events of the paradigm, the observed FMRI signals.

    Parameters
    ----------
    name : Identifier
        The name or of the FMRI experiment
    data : array, shape (t,x,y,z)
        The data of observations, where t is the number of scan
        cycles
    epi_code : int
        Coded EPI sequence, either one of -3, -2, -1, 1, 2, 3
    spacial_resolution : ndarray, shape (3,), dtype: float
        Spacial resolution of the image
    temporal_resolution : float
        Temporal resolution of the image
    reference : Affine or ndarray, shape (4,4), dtype: float
        The scanner reference
    """

    def __init__(self, name, data, epi_code, spacial_resolution,
            temporal_resolution, reference):
        assert type(name) is Identifier, 'name must be of type Identifier'
        assert epi_code <=  3, 'epi_code must be ≤ 3'
        assert epi_code >= -3, 'epi_code must be ≥-3'
        assert epi_code !=  0, 'epi_code must be different from 0'

        self.name = name
        self.raw  = data      # will always contain the raw data
        self.data = self.raw  # will (one day) become a numpy masked array

        self.epi_code = epi_code
        self.ep = abs(epi_code)-1

        self.spacial_resolution = spacial_resolution
        self.temporal_resolution = temporal_resolution

        self.shape = self.raw.shape[1:]
        self.numob = self.raw.shape[0]

        if type(reference) is Affine:
            self.reference = reference
        else:
            self.reference = Affine(reference)

    ####################################################################
    # Functions that act on the data or extract information from it
    ####################################################################

    def set_slice_timing(self, slice_timing=None, interleaved=False):
        """
        Set up the slice time information

        Parameters
        ----------
        slice_timing : ndarray of float, shape (n,) or (n,x,y,z)
            The time vector indicating the point in time at which a
            given slice of the acquisition grid has been measured.  If
            slice_timing has shape (n,), all slices of a scan cycle are
            assumed to have been measured during the same time.
        """
        if slice_timing is None:
            slice_timing = create_slice_timing(
                shape               = self.data.shape,
                epi_code            = self.epi_code,
                temporal_resolution = self.temporal_resolution,
                interleaved         = False)

        shape_test = (self.numob, self.shape[self.ep])
        assert slice_timing.shape == shape_test, \
                'slice times have shape: {}, expected {}'.format(
                        slice_timing.shape, shape_test)

        self.slice_timing = slice_timing

    def set_stimulus(self, stimulus):
        """
        Set up the stimulus
        """
        assert issubclass(type(stimulus), Stimulus) or type(stimulus) is Stimulus, \
                'type must be Stimulus or a subclass if the same'
        self.stimulus = stimulus

    ####################################################################
    # Foreground detection
    ####################################################################

    def fit_foreground(self):
        self.data = self.raw.astype(float).copy()
        self.thresholds = fit_foreground(self.data, ep=self.ep)

    def set_foreground(self, foreground, is_mask=False):
        assert foreground.shape == data.shape, 'foreground and data not of equal shape'
        assert type(is_mask) is bool, 'is_mask must be bool'

        if is_mask:
            self.data = self.raw.astype(float).copy()
            self.data [ ~(foreground > 0) ] = np.nan
        else:
            self.data = foreground

    ###################################################################
    # Descriptive statistics of this session
    ###################################################################

    def descriptive_statistics(self):
        """
        Descriptive statistics for assessing the quality of the
        MR-signal on scan cycle level.
        """
        #TODO: give nicer description of what this function does
        ct, pt = self.thresholds
        dat = np.moveaxis(self.raw, (0, self.ep+1), (-1,-2))

        foreground = (dat >= pt) | (dat >= ct)
        background = (dat <  pt) | (dat >= ct)

        # Mean signal in the foreground
        muf = dat [ foreground ].mean()

        # Mean signal in the background
        mub = dat [ background ].mean()

        # Mean background to foreground signal
        mfb = mub / muf

        # Standard deviation in the foreground
        sdf = dat [ foreground ].std()

        # Standard deviation in the background
        sdb = dat [ background ].std()

        # foreground noise to mean signal foreground (foreground coefficient of variation)
        cvf = sdf / muf

        # background noise to mean signal foreground (background coefficient of variation)
        cvb = sdb / muf

        return muf, mub, mfb, sdf, sdb, cvf, cvb

    def descriptive_statistics_detailed(self):
        """
        Descriptive statistics for assessing the quality of the
        MR-signal on scan level.
        """
        #TODO: give nicer description of what this function does
        ct, pt = self.thresholds
        dat = np.moveaxis(self.raw, (0, self.ep+1), (-1,-2))

        background = (dat <  pt)
        foreground = ~background

        sumstats = np.zeros((pt.shape+(7,)))
        sumstats.shape

        fore = dat.astype(float)
        back = fore.copy()
        fore [ background ] = np.nan
        back [ foreground ] = np.nan

        muf = np.nanmean(fore, axis=(0,1)) # mean foreground signal
        mub = np.nanmean(back, axis=(0,1)) # mean background signal
        mfb = mub / muf                    # mean background to mean foreground signal
        sdf = np.nanstd(fore, axis=(0,1))  # standard deviation of foreground signal
        sdb = np.nanstd(back, axis=(0,1))  # standard deviation of foreground signal
        cvf = sdf / muf                    # foreground coefficient of variation
        cvb = sdb / muf                    # background coefficient of variation

        return muf, mub, mfb, sdf, sdb, cvf, cvb

    def describe(self):
        description = """
        EPI code: {:d}"""
        return description.format(
                self.epi_code,
                self.spacial_resolution,
                self.temporal_resolution,
                )

    #######################################################################
    # Save instance to disk
    #######################################################################

    def save(self, file, **kwargs):
        """
        Save model instance to disk

        This will save the current model instance to disk for later use.

        Parameters
        ----------
        file : str
            File name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)
