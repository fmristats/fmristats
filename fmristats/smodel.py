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

Defines a class for the FMRI signal model and its fits.

"""

from .affines import Affine

from .name import Identifier

from .session import Session

from .reference import ReferenceMaps

from .diffeomorphisms import Image

from .pmap import PopulationMap

from .stimulus import Stimulus

from .fit import fit_field, extract_field, model_at, data_at, fit_at

import time

import numpy as np

from numpy.linalg import inv

from numpy.linalg import norm

import scipy.stats.distributions as dist

from patsy import dmatrix

from pandas import DataFrame

import pickle

import math

class SignalModel:
    """
    The signal model

    Parameters
    ----------
    session : Session
        A session (with fitted foreground).
    reference_maps : ReferenceMaps
        The estimated reference maps.
    population_map : PopulationMap
        A population map
    """
    def __init__(self, session, reference_maps, population_map):
        assert type(session) is Session, 'session must be of type Session'
        assert type(reference_maps) is ReferenceMaps, 'reference_maps must be of type ReferenceMaps'
        assert type(population_map) is PopulationMap, 'population_map must be of type PopulationMap'

        self.session = session
        self.stimulus = session.stimulus
        self.reference_maps = reference_maps
        self.population_map = population_map

        self.shape = (session.numob, session.shape[session.ep])
        self.slice_timing = session.slice_timing

        self.name = self.session.name
        self.epi_code = self.session.epi_code
        self.ep = abs(self.epi_code) - 1 # ep in [0,1,2]
        assert self.ep in [0,1,2], 'ep is not in [0,1,2]'

        # The references defined here are the affine maps that map from
        # the *index space* of the Session to the subject reference
        # space defined by the acquisition maps in the ReferenceMaps.

        self.references = reference_maps.acquisition_maps.dot(session.reference)

        # TODO: when having added the option: --reference-maps None,
        # replace this with something reasonable

    #######################################################################
    # Set hyperparameters for the fit
    #######################################################################

    def set_hyperparameters(self, scale_type='max', factor=3,
            scale=None, mass=None):
        """
        Set up hyperparameters for field estimation

        All parameters which control details the effect field estimation
        and which are not per se sample design specific are called
        Hyperparameters.

        Parameters
        ----------
        scale_type : str
            One of min, max (default), diagonal
        scale : float
            Standard deviation of the Gaussian used to weight
            observations with respect to distance to point.
        mass : float (optional)

        Notes
        -----
        The only hyperparameter you need to set is `scale`, which is set
        for you if you set the scale_type.

        Factor takes precedence to mass, as I will assume if you are
        setting the mass, you know what you are doing.
        """
        if scale is None:
            if scale_type == 'diagonal':
                self.scale_type = scale_type
                self.scale = 0.5 * self.session.reference.diagonal()
            elif scale_type == 'min':
                self.scale_type = scale_type
                self.scale = 0.5 * min(self.session.reference.resolution())
            elif scale_type == 'max':
                self.scale_type = scale_type
                self.scale = 0.5 * max(self.session.reference.resolution())
            else:
                raise ValueError('If SCALE is not set, SCALE_TYPE is compulsory.')
        else:
            self.scale_type = 'user'
            self.scale = scale

        if factor is not None:
            assert factor > 0, 'factor must be strictly positve'
            self.factor = factor
            self.radius = factor * self.scale
            self.mass   = dist.norm.cdf(factor) - dist.norm.cdf(-factor)

        if mass is not None:
            assert mass > 0, 'mass must be strictly positive'
            assert mass < 1, 'mass must be strictly less than 1'
            self.mass = mass
            self.factor = dist.norm.ppf(1 - (1-mass)/2)
            self.radius = self.factor * self.scale

    #######################################################################
    # The mean index affine
    #######################################################################

    def mean_index_affine(self):
        """
        Affine transformation from reference space to the index space of
        the acquisition grid

        Returns
        -------
        Affine : The affine transformation that maps a coordinate in
        reference space to the mean index in the index space of the
        acquisition grid of the session.
        """
        return self.references.inv().mean()

    #######################################################################
    # All about the observations
    #######################################################################

    def coordinates(self):
        """
        The coordinate grid of observations

        The coordinates in subject reference space of all points at
        which MR signals have been acquired.

        Returns
        -------
        ndarray, shape (n,x,y,z,3), dtype: float
            [...,:3] = coordinates of observation
        """
        n,x,y,z = self.session.data.shape
        indices = ((slice(0,x), slice(0,y), slice(0,z)))

        coordinates = np.empty((n,x,y,z,3))
        for t in range(n):
            coordinates[t] = self.references.index(t).apply_to_indices(indices)

        return coordinates

    def set_stimulus_design(self, **kwargs):
        self.stimulus_design = self.stimulus.design(
                slice_timing=self.slice_timing, **kwargs)
        return self.stimulus_design

    def get_observations(self, include_background=False):
        """
        Observation matrix

        Returns the array of (un-cleaned) observations

        Returns
        -------
        ndarray, shape (n,x,y,z,9), dtype: float
            [...,:3] = coordinates of observation
            [..., 3] = MR signal response
            [..., 4] = time of observation
            [..., 5] = task during time of observation
            [..., 6] = block number during time of observation
            [..., 7] = scan cycle
            [..., 8] = slice number
        """
        assert hasattr(self, 'stimulus_design'), 'first set the stimulus design'
        n,x,y,z = self.session.data.shape

        observations = np.zeros((n,x,y,z,9))
        observations[...,:3] = self.coordinates()

        if include_background:
            observations[..., 3] = self.session.raw
        else:
            observations[..., 3] = self.session.data

        tmp = np.moveaxis(observations, self.ep+1, 1)
        tmp[...,4]  = self.slice_timing[...,None,None]
        tmp[...,5:7] = self.stimulus_design[...,None,None,:]

        x = np.mgrid[:self.shape[0], :self.shape[1]]
        x = np.moveaxis(x, 0, -1)
        tmp[...,7:9] = x[...,None,None,:]

        return observations

    def set_data(self, burn_in=4, demean=True, dropna=True,
            include_background=False, verbose=True):
        """
        Set the observation matrix

        Notes
        -----

        The observations will be an array of the following shape:

            ndarray, shape (n,x,y,z,9), dtype: float
                [...,:3] = coordinates of observation
                [..., 3] = MR signal response
                [..., 4] = time of observation
                [..., 5] = task during time of observation
                [..., 6] = block number during time of observation
                [..., 7] = scan cycle
                [..., 8] = slice number

        The data in this array is the basis of all model fits.
        """
        self.burn_in = burn_in

        observations = self.get_observations(include_background)

        # remove outlying scans due to severe movements of the subject
        # in the scanner

        if hasattr(self.reference_maps, 'outlying_scans'):
            scans = self.reference_maps.outlying_scans
            if scans.any():
                tmp = np.moveaxis(observations, self.ep+1, 1)
                tmp [scans] = np.nan

                if verbose:
                    print('{}: Removed {} ({:.2f}%) outlying scans'.format(
                        self.name.name(), scans.sum(), 100*scans.mean()))

        # these observations have no response as they are missing
        # (they may lie outside of the brain)
        none = np.isnan(observations[...,3])

        # these observations have no response as they are zero
        null = np.isclose(observations[...,3], 0)

        # put them together and set all observations to NAN when
        # anything is missing or invalid
        missings = none | null
        observations[missings] = np.nan

        if dropna:
            observations[ np.isnan(observations).any(axis=-1) ] = np.nan

        # scan cycles that we do not need to process
        if burn_in:
            observations[:burn_in] = np.nan

        # it is more numerically stable to work with a demeaned time
        # vector. This has also the consequence that the intercept will
        # refer to the mean signal intensity at the midpoint of the FMRI
        # session. Again, this is the point at which the intercept will
        # have the least variance.

        if demean:
            self.midpoint = np.nanmean(observations[...,4])
            observations[...,4] = observations[...,4] - self.midpoint

        valid = np.isfinite(observations).all(axis=-1)

        self.observations = observations
        self.data = observations[valid]

        self.dataframe = DataFrame({
            'x'      : self.data[...,0],
            'y'      : self.data[...,1],
            'z'      : self.data[...,2],
            'signal' : self.data[...,3],
            'time'   : self.data[...,4],
            'task'   : self.data[...,5],
            'block'  : self.data[...,6],
            'cycle'  : self.data[...,7],
            'slice'  : self.data[...,8]})

        if verbose:
            if demean:
                print("""{}:
            Number of within brain observations: {:>10,d}
            Number of    non brain observations: {:>10,d}
            Intercept field refers to time:      {:>10.2f} s""".format(
                self.name.name(), valid.sum(), (~valid).sum(),
                self.midpoint))
            else:
                print("""{}:
            Number of within brain observations: {:>10,d}
            Number of    non brain observations: {:>10,d}""".format(
                self.name.name(), valid.sum(), (~valid).sum()))

    def set_design(self, formula='C(task)/C(block, Sum)',
            parameter=['intercept', 'task'], verbose=True):
        """
        Create the design matrix

        Returns
        -------
        DesignMatrix
        """
        assert hasattr(self, 'dataframe'), 'first run .set_data()'

        dmat = dmatrix(formula, self.dataframe, eval_env=-1)
        names = dmat.design_info.column_names
        parameter_dict = { p : [p in n.lower() for n in names].index(True) for p in parameter}

        self.design  = np.asarray(dmat)
        self.formula = formula
        self.parameter_dict = parameter_dict
        return dmat

    ####################################################################
    # Fit at one coordinate
    ####################################################################

    def data_at_subject_coordinate(self, x):
        """
        Fit the signal model to data at specified coordinates given
        with respect to the coordinate system of the subject reference
        space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        """
        assert hasattr(self, 'scale'), 'first set hyperparameters'
        assert hasattr(self, 'radius'), 'first set hyperparameters'
        assert hasattr(self, 'data'), 'first set data'

        return data_at(coordinate=x,
                data=self.data,
                epi_code=self.epi_code,
                scale=self.scale,
                radius=self.radius)

    def data_at_index(self, index):
        """
        Fit the signal model to data at specified coordinates given
        with respect to the coordinate system of the subject reference
        space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        """
        x=self.population_map.diffeomorphism.apply_to_index(index)
        return self.data_at_subject_coordinate(x)

    def data_at_coordinate(self, x):
        """
        Fit the signal model to data at specified coordinates given
        with respect to the coordinate system of the subject reference
        space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        """
        x=self.population_map.diffeomorphism.apply(x)
        return self.data_at_subject_coordinate(x)

    def model_at_subject_coordinate(self, x,
            formula='signal~C(task)/C(block, Sum)', also_return_data=False):
        """
        Fit the signal model to data at specified coordinates given
        with respect to the coordinate system of the subject reference
        space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        formula : str
            A formula.
        also_return_data : bool
            If the formula does not contain time, the design matrix
            (design) will not contain time as a covariate. If
            also_return_data is True, this will also return the data
            matrix.
        """
        assert hasattr(self, 'scale'), 'first set hyperparameters'
        assert hasattr(self, 'radius'), 'first set hyperparameters'
        assert hasattr(self, 'data'), 'first set data'

        return model_at(formula=formula,
                coordinate=x,
                data=self.data,
                epi_code=self.epi_code,
                scale=self.scale,
                radius=self.radius)

    def model_at_index(self, index, **kwargs):
        """
        Fit the signal model to data at specified coordinates

        Parameters
        ----------
        index : tuple
            The index at which to fit the model
        """
        x=self.population_map.diffeomorphism.apply_to_index(index)
        return self.model_at_subject_coordinate(x, **kwargs)

    def model_at_coordinate(self, x, **kwargs):
        """
        Fit the signal model to data at specified coordinates given with
        respect to the coordinate system of the population space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        """
        x=self.population_map.diffeomorphism.apply(x)
        return self.model_at_subject_coordinate(x, **kwargs)

    def fit_at_subject_coordinate(self, x,
            formula='signal~C(task)/C(block, Sum)'):
        """
        Fit the signal model to data at specified coordinates given
        with respect to the coordinate system of the subject reference
        space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        formula : str
            A formula.
        also_return_data : bool
            If the formula does not contain time, the design matrix
            (design) will not contain time as a covariate. If
            also_return_data is True, this will also return the data
            matrix.
        """
        assert hasattr(self, 'scale'), 'first set hyperparameters'
        assert hasattr(self, 'radius'), 'first set hyperparameters'
        assert hasattr(self, 'data'), 'first set data'

        return fit_at(formula=formula,
                coordinate=x,
                epi_code=self.epi_code,
                data=self.data,
                scale=self.scale,
                radius=self.radius)

    def fit_at_index(self, index, **kwargs):
        """
        Fit the signal model to data at specified coordinates

        Parameters
        ----------
        index : tuple
            The index at which to fit the model
        """
        x=self.population_map.diffeomorphism.apply_to_index(index)
        return self.fit_at_subject_coordinate(x, **kwargs)

    def fit_at_coordinate(self, x, **kwargs):
        """
        Fit the signal model to data at specified coordinates given with
        respect to the coordinate system of the population space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        """
        x=self.population_map.diffeomorphism.apply(x)
        return self.fit_at_subject_coordinate(x, **kwargs)

    ####################################################################
    # Fit at many coordinates
    ####################################################################

    def get_mask(self, verbose=True):
        """
        Creates a data mask
        """
        coordinates = self.population_map.diffeomorphism.coordinates()
        to_index = self.mean_index_affine()
        idx = to_index.apply(coordinates)
        idx = idx.round().astype(int)

        mdx = idx.reshape((-1,3))

        di, dj, dk = self.session.shape

        for x in iter(mdx):
            i,j,k = x
            if i < 0 or j < 0 or k < 0 or i >= di or j >= dj or k >= dk:
                x[...] = 0
            elif np.isfinite(self.session.data[:,i,j,k]).any():
                x[...] = 1
            else:
                x[...] = 0

        mask = (idx > 0).any(-1)

        if verbose:
            print("""{}:
            Coordinates which are in data mask: {:>10,d}
            Coordinates which are not:          {:>10,d}""".format(
                self.name.name(), mask.sum(), (~mask).sum()))

        return mask

    def get_roi(self, mask=True, verbose=True):
        """
        Coordinates and mask

        Combines the data mask with other masks to produce a suitable
        mask for fitting.

        Parameters
        ----------
        mask : None or bool or str or ndarray, dtype: bool
            Either one of 'vb_mask', 'vb', 'vb_background',
            'foreground', or 'vb_estimate'. If None defaults to 'data'.
            True will take precedence: 'vb_mask' > 'vb' >
            'vb_background'> 'vb_estimate' > 'foreground'.
        verbose : bool
            Increase output verbosity

        Returns
        -------
        coordinates : ndarray, dtype: float
            Coordinates in standard space.
        mask : ndarray, dtype: bool
            Mask in standard space. True when within a ROI.
        """
        coordinates = self.population_map.diffeomorphism.coordinates()

        datamask = self.get_mask()

        if (mask is None) or (mask is False):
            mask = None
            maskname = 'no mask being applied '
        elif mask is True:
            try:
                mask = self.population_map.vb_mask.get_mask()
                maskname = 'template mask (vb_mask)'
            except AttributeError:
                try:
                    mask = self.population_map.vb.get_mask()
                    maskname = 'template (vb)'
                except AttributeError:
                    mask = datamask
                    maskname = 'data driven (foreground/background)'
        elif type(mask) is str:
            if mask == 'vb':
                mask = self.population_map.vb.get_mask()
                massname = 'template (vb)'
            elif mask == 'vb_background':
                mask = self.population_map.vb_background.get_mask()
                maskname = 'template background (vb_background)'
            elif mask == 'vb_estimate':
                mask = self.population_map.vb_estimate.get_mask()
                maskname = 'template estimate (vb_estimate)'
            else:
                mask = datamask
                maskname = 'data mask (foreground/background)'
        else:
            maskname = 'user defined'

        if mask is not None:
            assert type(mask) is np.ndarray, 'mask must be an ndarray'
            assert mask.dtype == bool, 'mask must be of dtype bool'
            assert mask.shape == coordinates.shape[:-1], 'mask shape must match image shape'
            mask = mask & datamask

        if verbose:
            print('{}: Fit is restricted to: {}'.format(self.name.name(), maskname))

        return coordinates, mask

    def fit_at_subject_coordinates(self, coordinates, mask=None,
            verbose=True, backend='numba'):
        """
        Fit the signal model to data

        Parameters
        ----------
        coordinates : None or ndarray, shape (…,3), dtype: float
            The coordinates at which to fit the model
        verbose : bool
            increase output verbosity

        Returns
        -------
        Result
        """
        if not hasattr(self, 'scale'):
            if verbose:
                print('{}! Set hyperparameters to default'.format(self.name.name()))
            self.set_hyperparameters()

        if not hasattr(self, 'radius'):
            if verbose:
                print('{}! Set hyperparameters to default'.format(self.name.name()))
            self.set_hyperparameters()

        if not hasattr(self, 'observations'):
            if verbose:
                print('{}! Set observations to default'.format(self.name.name()))
            self.set_data()

        if not hasattr(self, 'data'):
            if verbose:
                print('{}! Set design to default'.format(self.name.name()))
            self.set_design()

        if not hasattr(self, 'design'):
            if verbose:
                print('{}! Set design to default'.format(self.name.name()))
            self.set_design()

        if verbose:
            if mask is None:
                print("""{}:
            Number of coordinates to fit: {:>10,d}""".format(
                self.name.name(), mask.sum()))
            else:
                print("""{}:
            Number of coordinates to fit: {:>10,d}
            Number of coordinates not to: {:>10,d}""".format(
                self.name.name(), mask.sum(), (~mask).sum()))

        old_settings = np.seterr(divide='raise', invalid='raise')
        time0 = time.time()

        statistics, parameter_dict, value_dict = fit_field(
                coordinates = coordinates,
                mask        = mask,
                data        = self.data,
                design      = self.design,
                epi_code    = self.epi_code,
                scale       = self.scale,
                radius      = self.radius,
                verbose     = verbose,
                backend     = backend)

        time1 = time.time()
        np.seterr(**old_settings)

        if verbose:
            time_spend = time1 - time0
            print('{}: Time needed for the fit: {:.2f} min'.format(
                self.name.name(), time_spend / 60))
            print('{}: Time needed for the fit: {:.2f} h'  .format(
                self.name.name(), time_spend / 60**2))

        parameter_dict.update(self.parameter_dict)

        result = Result(
                coordinates     = coordinates,
                statistics      = statistics,
                population_map  = self.population_map,
                hyperparameters = self.hyperparameters(),
                parameter_dict  = parameter_dict,
                value_dict      = value_dict)

        return result

    def fit_at_indices(self, indices, **kwargs):
        coordinates = self.population_map.diffeomorphism.apply_to_indices(indices)
        return self.fit_at_subject_coordinates(coordinates = coordinates, **kwargs)

    def fit(self, mask=True, verbose=True, backend='numba'):
        """
        Fit the signal model to data

        Parameters
        ----------
        mask : None or bool or str or ndarray, dtype: bool
            string can be one of 'vb', 'vb_background',
            'foreground', or 'vb_estimate'. None defaults to 'data'. True will
            take precedence: 'vb'> 'vb_background'> 'vb_estimate' >
            'foreground'.
        verbose : bool
            increase output verbosity

        Returns
        -------
        Result : Fitted field.
        """
        coordinates, mask = self.get_roi(mask=mask, verbose=verbose)

        return self.fit_at_subject_coordinates(
                coordinates = coordinates,
                mask        = mask,
                backend     = backend)

    ###################################################################
    # Descriptive statistics of this session
    ###################################################################

    def hyperparameters(self):
        if hasattr(self, 'scale_type'):
            return {
                    'scale_type':self.scale_type,
                    'scale':self.scale,
                    'factor':self.factor,
                    'mass':self.mass,
                    'radius':self.radius,
                    }
        else:
            print('Hyperparameters not set')
            return

    def describe(self):
        description = """
        Hyperparameters:
            Scale type:      {:>6s}
            Factor:          {:>6.2f}

        Resulting in:
            Mass:            {:>5.4f}
            Scale:           {:>6.2f} mm
            FWHM:            {:>6.2f} mm
            Radius:          {:>6.2f} mm
            Diagonal:        {:>6.2f} mm
        """
        hyperparameters = self.hyperparameters()
        return description.format(
                hyperparameters['scale_type'],
                hyperparameters['factor'],
                hyperparameters['mass'],
                hyperparameters['scale'],
                2*np.sqrt(2*np.log(2)) * hyperparameters['scale'],
                hyperparameters['radius'],
                2*hyperparameters['radius'],
                )

    def save(self, file, **kwargs):
        """
        Save instance to disk

        This will save the current instance to disk for later use.

        Parameters
        ----------
        file : str
            File name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)

#######################################################################
#######################################################################
#
# Store the result of a fit to the data of an FMRI experiment
#
#######################################################################
#######################################################################

class Result:
    """
    Result of a FMRI fitting

    Defines a class for the result, i.e., the estimated effect field of
    an FMRI experiments, fitted by the function attribute of the `FMRI`
    class.

    Parameters
    ----------
    name : Identifier
        An identifier
    population_map : PopulationMap
    field : ndarray, shape (…,3), dtype: float
        The effect or fitted model parameter field
    coordinates : ndarray, shape (…,3), dtype: float
    hyperparameters : dict
    parameters : dict {str: int}
        The names of the parameters which have been fitted by the
        signal model in the order of appearance in the model.
    """
    def __init__(self, coordinates, statistics, population_map,
            hyperparameters, parameter_dict, value_dict):
        assert type(coordinates) is np.ndarray, 'coordinates must be a ndarray'
        assert type(population_map) is PopulationMap, 'population_map must be a PopulationMap'
        assert type(hyperparameters) is dict, 'hyperparameters must be a dict'
        assert type(parameter_dict) is dict, 'parameter_dict must be a dict'
        assert type(value_dict) is dict, 'value_dict must be a dict'

        self.coordinates = coordinates
        self.statistics = statistics
        self.population_map = population_map
        self.hyperparameters = hyperparameters
        self.parameter_dict = parameter_dict
        self.value_dict = value_dict

        self.name = self.population_map.diffeomorphism.nb

    ####################################################################
    # Norm to ATI
    ####################################################################

    def norm_to(self, rf):
        """
        Norm the part of the statistics field that has a unit, to a
        reference field.

        Parameters
        ----------
        rf : Image
            The reference field for the unit.
        """

        has_unit = self.statistics[...,:2,:]
        self.statistics[...,:2,:] = self.statistics[...,:2,:] * \
                (rf.data / self.statistics[...,0,0])[...,None,None]

    def norm_to_ati(self):
        """
        Norm the part of the statistics field that has a unit, to the
        ati-reference field that is stored in the PopulationMap.
        """
        self.norm_to(self.population_map.vb_ati)

    ####################################################################
    # Extract summary statistics
    ####################################################################

    def get_field(self, param, value=None):
        """
        Extract scalar field

        Parameters
        ----------
        param : str
            Check the attribute parameter_dict for possible values of
            `param`.
        value : str
            Check the attribute value_dict for possible values of
            `value`. Note that 'tstatistic' is also valid and will be
            calculated from the key point and stderr.

        Returns
        -------
        Image
            The queried field

        Notes
        -----
        It only makes sense to combine intercept or activation
        with point, stderr, and tstatistic, and it only makes sense to
        combine other with mse (an estimate of sigma-squared) or df
        (residual degrees of freedom)
        """
        if value is None:
            value = param

        field = extract_field(
                field=self.statistics,
                param=param,
                value=value,
                parameter_dict=self.parameter_dict,
                value_dict=self.value_dict)

        return Image(
                reference=self.population_map.diffeomorphism.reference,
                data=field, name=self.name)

    def volume(self):
        return np.isfinite(self.statistics).sum() * \
                self.population_map.diffeomorphism.reference.volume()

    ###################################################################
    # Descriptive statistics of this session
    ###################################################################

    def mask(self, mask=True, verbose=False):
        """
        Apply mask to parameter fields
        """
        if (mask is None) or (mask is False):
            mask = None
            maskname = 'no mask is being applied '
        elif mask is True:
            try:
                mask0 = self.population_map.vb_mask.get_mask()
                maskname0 = 'template mask (vb_mask)'
            except AttributeError:
                mask0 = True
                maskname0 = 'no vb to apply as mask'

            try:
                mask1 = self.population_map.vb.get_mask()
                maskname1 = 'template (vb)'
            except AttributeError:
                mask1 = True
                maskname1 = 'no vb mask to apply'

            mask = mask0 & mask1
            maskname = maskname0 + ' and ' + maskname1

        elif type(mask) is str:
            if mask == 'vb':
                mask = self.population_map.vb.get_mask()
                massname = 'template (vb)'
            elif mask == 'vb_background':
                mask = self.population_map.vb_background.get_mask()
                maskname = 'template background (vb_background)'
            elif mask == 'vb_estimate':
                mask = self.population_map.vb_estimate.get_mask()
                maskname = 'template estimate (vb_estimate)'
            elif mask == 'vb_mask':
                mask = self.population_map.vb_mask.get_mask()
                maskname = 'template mask (vb_mask)'
            else:
                mask = None
                maskname = 'no mask to apply'
        else:
            maskname = 'user defined'

        if mask is not None:
            assert type(mask) is np.ndarray, 'mask must be an ndarray'
            assert mask.dtype == bool, 'mask must be of dtype bool'
            assert mask.shape == self.statistics.shape[:-2], 'mask shape must match image shape'

        if verbose:
            print('Statistics field is restricted to: {}'.format(maskname))

        if mask is not None:
            assert type(mask) is np.ndarray, 'mask must be an ndarray'
            assert mask.dtype == bool, 'mask must be of dtype bool'
            assert mask.shape == self.statistics.shape[:-2], 'mask shape must match image shape'
            self.statistics [ ~mask ] = np.nan

    def descriptive_statistics(self):
        mask = np.isnan(self.statistics).any(axis=(-1,-2))
        nooc = np.isnan(mask).all(axis=(0,1))
        return self.volume(), (~nooc).sum()

    def describe(self):
        description = """
        Cohort:   {}
        Subject:  {}
        Paradigm: {}

        Hyperparameters:
            Scale type:      {:>6s}
            Factor:          {:>6.2f}

        Resulting in:
            Mass:            {:>5.4f}
            Scale:           {:>6.2f} mm
            FWHM:            {:>6.2f} mm
            Radius:          {:>6.2f} mm
            Diagonal:        {:>6.2f} mm

        Fitted statistics field:
        Shape:   {}
        Volume:  {:.2f} mm^3"""
        return description.format(
                self.population_map.diffeomorphism.nb.cohort,
                self.population_map.diffeomorphism.nb.j,
                self.population_map.diffeomorphism.nb.paradigm,

                self.hyperparameters['scale_type'],
                self.hyperparameters['factor'],

                self.hyperparameters['mass'],
                self.hyperparameters['scale'],
                2*np.sqrt(2*np.log(2)) * self.hyperparameters['scale'],
                self.hyperparameters['radius'],
                2*self.hyperparameters['radius'],

                self.statistics.shape,
                self.volume())

    def __str__(self):
        return self.describe()

    #######################################################################
    # Save instance to disk
    #######################################################################

    def save(self, file, **kwargs):
        """
        Save instance to disk

        This will save the current instance to disk for later use.

        Parameters
        ----------
        file : str
            File name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)
