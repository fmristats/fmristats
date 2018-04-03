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

from .irritation import Irritation

from .fit import fit_field, extract_field, model_at, data_at

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
        self.irritation = session.irritation
        self.reference_maps = reference_maps
        self.population_map = population_map

        self.shape = (session.numob, session.shape[session.ep])
        self.slice_timing_design = session.slice_time

        self.name = self.session.name
        self.epi_code = self.session.epi_code
        self.ep = abs(self.epi_code) - 1

        # Affine maps from the index space of the Session to
        # the reference space of the ReferenceMaps.

        # TODO: when having added the option: --reference-maps None,
        # replace this with something reasonable
        self.references = reference_maps.scan_references.inv().dot(session.reference)

        self.tsr = self.session.temporal_resolution / self.shape[1]

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
        The affine transformation that maps a coordinate in reference
        space to the mean index in the index space of the acquisition
        grid of the session.
        """
        return self.references.inv().mean()

    #######################################################################
    # All about the observations
    #######################################################################

    def coordinates(self):
        """
        The coordinate grid of observations

        The coordinates in subject reference space of all points at
        which MR signals have been aquired.

        Returns
        -------
        ndarray, shape (n,x,y,z,3), dtype: float
        """
        n,x,y,z = self.session.data.shape
        indices = ((slice(0,x), slice(0,y), slice(0,z)))

        coordinates = np.empty((n,x,y,z,3))
        for t in range(n):
            coordinates[t] = self.references.index(t).apply_to_indices(indices)

        return coordinates

    def create_slice_timing_design(self):
        return self.slice_timing_design

    def create_irritation_design(self, **kwargs):
        assert hasattr(self, 'slice_timing_design'), 'first create the slice timing design'
        self.irritation_design = self.irritation.design(
                slice_timing=self.slice_timing_design, **kwargs)
        return self.irritation_design

    def observations(self):
        """
        Create observation matrix

        Returns
        -------
        ndarray, shape (n,x,y,z,6), dtype: float
            [...,:3] = coordinates of observation
            [..., 3] = MR signal response
            [..., 4] = time of observation
            [..., 5] = task during time of observation
            [..., 6] = block number during time of observation
            [..., 7] = scan cycle
            [..., 8] = slice number
        """
        assert hasattr(self, 'slice_timing_design'), 'first set the slice time design'
        assert hasattr(self, 'irritation_design'), 'first set the irritation design'
        n,x,y,z = self.session.data.shape

        observations = np.zeros((n,x,y,z,9))
        observations[...,:3] = self.coordinates()
        observations[..., 3] = self.session.data

        tmp = np.moveaxis(observations, self.ep+1, 1)
        tmp[...,4]  = self.slice_timing_design[...,None,None]
        tmp[...,5:7] = self.irritation_design[...,None,None,:]

        x = np.mgrid[:self.shape[0], :self.shape[1]]
        x = np.moveaxis(x, 0, -1)
        tmp[...,7:9] = x[...,None,None,:]

        return observations

    def create_data_matrix(self, burn_in=None, verbose=True):
        """
        Create the data matrix

        This is the essentially the same matrix as returned by
        .observations(). However, scan cycles with no subject
        irritation or missing data are removed. However, time is centred
        to its mean.

        Returns
        -------
        ndarray, shape (n',x,y,z,6), dtype: float
            [...,:3] = coordinates of observation
            [..., 3] = MR signal response
            [..., 4] = time of observation
            [..., 5] = task during time of observation
            [..., 6] = block number during time of observation
        """
        self.burn_in = burn_in

        observations = self.observations()

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

        # these observations have either no response (e.g. they are
        # outside of the brain) or there are not all covariates
        # available (e.g. they have been acquired outside of task blocks
        # or the position of the slice in space could not be detected
        # with success.
        none = np.isnan(observations).any(axis=-1)

        # these observations have no response as they are identical
        # zero.
        null = np.isclose(observations[...,3], 0)

        # put them together and set all observations to NAN when
        # anything is missing or invalid
        missings = none | null
        observations[missings] = np.nan

        # scan cycles that we do not need to process at all
        if burn_in:
            observations[:burn_in] = np.nan

        # it is more numerically stable to work with a demeaned time
        # vector. This has also the consequence that the intercept will
        # refer to the mean signal intensity at the midpoint of the FMRI
        # session. Again, this is the point at which the intercept will
        # have the least variance.
        self.midpoint = np.nanmean(observations[...,4])
        observations[...,4] = observations[...,4] - self.midpoint

        if verbose:
            valid = np.isfinite(observations).all(axis=-1)
            print("""{}:
            …voxels which are part of the model: {:>10,d}
            …voxels which are not:               {:>10,d}
            …intercept field refers to:          {:>10.2f} s
            """.format(self.name.name(), valid.sum(), (~valid).sum(), self.midpoint))

        data_matrix = DataFrame({
            'i'     : observations[...,0].ravel(),
            'j'     : observations[...,1].ravel(),
            'k'     : observations[...,2].ravel(),
            'y'     : observations[...,3].ravel(),
            'time'  : observations[...,4].ravel(),
            'task'  : observations[...,5].ravel(),
            'block' : observations[...,6].ravel(),
            'cycle' : observations[...,7].ravel(),
            'slice' : observations[...,8].ravel()})

        data_matrix.dropna(inplace=True)
        self.data_matrix = data_matrix
        return self.data_matrix

    def create_design_matrix(self, formula='C(task)/C(block, Sum)',
            parameter=['intercept', 'task'], verbose=True):
        """
        Create the design matrix

        Returns
        -------
        (DesignMatrix, DesignMatrix)
        """
        assert hasattr(self, 'data_matrix'), 'first set the design using .create_data_matrix()'

        dmat = dmatrix(formula, self.data_matrix, eval_env=-1)
        names = dmat.design_info.column_names
        parameter_dict = { p : [p in n.lower() for n in names].index(True) for p in parameter}

        self.exog  = np.asarray(dmat)
        self.formula = formula
        self.parameter_dict = parameter_dict
        return dmat

    ####################################################################
    # Fit at one coordinate
    ####################################################################

    def data_at_subject_coordinate(self, x):
        """
        Fit the signal model to data at specificed coordinates given
        with respect to the coordinate system of the subject reference
        space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        """
        assert hasattr(self, 'scale'), 'first set hyperparameters'
        assert hasattr(self, 'radius'), 'first set hyperparameters'
        assert hasattr(self, 'data_matrix'), 'first set data_matrix'

        return data_at(coordinate=x,
                data=self.data_matrix,
                scale=self.scale,
                radius=self.radius)

    def data_at_index(self, index):
        """
        Fit the signal model to data at specificed coordinates given
        with respect to the coordinate system of the subject reference
        space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        """
        x=self.population_map.diffeomorphism.apply_to_index(index)
        return self.data_at_subject_coordinate(x)

    def model_at_subject_coordinate(self, x, formula='y~C(task)/C(block, Sum)'):
        """
        Fit the signal model to data at specificed coordinates given
        with respect to the coordinate system of the subject reference
        space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        """
        assert hasattr(self, 'scale'), 'first set hyperparameters'
        assert hasattr(self, 'radius'), 'first set hyperparameters'
        assert hasattr(self, 'data_matrix'), 'first set data_matrix'

        return model_at(coordinate=x, formula=formula,
                data=self.data_matrix,
                scale=self.scale,
                radius=self.radius)

    def model_at_index(self, index, **kwargs):
        """
        Fit the signal model to data at specificed coordinates

        Parameters
        ----------
        index : tuple
            The index at which to fit the model
        """
        x=self.population_map.diffeomorphism.apply_to_index(index)
        return self.model_at_subject_coordinate(x, **kwargs)

    def model_at_coordinate(self, x, **kwargs):
        """
        Fit the signal model to data at specificed coordinates given
        with respect to the coordinate system of the population space.

        Parameters
        ----------
        x : ndarray, shape (3,), dtype: float
            The coordinates at which to fit the model
        """
        x=self.population_map.diffeomorphism.apply(x)
        return self.model_at_subject_coordinate(x, **kwargs)

    ####################################################################
    # Fit at many coordinates
    ####################################################################

    def get_data_mask(self, verbose=True):
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
            …coordinates which are in data mask: {:>10,d}
            …coordinates which are not:          {:>10,d}""".format(
                self.name.name(), mask.sum(), (~mask).sum()))

        return mask

    def create_estimation_matrix(self, mask=True, verbose=True):
        """
        Parameters
        ----------
        mask : None or bool or str or ndarray, dtype: bool
            string can be one of 'template', 'template_mask', or 'target'
        verbose : bool
            increase output verbosity

        Returns
        -------
        ndarray
        """
        coordinates = self.population_map.diffeomorphism.coordinates()

        if (mask is False) or (mask is None):
            mask = None

        elif mask is True:
            if hasattr(self.population_map, 'template_mask'):
                mask = self.population_map.template_mask.get_mask()
                maskname = 'template mask'
            elif hasattr(self.population_map, 'template'):
                mask = self.population_map.template.get_mask()
                maskname = 'template'
            else:
                mask = None
                maskname = 'data mask'

        elif type(mask) is str:
            if mask == 'template_mask':
                mask = self.population_map.template_mask.get_mask()
                massname = 'template mask'
            elif mask == 'template':
                mask = self.population_map.template.get_mask()
                maskname = 'template'
            elif mask == 'target':
                mask = self.population_map.target.get_mask()
                masskname = 'target'
            else:
                mask = None
                maskname = 'data mask'

        else:
            maskname = 'user defined'

        if mask is not None:
            assert type(mask) is np.ndarray, 'mask must be an ndarray'
            assert mask.dtype == bool, 'mask must be of dtype bool'
            assert mask.shape == coordinates.shape[:-1], 'mask shape must match image shape'
            mask = mask & self.get_data_mask()
        else:
            mask = self.get_data_mask()

        if verbose:
            print('{}: Fit is restricted to: {}'.format(self.name.name(), maskname))

        self.estimation_matrix = coordinates
        self.estimation_mask   = mask
        return self.estimation_matrix, self.estimation_mask

    def fit_at_subject_coordinates(self, coordinates, mask=None, verbose=True):
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
        assert hasattr(self, 'scale'), 'first set hyperparameters using .set_hyperparameters()'
        assert hasattr(self, 'radius'), 'first set hyperparameters using .set_hyperparameters()'
        assert hasattr(self, 'data_matrix'), 'first set the design using .create_data_matrix()'
        assert hasattr(self, 'exog'), 'first set the design using .create_design_matrix()'

        if verbose:
            print('{}: Start fit'.format(self.name.name()))

        old_settings = np.seterr(divide='raise', invalid='raise')
        time0 = time.time()

        statistics, parameter_dict, value_dict = fit_field(
                coordinates = coordinates,
                mask        = mask,
                endog       = self.data_matrix.y.values,
                exog        = self.exog,
                agc         = self.data_matrix[['i','j','k']].values,
                scale       = self.scale,
                radius      = self.radius,
                verbose     = verbose,
                name        = self.name.name())

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
        coordinates = self.diffeomorphism.apply_to_indices(indices)
        return self.fit_at_subject_coordinates(coordinates = coordinates, **kwargs)

    def fit(self, mask=True, verbose=True):
        """
        Fit the signal model to data

        Parameters
        ----------
        mask : None or bool or str or ndarray, dtype: bool
            string can be one of 'template', 'template_mask', or 'target'
        verbose : bool
            increase output verbosity

        Returns
        -------
        Result : Fitted field.
        """
        if not hasattr(self, 'scale'):
            if verbose:
                print('{}! Set hyperparameters to default'.format(self.name.name()))
            self.set_hyperparameters()

        if not hasattr(self, 'radius'):
            if verbose:
                print('{}! Set hyperparameters to default'.format(self.name.name()))
            self.set_hyperparameters()

        if not hasattr(self, 'data_matrix'):
            if verbose:
                print('{}! Set data matrix to default (no burn in)'.format(self.name.name()))
            self.create_data_matrix()

        if not hasattr(self, 'exog'):
            if verbose:
                print('{}! Set design matrix / model to default'.format(self.name.name()))
            self.create_design_matrix()

        if verbose:
            print('{}: Create estimation matrix and mask'.format(self.name.name()))
        self.create_estimation_matrix(mask=mask, verbose=verbose)

        return self.fit_at_subject_coordinates(
                coordinates = self.estimation_matrix,
                mask        = self.estimation_mask)


    ###################################################################
    # Descriptive statistics of this session
    ###################################################################

    def hyperparameters(self):
        return {
                'scale_type':self.scale_type,
                'scale':self.scale,
                'factor':self.factor,
                'mass':self.mass,
                'radius':self.radius,
                #'half_edge':self.half_edge,
                }

    def describe(self):
        description = """
        Hyperparameters were set to:
            Scale type:   {:s}
            Scale:        {:.2f} mm
            Factor:       {:.2f}
        Resulting in:
            Mass:         {:3>,.5f}
            FWHM:         {:3>,.2f} mm
            Radius:       {:3>,.2f} mm
            Diagonal:     {:3>,.2f} mm"""
        hyperparameters = self.hyperparameters()
        return description.format(
                hyperparameters['scale_type'],
                hyperparameters['scale'],
                hyperparameters['factor'],

                hyperparameters['mass'],
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
        The effect field, i.e., result of fitting the RSM model.
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

        self.name = self.population_map.name

    ####################################################################
    # Extract summary statistics
    ####################################################################

    def get_field(self, param, value='other'):
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
        #if value == 'tstatistic':
        #    pfield = extract_field(
        #            field=self.statistics,
        #            param=param,
        #            value='point',
        #            parameter_dict=self.parameter_dict,
        #            value_dict=self.value_dict)

        #    sfield = extract_field(
        #            field=self.statistics,
        #            param=param,
        #            value='stderr',
        #            parameter_dict=self.parameter_dict,
        #            value_dict=self.value_dict)

        #    field = pfield / sfield

        #else:
        #    field = extract_field(
        #            field=self.statistics,
        #            param=param,
        #            value=value,
        #            parameter_dict=self.parameter_dict,
        #            value_dict=self.value_dict)

        field = extract_field(
                field=self.statistics,
                param=param,
                value=value,
                parameter_dict=self.parameter_dict,
                value_dict=self.value_dict)

        return Image(
                reference=self.population_map.diffeomorphism.reference,
                data=field,
                name='-'.join((param,value)))

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
        if (mask is False) or (mask is None):
            mask = None

        elif mask is True:
            if hasattr(self.population_map, 'template_mask'):
                mask = self.population_map.template_mask.get_mask()
                maskname = 'template mask'
            elif hasattr(self.population_map, 'template'):
                mask = self.population_map.template.get_mask()
                maskname = 'template'
            else:
                mask = None
                maskname = 'none'

        elif type(mask) is str:
            if mask == 'template_mask':
                mask = self.population_map.template_mask.get_mask()
                massname = 'template mask'
            elif mask == 'template':
                mask = self.population_map.template.get_mask()
                maskname = 'template'
            elif mask == 'target':
                mask = self.population_map.target.get_mask()
                masskname = 'target'
            else:
                mask = None
                maskname = 'none'

        else:
            maskname = 'user defined'

        if verbose:
            print('{}: Fit is restricted to: {}'.format(self.name.name(), maskname))

        if mask is not None:
            assert type(mask) is np.ndarray, 'mask must be an ndarray'
            assert mask.dtype == bool, 'mask must be of dtype bool'
            assert mask.shape == self.statistics.shape[:-2], 'mask shape must match image shape'
            self.statistics [ ~mask ] = np.nan

    def summary(self, quantile_t=None):
        tf = self.get_field('activation', 'tstatistic')

        if quantile_t:
            assert quantile_t <= 1, 'quantile_t must be less than 1'
            assert quantile_t >= 0, 'quantile_t must be greater than 0'
            threshold_t = np.nanpercentile(tf.data, q=quantile_t*100)
            neuronal_active = tf.data > threshold_t
        else:
            neuronal_active = np.isfinite(tf.data)

        intercept  = self.parameter_dict ['intercept']
        activation = self.parameter_dict ['activation']
        degfreedom = self.parameter_dict ['df']

        point      = self.value_dict ['point']
        stderr     = self.value_dict ['stderr']
        tstatistic = self.value_dict ['tstatistic']
        other      = self.value_dict ['other']

        statistics = self.statistics [neuronal_active]

        return DataFrame({
            'cohort'   : self.name.cohort,
            'id'       : self.name.j,
            'paradigm' : self.name.paradigm,
            'date'     : self.name.datetime,
            'x' : self.coordinates[neuronal_active,0],
            'y' : self.coordinates[neuronal_active,1],
            'z' : self.coordinates[neuronal_active,2],
            'intercept'  : statistics[...,point,intercept],
            'activation' : statistics[...,point,activation],
            'stderr'     : statistics[...,stderr,activation],
            'tstatistic' : statistics[...,tstatistic,activation],
            'df' : statistics[...,other,degfreedom],
            })

    def descriptive_statistics(self):
        mask = np.isnan(self.statistics).any(axis=(-1,-2))
        nooc = np.isnan(mask).all(axis=(0,1))
        return self.volume(), (~nooc).sum()

    def describe(self):
        description = """
        Cohort:   {}
        Subject:  {}
        Paradigm: {}
        Hyperparameters were set to:
            Scale type:      {:s}
            Factor:          {:.2f}
        Resulting in:
            Scale:           {:.2f} mm
            Mass:            {:.5f}
            FWHM:            {:.2f} mm
            Radius:          {:.2f} mm
            Diagonal:        {:.2f} mm
        Shape:               {}
        Volume:              {:.2f} mm^3
        """
        return description.format(
                self.population_map.diffeomorphism.nb.cohort,
                self.population_map.diffeomorphism.nb.j,
                self.population_map.diffeomorphism.nb.paradigm,

                self.hyperparameters['scale_type'],
                self.hyperparameters['factor'],

                self.hyperparameters['scale'],
                self.hyperparameters['mass'],
                2*np.sqrt(2*np.log(2)) * self.hyperparameters['scale'],
                self.hyperparameters['radius'],
                2*self.hyperparameters['radius'],

                self.statistics.shape,
                self.volume())

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
