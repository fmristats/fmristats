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

Defines classes for diffeomorphisms which are not necessarily affine
transformations.

"""

from .name import Identifier

from .affines import Affine

import pickle

import numpy as np

import numpy.ma as ma

from numpy.linalg import inv

from scipy.ndimage import label, map_coordinates, \
        maximum_filter, generate_binary_structure, binary_erosion, \
        maximum_position

class Diffeomorphism:
    """
    A diffeomorphism ψ from standard space to subject space

    Parameters
    ----------
    reference : Affine or ndarray, shape (4,4)
        An affine transformation
    shape : tuple
        Shape of template data in standard space.
    vb : None or str or Identifier
        An object that will be used to identify the domain (vb) of this
        diffeomorphism.
    nb : None or Identifier
        An object that will be used to identify the image (nb) of this
        diffeomorphism.
    name : None or str
        Method that has produced this diffeomorphism.
    metadata : dict
        If meta data is provided which give information on which data
        had been used for the fit: which had provided, should be a dict
        with at least the following fields: ``{'vb_file':
        path/which/defined/the/domain, 'nb_file':
        path/which/defined/the/image,}``.

    Notes
    -----
    A reasonable subclass of Diffeomorphism must define the attributes
    :func:`apply_to_index`, :func:`apply_to_indices`, and :func:`apply`.
    """
    def __init__(self, reference, shape, vb=None, nb=None, name=None, metadata=None):
        self.reference = reference
        self.shape = shape
        self.vb = vb
        self.nb = nb
        self.name = nb.name

        if metadata is not None:
            self.metadata = metadata

    def apply_to_index(self, index):
        """
        Apply diffeomorphism to the point with given index

        Parameters
        ----------
        index : tuple
            Should be a tuple of integers

        Returns
        -------
        ndarray
        """
        pass

    def apply_to_indices(self, indices):
        """
        Apply diffeomorphism to the points with given indices

        Parameters
        ----------
        index : tuple
            Should be a tuple of slices or ellipses

        Returns
        -------
        ndarray
        """
        pass

    def apply(self, coordinates):
        """
        Apply diffeomorphism to the point at given coordinate

        This is the identity; it will return the same coordinates

        Parameters
        ----------
        index : tuple
            Should be a tuple of integers

        Returns
        -------
        ndarray
        """
        pass

    def coordinates(self):
        """
        Coordinates in the image

        Returns
        -------
        ndarray
        """
        x,y,z = self.shape
        indices = ((slice(0,x), slice(0,y), slice(0,z)))
        return self.apply_to_indices(indices)

    def coordinates_domain(self):
        """
        Coordinates in the domain

        Returns
        -------
        ndarray
        """
        x,y,z = self.shape
        indices = ((slice(0,x), slice(0,y), slice(0,z)))
        return self.reference.apply_to_indices(indices)

    def describe(self):
        if type(self.nb) is Identifier:
            describe_subject_space = """
        Subject reference space
        -----------------------
        Cohort:   {}
        Subject:  {}
        Paradigm: {}""".format(
                self.nb.cohort,
                self.nb.j,
                self.nb.paradigm,
                )
        else:
            describe_subject_space = """
        Subject reference space
        -----------------------
        Name:   {}""".format(self.nb)

        description = """
        Population space
        ----------------
        Name:  {}
        Shape: {}
        {:s}
        {:s}"""

        return description.format(
                self.vb,
                self.shape,
                self.reference.describe(),
                describe_subject_space)

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

    def __str__(self):
        return self.describe()

class Identity(Diffeomorphism):
    """
    The identity diffeomorphism ψ.
    """
    def __init__(self, reference, shape, vb=None, nb=None, name=None):
        self.reference = reference
        self.shape = shape
        self.vb = vb
        self.nb = nb
        self.name = name

    def apply_to_index(self, index):
        """
        Apply diffeomorphism to the point with given index

        Parameters
        ----------
        index : tuple
            Should be a tuple of integers.

        Returns
        -------
        ndarray
        """
        return self.reference.apply_to_index(index)

    def apply_to_indices(self, indices):
        """
        Apply diffeomorphism to the points with given indices

        Parameters
        ----------
        index : tuple
            Should be a tuple of slices or ellipses.

        Returns
        -------
        ndarray
        """
        return self.reference.apply_to_indices(indices)

    def apply(self, coordinates):
        """
        Apply diffeomorphism to the point at given coordinate

        This is the identity; it will return the same coordinates

        Parameters
        ----------
        index : tuple
            Should be a tuple of integers.

        Returns
        -------
        ndarray
        """
        return coordinates

# TODO: why not treat an Image as a masked array with values 0 and nan.

class Image(Identity):
    """
    Holds all information of a 3D-image.

    Parameters
    ----------
    reference : Affine or ndarray, shape (4,4)
        An affine transformation that maps indices in the data
        array to coordinates in Euclidean space.
    data : np.ndarray
        The data array.
    name : str or Identifier
        An identifier or name for this image.
    """
    def __init__(self, reference, data, name=None):
        if type(reference) is Affine:
            self.reference = reference
        else:
            self.reference = Affine(reference)

        self.data = data
        self.shape = self.data.shape

        if name:
            self.name = name
            self.vb = name
            self.nb = name
        else:
            self.name = None
            self.vb = None
            self.nb = None

    def ravel(self):
        """
        Return a contiguous flattened array of the image

        A 1-D array, containing the non-null, and non-nan elements of
        the input, is returned. Array is a copy of the data.
        """
        data = self.mask(inplace=False).data
        data = data[np.isfinite(data)]
        return data.ravel()

    def flatten(self, epi_code, xpic, ypic, slices=None):
        """
        Flatten a 3D image into a 2D image

        Parameters
        ----------
        epi_code : int

        Returns
        -------
        ndarray

        Notes
        -----
        Flattens by slicing the 3D image along the axis of the respected
        EPI code, and stack the slices into a xpic x ypic matrix
        """
        assert epi_code <=  3, 'epi code must be ≤ 3'
        assert epi_code >= -3, 'epi code must be ≥-3'
        assert epi_code !=  0, 'epi code must be different from 0'

        npic = xpic * ypic
        data = np.moveaxis(self.data,abs(epi_code)-1,0)
        filled, = np.where(~np.isnan(data).all(axis=(1,2)))
        if slices is None:
            select = np.linspace(0,len(filled)-1,npic,dtype='int')
            slices = filled[select]
        else:
            assert len(slices) == npic, 'number of slices must equal xpic*ypic'

        data = data[slices]

        imatrix = np.vstack(
                [np.hstack([data[i].T for i in range(j,xpic+j)])
                    for j in range(0,npic,xpic)])

        xticks = [i*data.shape[2] for i in range(xpic)]
        yticks = [i*data.shape[1] for i in range(ypic)]
        return imatrix, xticks, yticks, slices

    def mask(self, mask=None, inplace=True):
        """
        Mask the data array by mask

        Parameters
        ----------
        mask : ndarray, dtype(bool)
            An array of bool indicating the foreground

        Returns
        -------
        Image
            The same image but with data set to nan outside the mask
        """
        data = self.data.astype(float)

        if mask is None:
            mask = self.get_mask()

        data [~mask] = np.nan

        if inplace:
            self.data = data
        else:
            return Image(reference=self.reference, data=data, name=self.name)

    def get_mask(self, **kwargs):
        """
        Create a mask

        Returns
        -------
        ndarray, dtype(bool)
            boolean array which is *True* whenever values in the data are
            *valid*.
        """
        return ~( (np.isnan(self.data)) | np.isclose(self.data, 0, **kwargs) )

    def round(self):
        """
        Round data in the image to integer

        Returns
        -------
        Image
        """
        data = self.data.copy()
        data [ ~self.get_mask() ] = 0
        data = data.round().astype(int)
        return Image(reference=self.reference, data=data, name=self.name)

    def volume(self):
        """
        Calculates the volume occupied by the image

        Returns
        -------
        float : volume of non-zero values in the image
        """
        return self.get_mask().sum() * self.reference.volume()

    def detect_peaks(self, connectivity=3, **kwargs):
        """
        Detects the peaks in an image

        Takes an image and detect the peaks using a local maximum
        filter. Returns a boolean mask of the peaks (i.e. 1 when the
        pixel's value is the neighbourhood of a maximum and 0 otherwise)

        Parameters
        ----------


        """
        if self.data is None:
            return

        image = self.data

        # define an connected neighbourhood
        neighbourhood = generate_binary_structure(3,
                connectivity=connectivity, **kwargs)

        # apply the local maximum filter; all pixel of maximal value
        # in their neighbourhood are set to 1:
        local_max = maximum_filter(image, footprint=neighbourhood)==image

        # local_max is a mask that contains the peaks we are looking
        # for, but also the background. In order to isolate the peaks we
        # must remove the background from the mask.

        # we create the mask of the background
        background = (image==0) | np.isclose(image, 0)

        # a little technicality: we must erode the background in order
        # to successfully subtract it form local_max, otherwise a line
        # will appear along the background border (artefact of the local
        # maximum filter).

        eroded_background = binary_erosion(background,
                structure=neighbourhood, border_value=1)

        # we obtain the final mask, containing only peaks, by removing
        # the background from the local_max mask (xor operation)
        detected_peaks = local_max ^ eroded_background

        return Image(self.reference, detected_peaks)

    def components(self, threshold, component_size=None, direction='lower'):
        """
        Components of an image

        Enumerates (more precisely: labels) the components above or below
        the given threshold in the image.

        Parameters
        ----------
        threshold : float
            The threshold
        component_size : float
            Minimal volume for a component to be labelled
        direction : str, 'lower' or 'upper'
            Direction of the cut. Default is 'lower'.

        Returns
        -------
        labelled : Image
            The labelled image
        labels : list
            The labels in the image

        Notes
        -----
        Please note that the `component_size` is expected to be given in
        mm^3 and **not** as the number of voxels.
        """
        if direction == 'upper':
            label_map, n_labels = label(self.data < threshold)
        else:
            label_map, n_labels = label(self.data > threshold)

        if component_size is not None:
            for l in range(1,n_labels+1):
                x = label_map == l
                if self.reference.volume() * x.sum() < component_size:
                    label_map[x] = 0

        label_image = Image(self.reference, label_map)
        labels = np.unique(label_image.data)
        labels = labels[labels != 0]
        return label_image, labels

    def mean(self):
        """
        Calculates the mean value in the image

        Notes
        -----
        Values of 0 or nan will not contribute to the mean signal.  If
        your image ``img`` contains 0 values, then ``im.data.mean()``
        and ``im.mean()`` will differ and this is intended.

        Returns
        -------
        float : mean value of the image
        """
        data = ma.array(self.data, mask=~self.get_mask(), fill_value=0)
        return data.mean()

    def range(self):
        """
        Return the range of values in the image

        Notes
        -----
        Values of 0 or nan will not contribute to the range.

        Returns
        -------
        tuple : (minimum, maximum)
        """
        return np.nanmin(self.data), np.nanmax(self.data)

    def describe(self):
        """
        Describe the image

        Returns
        -------
        str
        """
        description = """
        Name: {}
        Shape: {}
        Volume of one pixel: {:.2f}
        Volume of image: {:.2f}
        """
        return description.format(
                self.name,
                self.shape,
                self.reference.volume(),
                self.volume())

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

class AffineTransformation(Diffeomorphism):
    """
    An affine transformation ψ mapping from a standard space :math:`M`
    to a subject reference space :math:`R`. It is

    .. math::

        ψ (x) = A⋅x.

    for :math:`x∈M` and :math:`A` an affine transformation.

    A point [i,j,k] in the index space of the template in :math:`M` maps
    to the coordinates

    .. math::

        ψ(reference[i,j,k]) = A ⋅ reference[i,j,k]

    in the subject reference space.

    Parameters
    ----------
    reference : Affine or ndarray, shape (4,4)
        An affine transformation
    affine : Affine or ndarray, shape (4,4)
        An affine transformation
    """

    def __init__(self, reference, affine, shape, vb=None, nb=None, name=None):
        if type(reference) is Affine:
            self.reference = reference
        else:
            self.reference = Affine(reference)

        if type(affine) is Affine:
            self.affine = affine
        else:
            self.affine = Affine(affine)

        self.shape = shape
        self.vb = vb
        self.nb = nb
        self.name = name

    def apply_to_index(self, index):
        return self.affine.dot(self.reference).apply_to_index(index)

    def apply_to_indices(self, indices):
        return self.affine.dot(self.reference).apply_to_indices(indices)

    def apply(self, coordinates):
        """
        Apply diffeomorphism to the point at given coordinate

        Parameters
        ----------
        index : tuple
            Should be a tuple of integers

        Returns
        -------
        ndarray
        """
        return self.affine.apply(coordinates)

class Warp(Diffeomorphism):
    """
    A diffeomorphism ψ mapping from a population space :math:`M` to a
    subject space $R$. Here, it is

    .. math::

        ψ (x) = δ(x).

    for a warp field :math:`δ`.

    A point [i,j,k] in the index space of the template of :math:`M` maps
    to the coordinates

    .. math::

        ψ(reference[i,j,k]) = warp[i,j,k]

    An example: The FSL FNIRT warp coefficients file is saved as a Warp
    instance.

    Parameters
    ----------
    reference : Affine or ndarray, shape (4,4)
        An affine transformation
    warp : ndarray (...,3)
        A warp field
    vb : str or Identifier
        An object that will be used to identify the domain (vb) of this
        diffeomorphism.
    nb : str or Identifier
        An object that will be used to identify the image (nb) of this
        diffeomorphism.
    metadata : dict
        If meta data is provided which give information on which data
        had been used for the fit: which had provided, should be a dict
        with at least the following fields: ``{'vb_file':
        path/which/defined/the/domain, 'nb_file':
        path/which/defined/the/image,}``.
    """
    def __init__(self, reference, warp, vb=None, nb=None, name=None, metadata=None):
        assert type(warp) is np.ndarray, 'warp must be numpy.ndarray'
        assert 3 == warp.shape[-1], 'last dimension of warp must be 3'

        if type(reference) is Affine:
            self.reference = reference
        else:
            self.reference = Affine(reference)

        self.warp = warp
        self.shape = self.warp.shape[:-1]
        self.vb = vb
        self.nb = nb
        self.name = name

        if metadata is not None:
            self.metadata = metadata

    def apply_to_index(self, index):
        return self.warp[index]

    def apply_to_indices(self, indices):
        return self.warp[indices]

    def apply(self, coordinates):
        """
        Apply diffeomorphism to the point at given coordinate

        Parameters
        ----------
        index : tuple
            Should be a tuple of integers

        Returns
        -------
        ndarray

        Notes
        -----
        There are interpolations at place here, which makes this
        potentially slow for large queries.
        """
        indices = self.reference.inv().apply(coordinates)
        newx = map_coordinates(self.warp[...,0], indices.T.reshape(3,-1))
        newy = map_coordinates(self.warp[...,1], indices.T.reshape(3,-1))
        newz = map_coordinates(self.warp[...,2], indices.T.reshape(3,-1))
        return np.vstack((newx, newy, newz)).T

class Displacement(Diffeomorphism):
    """
    A diffeomorphism ψ mapping from a population space :math:`M` to a
    subject space $R$. Here, it is

    .. math::

        ψ (x) = x + δ(x).

    for a displacement field :math:`δ`.

    A point [i,j,k] in the index space of the template of :math:`M` maps
    to the coordinates

    .. math::

        ψ(reference[i,j,k]) = reference[i,j,k] + displacement[i,j,k]

    Parameters
    ----------
    reference : Affine or ndarray, shape (4,4)
        An affine transformation.
    displacement : ndarray (?,?,?,3)
        A displacement field
    vb : str or Identifier
        An object that will be used to identify the domain (vb) of this
        diffeomorphism.
    nb : str or Identifier
        An object that will be used to identify the image (nb) of this
        diffeomorphism.
    metadata : dict
        If meta data is provided which give information on which data
        had been used for the fit: which had provided, should be a dict
        with at least the following fields: ``{'vb_file':
        path/which/defined/the/domain, 'nb_file':
        path/which/defined/the/image,}``.
    """

    def __init__(self, reference, displacement, vb=None, nb=None, name=None):
        assert type(displacement) is np.ndarray, 'displacement must be numpy.ndarray'
        assert 3 == displacement.shape[-1], 'last dimension of displacement must be 3'

        if type(reference) is Affine:
            self.reference = reference
        else:
            self.reference = Affine(reference)

        self.displacement = displacement
        self.shape = self.displacement.shape[:-1]
        self.vb = vb
        self.nb = nb
        self.name = name

        if metadata is not None:
            self.metadata = metadata

    def apply_to_index(self, index):
        return self.reference.apply_to_index(index) + self.displacement[index]

    def apply_to_indices(self, indices):
        return self.reference.apply_to_indices(indices) + self.displacement[indices]

    # TODO: test this function
    def apply(self, coordinates):
        """
        Apply diffeomorphism to the point at given coordinate

        Parameters
        ----------
        index : tuple
            Should be a tuple of integers

        Returns
        -------
        ndarray

        Notes
        -----
        There are interpolations at place here, which makes this
        potentially slow for large queries.
        """
        indices = self.reference.inv().apply(coordinates)
        newx = map_coordinates(self.displacement[...,0], indices.T.reshape(3,-1))
        newy = map_coordinates(self.displacement[...,1], indices.T.reshape(3,-1))
        newz = map_coordinates(self.displacement[...,2], indices.T.reshape(3,-1))
        return coordinates + np.vstack((newx, newy, newz)).T
