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

A population map is a diffeomorphism from a population space to the
subject reference space of a FMRI session.

"""

from .name import Identifier

from .affines import Affine, isclose, from_cartesian

from .diffeomorphisms import Diffeomorphism, Image, Identity, AffineTransformation

import pickle

import numpy as np

#######################################################################
#
# PopulationMap
#
#######################################################################

class PopulationMap:
    """
    A population map is a diffeomorphism from a population space to the
    reference space of a FMRI session.

    Parameters
    ----------
    diffeomorphism : subclass of Diffeomorphism
        The diffeomorphism ψ that maps from vb (the population space) to
        nb (the subject space).
    vb : None, Image or ndarray
        An image in population space.
    nb : None, Image or ndarray
        An image in subject space.
    vb_background : None, Image or ndarray
        An image in population space.
    nb_background : None, Image or ndarray
        An image in subject space.
    vb_estimate : None, Image or ndarray
        An image in population space.
    nb_estimate : None, Image or ndarray
        An image in subject space.
    vb_mask : None, Image or array
        An image in population space.
    nb_mask : None, Image or array
        An image in subject space.
    vb_ati : None, Image or array
        An image in population space that will serve as the reference
        field for the BOLD unit.

    Notes
    -----
    If vb (or any of the other arguments) is an Image, it will be
    checked whether its *reference affine* agrees with the reference
    affine of the diffeomorphism ψ.  If vb is an ndarray, its shape must
    be equal to the shape of the diffeomorphism, as it will then be
    assumed that they share the same reference affine: the vb's
    reference will be set identical to the reference of the
    diffeomorphism.

    If the software which has been used to fit the diffeomorphism
    creates a warped images of its inputs, they can be stored in
    `vb_estimate` and `nb_estimate` respectively.

    There appears to be no clear naming convention for images in `vb`
    and `nb`. Examples for images in `vb` and `nb` are `input`,
    `target`, `template`, `reference`, `moving`, `fixed`, and more: The
    intuitive meaning of what these images shall represent is often (and
    not always) clear from the context, but the question which of these
    images live in `vb` or `nb` is ambiguous at best. Say, when fitting
    a diffeomorphism providing an algorithm with the arguments
    `reference` and `input` (FNIRT) or `fixed` and `moving` (ANTS): does
    the diffeomorphism go from `reference` to `input` or vice versa?
    """

    def __init__(self, diffeomorphism,
            vb=None,
            nb=None,
            vb_background=None,
            nb_background=None,
            vb_estimate=None,
            nb_estimate=None,
            vb_mask=None,
            nb_mask=None):

        assert issubclass(type(diffeomorphism), Diffeomorphism), \
                'diffeomorphism must be a subclass of Diffeomorphism'
        assert type(diffeomorphism.nb) is Identifier, \
                'nb of diffeomorphism must be Identifier'

        self.diffeomorphism = diffeomorphism
        self.name = self.diffeomorphism.nb

        if vb:
            self.set_vb(image=vb)

        if vb_background:
            self.set_vb_background(image=vb_background)

        if vb_estimate:
            self.set_vb_estimate(image=vb_estimate)

        if vb_mask:
            self.set_vb_mask(image=vb_mask)

        if nb:
            self.set_nb(image=nb)

        if nb_background:
            self.set_nb_background(image=nb_background)

        if nb_estimate:
            self.set_nb_estimate(image=nb_estimate)

        if nb_mask:
            self.set_nb_mask(image=nb_mask)

    def set_vb(self, image):
        """
        Set or reset the image in vb

        Parameters
        ----------
        image : Image or ndarray
            The image or data array.

        Notes
        -----
        If image is an Image, it will be checked whether its reference
        agrees with the reference of the diffeomorphism.  If image is an
        ndarray, it must be equal to the shape stored with the
        diffeomorphism, as it will then be assumed that they share the
        same reference: the image's reference will be set identical to
        the reference of the diffeomorphism.
        """
        if type(image) is Image:
            assert image.shape == self.diffeomorphism.shape, \
                    'shapes of image and diffeomorphism must match'
            assert isclose(image.reference, self.diffeomorphism.reference), \
                    'references of image and diffeomorphism must match'
            image.reference = self.diffeomorphism.reference
            self.vb = image
        else:
            self.vb = Image(
                    reference=self.diffeomorphism.reference,
                    data=image,
                    name=self.diffeomorphism.vb)

    def set_vb_background(self, image):
        """
        Set or reset the background image in vb

        Parameters
        ----------
        image : Image or ndarray
            The image or data array.

        Notes
        -----
        If image is an Image, it will be checked whether its reference
        agrees with the reference of the diffeomorphism.  If image is an
        ndarray, it must be equal to the shape stored with the
        diffeomorphism, as it will then be assumed that they share the
        same reference: the image's reference will be set identical to
        the reference of the diffeomorphism.
        """
        if type(image) is Image:
            assert image.shape == self.diffeomorphism.shape, \
                    'shapes of image and diffeomorphism must match'
            assert isclose(image.reference, self.diffeomorphism.reference), \
                    'references of image and diffeomorphism must match'
            image.reference = self.diffeomorphism.reference
            self.vb_background = image
        else:
            self.vb_background = Image(
                    reference=self.diffeomorphism.reference,
                    data=image,
                    name=self.diffeomorphism.vb)

    def set_vb_estimate(self, image):
        """
        Set or reset the estimate for the image in vb

        Parameters
        ----------
        image : Image or ndarray
            The image or data array.

        Notes
        -----
        If image is an Image, it will be checked whether its reference
        agrees with the reference of the diffeomorphism.  If image is an
        ndarray, it must be equal to the shape stored with the
        diffeomorphism, as it will then be assumed that they share the
        same reference: the image's reference will be set identical to
        the reference of the diffeomorphism.
        """
        if type(image) is Image:
            assert image.shape == self.diffeomorphism.shape, \
                    'shapes of image and diffeomorphism must match'
            assert isclose(image.reference, self.diffeomorphism.reference), \
                    'references of image and diffeomorphism must match'
            image.reference = self.diffeomorphism.reference
            self.vb_estimate = image
        else:
            self.vb_estimate = Image(
                    reference=self.diffeomorphism.reference,
                    data=image,
                    name=self.diffeomorphism.vb)

    def set_vb_mask(self, image):
        """
        Set or reset the estimate for the image in vb

        Parameters
        ----------
        image : Image or ndarray
            The image or data array.

        Notes
        -----
        If image is an Image, it will be checked whether its reference
        agrees with the reference of the diffeomorphism.  If image is an
        ndarray, it must be equal to the shape stored with the
        diffeomorphism, as it will then be assumed that they share the
        same reference: the image's reference will be set identical to
        the reference of the diffeomorphism.
        """
        if type(image) is Image:
            assert image.shape == self.diffeomorphism.shape, \
                    'shapes of image and diffeomorphism must match'
            assert isclose(image.reference, self.diffeomorphism.reference), \
                    'references of image and diffeomorphism must match'
            image.reference = self.diffeomorphism.reference
            self.vb_mask = image
        else:
            self.vb_mask = Image(
                    reference=self.diffeomorphism.reference,
                    data=image,
                    name=self.diffeomorphism.vb)

    def set_vb_ati(self, image):
        """
        Set or reset the reference field for the BOLD unit.

        Parameters
        ----------
        image : Image or ndarray
            The image or data array.

        Notes
        -----
        If image is an Image, it will be checked whether its reference
        agrees with the reference of the diffeomorphism.  If image is an
        ndarray, it must be equal to the shape stored with the
        diffeomorphism, as it will then be assumed that they share the
        same reference: the image's reference will be set identical to
        the reference of the diffeomorphism.
        """
        if type(image) is Image:
            assert image.shape == self.diffeomorphism.shape, \
                    'shapes of image and diffeomorphism must match'
            assert isclose(image.reference, self.diffeomorphism.reference), \
                    'references of image and diffeomorphism must match'
            image.reference = self.diffeomorphism.reference
            self.vb_ati = image
        else:
            self.vb_ati = Image(
                    reference=self.diffeomorphism.reference,
                    data=image,
                    name=self.diffeomorphism.vb)

    def set_nb(self, image):
        """
        Set or reset the image in nb

        Parameters
        ----------
        image : Image or ndarray
            The image or data array.
        """
        assert type(image) is Image, 'image must by an Image'
        self.nb = image

    def set_nb_background(self, image):
        """
        Set or reset the background image in nb

        Parameters
        ----------
        image : Image or ndarray
            The image or data array.
        """
        assert type(image) is Image, 'image must by an Image'
        self.nb_background = image

    def set_nb_estimate(self, image):
        """
        Set or reset the estimate of the image in nb

        Parameters
        ----------
        image : Image or ndarray
            The image or data array.
        """
        assert type(image) is Image, 'image must by an Image'
        self.nb_estimate = image

    def set_nb_mask(self, image):
        """
        Set or reset the estimate of the image in nb

        Parameters
        ----------
        image : Image or ndarray
            The image or data array.
        """
        assert type(image) is Image, 'image must by an Image'
        self.nb_mask = image

    def describe(self):
        description = """
        Population map
        --------------
        diffeomorphism: {:s}
        vb:             {:s}
        nb:             {:s}
        vb_background:  {:s}
        nb_background:  {:s}
        vb_estimate:    {:s}
        nb_estimate:    {:s}
        vb_mask:        {:s}
        nb_mask:        {:s}
        vb_ati:         {:s}"""

        try:
            vb = self.vb.name
        except:
            vb = '--'

        try:
            nb = self.nb.name.name()
        except:
            try:
                nb = self.nb.name
            except:
                nb = '--'

        try:
            vb_background = self.vb_background.name
        except:
            vb_background = '--'

        try:
            nb_background = self.nb_background.name
        except:
            nb_background = '--'

        try:
            vb_estimate = self.vb_estimate.name
        except:
            vb_estimate = '--'

        try:
            nb_estimate = self.nb_estimate.name
        except:
            nb_estimate = '--'

        try:
            vb_mask = self.vb_mask.name
        except:
            vb_mask = '--'

        try:
            nb_mask = self.nb_mask.name
        except:
            nb_mask = '--'

        try:
            vb_ati = self.vb_ati.name
        except:
            vb_ati = '--'

        return description.format(
                self.diffeomorphism.name,
                vb, nb, vb_background, nb_background,
                vb_estimate, nb_estimate,
                vb_mask, nb_mask, vb_ati,
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

#######################################################################
#
# Functions for creating canonical population maps
#
#######################################################################

def pmap_scanner(session, reference_maps=None, scan_cycle=None,
        resolution=2., name='scanner'):
    """
    The identity map

    Sets the standard space equal to the coordinate system of the
    scanner and the diffeomorphism from standard space to subject
    reference space to different rigid body transformation.

    Parameters
    ----------
    session : Session
        A FMRI session
    reference_maps : None or ReferenceMaps
        Reference Maps
    scan_cycle : None or int
        Reference scan cycle.
    resolution : str or float
        Resolution (in units of the session reference, this will likely
        be milli meter), default 2. If resolution is set to 'native' the
        resolution of the scanner is used.
    name = str
        A name for this population map.

    Notes
    -----
    When only studying a single subject, it is possible to set the
    standard space to an isometric image of the subject's head, i.e. to
    a space which preserves distances with respect to the subject.

    In other words, for a single subject, the standard space can be
    set equal to the subject reference space :math:`R` or to
    :math:`ρ_t[R]` for any time point :math:`t` in the Session.

    You may even set :math:`R` identical to :math:`S` if scanner and
    subject reference space are reasonable »close«. This may be archived
    by calling :func:`reset_reference_space` of the reference maps. If
    called without arguments, this will have the effect that the
    population space is then equal to the mean location of the subject
    in the scanner.
    """
    if (resolution is None) or (resolution == 'native') or np.isclose(resolution, 0):
        reference = session.reference
        shape = session.shape
    else:
        ref = session.reference
        mat = ref.affine[:3,:3] / (ref.resolution() / resolution)
        reference = from_cartesian(mat, ref.affine[:3,3])
        shape = tuple(np.ceil(reference.inv().apply(ref.apply_to_index(session.shape))).astype(int))

    if reference_maps is None:
        diffeomorphism = Identity(
                reference=reference,
                shape=shape,
                vb=session.name.name(),
                nb=session.name,
                name=name)
    else:
        if scan_cycle is None:
            diffeomorphism = AffineTransformation(
                    reference=reference,
                    affine=reference_maps.acquisition_maps.mean_rigid().inv(),
                    shape=session.shape,
                    vb=session.name.name(),
                    nb=session.name,
                    name=name)
        else:
            diffeomorphism = AffineTransformation(
                    reference=reference,
                    affine=Affine(reference_maps.acquisition_maps.affines[scan_cycle]).inv(),
                    shape=session.shape,
                    vb=session.name.name(),
                    nb=session.name,
                    name=name)

    return PopulationMap(diffeomorphism)
