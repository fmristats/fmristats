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

    Naming convetions for images in `vb` and `nb` seem to be
    non-existant, ranging from `input`, `target`, `template`,
    `reference`, `moving`, `fixed`, and more: Although the induitive
    meaning of what these images shall represent is typically clear from
    teh context, the question which of these images live in `vb` or `nb`
    is not clear and ambigious at best. When fitting a diffeomorphism
    providing the algorithm with the arguments `reference` and `input`
    (FNIRT) or `fixed` and `moving` (ANTS): does the diffeomorphism go
    from `reference` to `input` or vice versa?
    """

    def __init__(self, diffeomorphism,
            vb=None,
            nb=None,
            vb_background=None,
            nb_background=None,
            vb_estimate=None,
            nb_estimate=None,
            name=None):

        assert issubclass(type(diffeomorphism), Diffeomorphism), \
                'diffeomorphism must be a subclass of Diffeomorphism'
        assert type(diffeomorphism.nb) is Identifier, \
                'nb of diffeomorphism must be Identifier'

        self.diffeomorphism = diffeomorphism
        self.name = name

        if vb:
            self.set_vb(image=vb)

        if vb_background:
            self.set_vb_background(image=vb_background)

        if vb_estimate:
            self.set_vb_estimate(image=vb_estimate)

        if nb:
            self.set_nb(image=nb)

        if nb_background:
            self.set_nb_background(image=nb_background)

        if nb_estimate:
            self.set_nb_estimate(image=nb_estimate)

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

    def describe(self):
        description = """
        Population map
        --------------
        name:          {:s}
        vb:            {:s}
        nb:            {:s}
        vb_background: {:s}
        nb_background: {:s}
        vb_estimate:   {:s}
        nb_estimate:   {:s}"""

        try:
            vb = self.vb.name
        except:
            vb = '--'

        try:
            nb = self.nb.name.name()
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

        return description.format(
                self.name,
                vb, nb, vb_background, nb_background, vb_estimate, nb_estimate,
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

def pmap_scanner(session):
    """
    The identity map

    Sets the population space equal to the coordinate system of the
    scanner and the population map equal to the identity, and the
    template grid to native resolution.

    Parameters
    ----------
    session : Session
        A FMRI session

    Notes
    -----
    When only studying a single subject, it makes perfect sense to
    set the population space equal to an isometric image of the
    subject's head, i.e. to a population space which preserves distances
    with respect to the subject.

    In other words, for a single subject, the population space can be
    set equal to the subject reference space :math:`R` or to
    :math:`ρ_s[R]` for some scan reference :math:`ρ_s`. This function
    will do the former and sets the resolution of the template equal to
    the resolution of the acquisition grid. (When fitting on this grid,
    the signal model is estimated in *native* resolution).

    It only makes sense, though, to set the population space equal to
    the coordinate system of the FMRI session, if this space is *close*
    to the subject reference space of the FMRI session. This may be
    archived by calling :func:`reset_reference_space` of
    the reference_maps: This will have the effect that the population
    space is then equal to the mean location of the subject in the
    scanner.


    See also
    --------
    `pmap_reference`
    `pmap_scan`
    """
    diffeomorphism = Identity(
            reference=session.reference,
            shape=session.shape,
            vb=session.name.name(),
            nb=session.name,
            name='scanner')

    return PopulationMap(diffeomorphism, name='self')

def pmap_reference(session, resolution=2.):
    """
    The identity map

    Sets the population space equal to the coordinate system of the
    scanner and the population map equal to the identity, and the
    template grid to given resolution.

    Parameters
    ----------
    session : Session
        A FMRI session
    resolution : float
        Resolution in milli meter, default 2.

    Notes
    -----
    When only studying a single subject, it makes perfect sense to
    set the population space equal to an isometric image of the
    subject's head, i.e. to a population space which preserves distances
    with respect to the subject.

    In other words, for a single subject, the population space can be
    set equal to the subject reference space :math:`R` or to
    :math:`ρ_s[R]` for some scan reference :math:`ρ_s`. This function
    will do the former and sets the resolution to the given value.

    It only makes sense, though, to set the population space equal to
    the coordinate system of the FMRI session, if this space is *close*
    to the subject reference space of the FMRI session. This may be
    archived by calling :func:`reset_reference_space` of
    the reference maps. This will have the effect that the population
    space is then equal to the mean location of the subject in the
    scanner.

    See also
    --------
    :func:`pmap_scanner`
    :func:`pmap_scan`
    """

    ref = session.reference
    mat = ref.affine[:3,:3] / (ref.resolution() / resolution)

    reference = from_cartesian(mat, ref.affine[:3,3])
    shape = tuple(np.ceil(reference.inv().apply(ref.apply_to_index(session.shape))).astype(int))

    diffeomorphism = Identity(
            reference=reference,
            shape=shape,
            vb=session.name.name(),
            nb=session.name,
            name='reference')

    return PopulationMap(diffeomorphism, name='self')

def pmap_scan(session, reference_maps, scan_cycle):
    """
    Pick a scan reference as the population map

    Parameters
    ----------
    session : Session
        Session instance
    reference_maps : ReferenceMaps
        Reference Maps
    scan_cycle : int
        Reference scan cycle.

    Notes
    -----
    When only studying a single subject, it makes perfect sense to
    set the population space equal to an isometric image of the
    subject's head, i.e. to a population space which preserves distances
    with respect to the subject.

    In other words, for a single subject, the population space can be
    set equal to the subject reference space :math:`R` or to
    :math:`ρ_s[R]` for some scan reference :math:`ρ_s`. This function
    will do the latter for the given scan, and sets the resolution of
    the template equal to the resolution of the acquisition grid. (When
    fitting on this grid, the signal model is estimated in *native*
    resolution).

    See also
    --------
    :func:`pmap_reference`
    :func:`pmap_scanner`
    """
    diffeomorphism = AffineTransformation(
            reference=session.reference,
            affine=reference_maps.scan_references.inv().affines[scan_cycle],
            shape=session.shape,
            vb=session.name.name(),
            nb=session.name,
            name='cycle{:d}'.format(scan_cycle))

    return PopulationMap(diffeomorphism, name='self')
