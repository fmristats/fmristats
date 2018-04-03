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
        The diffeomorphism ψ that maps from each point in the data field
        of the template brain to the subject reference space of an FMRI
        session.
    template : None, Image or ndarray
        A template image in population space.
    template_mask : None, Image or ndarray
        A template mask image in population space. By default the
        template mask image will be used in many operations. The idea:
        You provide a detailed image, say also containing skull, in
        template but only the brain as in the template mask image.
    target : None, Image or ndarray
        A target image in population space.

    Notes
    -----
    If template is an Image, it will be checked whether its reference
    agrees with the reference of the diffeomorphism.  If template is an
    ndarray, it must be equal to the shape stored with the
    diffeomorphism, as it will then be assumed that they share the same
    reference: the template's reference will be set identical to the
    reference of the diffeomorphism.

    If the software which has been used to fit the diffeomorphism
    creates a warped image of its input, it may be stored in `target`.
    """

    def __init__(self, diffeomorphism, template=None,
            template_mask=None, target=None, metadata=None):

        assert issubclass(type(diffeomorphism), Diffeomorphism), \
                'diffeomorphism must be a subclass of Diffeomorphism'
        assert type(diffeomorphism.nb) is Identifier, \
                'nb of diffeomorphism must be Identifier'

        self.diffeomorphism = diffeomorphism
        self.name = self.diffeomorphism.nb

        if template:
            self.set_template(template=template)

        if template_mask:
            self.set_template_mask(template_mask=template_mask)

        if target:
            self.set_target(target=target)

        self.metadata = metadata

    def set_template(self, template):
        """
        Set or reset template image

        Parameters
        ----------
        template : Image or ndarray
            The data array.

        Notes
        -----
        If template is an Image, it will be checked whether its
        reference agrees with the reference of the diffeomorphism.  If
        template is an ndarray, it must be equal to the shape stored
        with the diffeomorphism, as it will then be assumed that they
        share the same reference: the template's reference will be set
        identical to the reference of the diffeomorphism.
        """
        assert template.shape == self.diffeomorphism.shape, \
                'shapes of template and diffeomorphism must match'

        if type(template) is Image:
            assert isclose(template.reference, self.diffeomorphism.reference), \
                    'references of template and diffeomorphism must match'
            template.reference = self.diffeomorphism.reference
            self.template = template
        else:
            self.template = Image(
                    reference=self.diffeomorphism.reference,
                    data=template,
                    name=self.diffeomorphism.vb)

    def set_template_mask(self, template_mask):
        """
        Set or reset template mask image

        Parameters
        ----------
        template_mask : Image or ndarray
            The data array.

        Notes
        -----
        If template_mask is an Image, it will be checked whether its
        reference agrees with the reference of the diffeomorphism.  If
        template_mask is an ndarray, it must be equal to the shape
        stored with the diffeomorphism, as it will then be assumed that
        they share the same reference: the template_mask's reference
        will be set identical to the reference of the diffeomorphism.
        """
        assert template_mask.shape == self.diffeomorphism.shape, \
                'shapes of template_mask and diffeomorphism must match'

        if type(template_mask) is Image:
            assert isclose(template_mask.reference, self.diffeomorphism.reference), \
                    'references of template_mask and diffeomorphism must match'
            template_mask.reference = self.diffeomorphism.reference
            self.template_mask = template_mask
        else:
            self.template_mask = Image(
                    reference=self.diffeomorphism.reference,
                    data=template_mask,
                    name=self.diffeomorphism.vb)

    def set_target(self, target):
        """
        Set or reset target image

        Parameters
        ----------
        target : Image or ndarray, shape (x,y,z)
            The data array.

        Notes
        -----
        If target is an Image, it will be checked whether its
        reference agrees with the reference of the diffeomorphism.  If
        target is an ndarray, it must be equal to the shape stored
        with the diffeomorphism, as it will then be assumed that they
        share the same reference: the target's reference will be set
        identical to the reference of the diffeomorphism.
        """
        assert target.shape == self.diffeomorphism.shape, \
                'shapes of target and diffeomorphism must match'

        if type(target) is Image:
            assert isclose(target.reference, self.diffeomorphism.reference), \
                    'references of target and diffeomorphism must match'
            target.reference = self.diffeomorphism.reference
            self.target = target
        else:
            self.target = Image(
                    reference=self.diffeomorphism.reference,
                    data=target,
                    name=self.diffeomorphism.vb)

    def describe(self):
        description = """
            Available images in population space
            ------------------------------------
            Template:      {}
            Template mask: {}
            Target:        {}
        """
        return description.format(
                hasattr(self, 'template'),
                hasattr(self, 'template_mask'),
                hasattr(self, 'target'),
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
    the coordinate system of the FMRI session, if this space is indeed
    *close* to the subject reference space of the FMRI session This may
    be archived by calling the :func:`reset_reference_space` attribute
    of the session first. .

    See also
    --------
    `pmap_reference`
    `pmap_scan`
    """
    diffeomorphism = Identity(
            reference=session.reference,
            shape=session.shape,
            vb='scanner',
            nb=session.name)

    return PopulationMap(diffeomorphism)

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
    the coordinate system of the FMRI session, if this space is indeed
    *close* to the subject reference space of the FMRI session This may
    be archived by calling the :func:`reset_reference_space` attribute
    of the session first. .

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
            vb='reference',
            nb=session.name)

    return PopulationMap(diffeomorphism)

def pmap_scan(session, scan_cycle):
    """
    Pick a scan reference as the population map

    Parameters
    ----------
    session : Session
        Session instance
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
            affine=session.scan_inverses[scan_cycle],
            shape=session.shape,
            vb='reference',
            nb=session.name)

    return PopulationMap(diffeomorphism)
