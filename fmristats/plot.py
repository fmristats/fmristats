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

"""
Plot effect fields
------------------

Some handy functions for plotting fMRI images and effect fields.
"""

import numpy as np

import matplotlib
import matplotlib.pyplot as pt
import matplotlib.cm as cm

def shift_colour_map(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Shift a colour map

    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Parameters
    ----------
      cmap : cmap
          The matplotlib colormap to be altered
      start : float
          Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : float
          The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the centre of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : float
          Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.

    Returns
    -------
    cmap : The shifted colour map.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    pt.register_cmap(cmap=newcmap)

    return newcmap

def cmap_shifted (dat, orig_cmap=cm.bwr):
    vmaxes = np.array((np.nanmin(dat), np.nanmax(dat)))
    vmax = max(vmaxes)
    vmin = min(vmaxes)
    midp = 1 - vmax/(vmax + abs(vmin))
    cmap_new = shift_colour_map(orig_cmap, midpoint=midp, name='shifted')
    return(cmap_new)

def picture(self, epi_code, nx, ny, slices=None, cmap=None,
        interpolation='none', background=None, mark_peak=False):

    if abs(epi_code) == 3:
        origin = 'lower'
        xlabel = 'Left -> Right'
        ylabel = 'Posterior -> Anterior'
        aspect = self.reference.aspect(1,0)
    if abs(epi_code) == 2:
        origin = 'lower'
        xlabel = 'Left -> Right'
        ylabel = 'Inferior -> Superior'
        aspect = self.reference.aspect(2,0)
    if abs(epi_code) == 1:
        origin = 'lower'
        xlabel = 'Posterior -> Anterior'
        ylabel = 'Inferior -> Superior'
        aspect = self.reference.aspect(2,1)

    imatrix, ixticks, iyticks, slices = self.flatten(epi_code,nx,ny,slices=slices)

    if (imatrix[np.isfinite(imatrix)] < 0).any():
        if cmap == None:
            cmap = cm.bwr
        cmap = cmap_shifted(imatrix, orig_cmap=cmap)
    else:
        if cmap == None:
            cmap = cm.magma

    if background:
        bmatrix, _, _, _ = background.flatten(epi_code,nx,ny,slices=slices)
        pt.imshow(bmatrix,interpolation=interpolation,cmap='Greys',origin=origin,aspect=aspect)

    im = pt.imshow(imatrix,interpolation=interpolation,cmap=cmap,origin=origin,aspect=aspect)
    #pt.colorbar()
    pt.box(on='off')
    pt.tick_params('off')
    pt.xlabel(xlabel)
    pt.ylabel(ylabel)

    ax = pt.axes()
    pt.setp(ax.get_yticklines(),visible=False)
    pt.setp(ax.get_xticklines(),visible=False)
    pt.setp(ax.get_xticklabels(),visible=False)
    pt.setp(ax.get_yticklabels(),visible=False)
    #ax.set_xticks(ixticks, minor=False)
    #ax.set_yticks(iyticks, minor=False)

    if mark_peak:
        peak  = np.nanargmax(imatrix)
        index = np.unravel_index(peak, imatrix.shape)
        pt.axvline(index[1], c='k', linewidth=.8)
        pt.axhline(index[0], c='k', linewidth=.8)

    return im, slices, imatrix

def contour(self, epi_code, nx, ny, slices=None, cmap=None,
        interpolation='none', background=None, mark_peak=False, **kwargs):

    if abs(epi_code) == 3:
        origin = 'lower'
        xlabel = 'Left -> Right'
        ylabel = 'Posterior -> Anterior'
        aspect = self.reference.aspect(1,0)
    if abs(epi_code) == 2:
        origin = 'lower'
        xlabel = 'Left -> Right'
        ylabel = 'Inferior -> Superior'
        aspect = self.reference.aspect(2,0)
    if abs(epi_code) == 1:
        origin = 'lower'
        xlabel = 'Posterior -> Anterior'
        ylabel = 'Inferior -> Superior'
        aspect = self.reference.aspect(2,1)

    imatrix, ixticks, iyticks, slices = self.flatten(epi_code,nx,ny,slices=slices)

    if (imatrix[np.isfinite(imatrix)] < 0).any():
        if cmap == None:
            cmap = cm.bwr
        cmap = cmap_shifted(imatrix, orig_cmap=cmap)
    else:
        if cmap == None:
            cmap = cm.magma

    if background:
        bmatrix, _, _, _ = background.flatten(epi_code,nx,ny,slices=slices)
        pt.imshow(bmatrix,interpolation=interpolation,cmap='Greys',origin=origin,aspect=aspect)

    cs=pt.contour(imatrix, **kwargs)
    pt.box(on='off')
    pt.tick_params('off')
    pt.xlabel(xlabel)
    pt.ylabel(ylabel)

    ax = pt.axes()
    pt.setp(ax.get_yticklines(),visible=False)
    pt.setp(ax.get_xticklines(),visible=False)
    pt.setp(ax.get_xticklabels(),visible=False)
    pt.setp(ax.get_yticklabels(),visible=False)
    #ax.set_xticks(ixticks, minor=False)
    #ax.set_yticks(iyticks, minor=False)
    #pt.colorbar()

    if mark_peak:
        peak  = np.nanargmax(imatrix)
        index = np.unravel_index(peak, imatrix.shape)
        pt.axvline(index[1], c='k', linewidth=.8)
        pt.axhline(index[0], c='k', linewidth=.8)

    return cs, slices, imatrix
