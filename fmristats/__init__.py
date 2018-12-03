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

from .name import Identifier

from .affines import Affine, Affines

from .diffeomorphisms import Image, Diffeomorphism

from .stimulus import Stimulus, Block

from .session import Session

from .reference import ReferenceMaps

from .pmap import PopulationMap

from .smodel import SignalModel, Result

from .sample import Sample

from .pmodel import PopulationModel, PopulationResult

from .tau import tau

from .load import load
