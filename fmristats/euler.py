"""
Euler angels
------------

Code for calculating Euler angels from rotation matrices
"""

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

# The following code has been almost identically taken from:
#   http://www.learnopencv.com/rotation-matrix-to-euler-angles/

import numpy as np
import math

def is_rotation_matrix(a):
    """
    Checks if a matrix is a valid rotation matrix.
    """
    shouldBeIdentity = np.dot(np.transpose(a), a)
    I = np.identity(3, dtype = a.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6

def rotation_matrix_to_euler_angles(a) :
    """
    Calculates rotation matrix to Euler angles

    The result is the same as MATLAB except the order
    of the Euler angles (x and z are swapped).
    """
    assert(is_rotation_matrix(a))
    sy = math.sqrt(a[0,0] * a[0,0] +  a[1,0] * a[1,0])
    singular = sy < 1e-6

    if  not singular :
        x = math.atan2(a[2,1] , a[2,2])
        y = math.atan2(-a[2,0], sy)
        z = math.atan2(a[1,0], a[0,0])
    else :
        x = math.atan2(-a[1,2], a[1,1])
        y = math.atan2(-a[2,0], sy)
        z = 0

    return np.array([x, y, z])
