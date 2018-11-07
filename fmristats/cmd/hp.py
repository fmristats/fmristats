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

mat = """a Matlab coded stimulus design of on- and offsets."""

########################################################################
# General arguments
########################################################################

session = """path to a session file or a template which defines the path
where to find or save the respected session files when using a protocol
file."""

sfit = """path to a result file or a template which defines the path
where to find or save the respected result files when using a protocol
file."""

vb = """name of a population space."""

vb_nii = """name of a population space."""

vb_name = """name of a population space."""

vb_path = """name of a population space."""

vb_background = """name of a population space."""

vb_background_nii = """name of a population space."""

vb_background_name = """name of a population space."""

vb_background_path = """name of a population space."""

nb = """name of a population space."""

nb_nii = """name of a population space."""

nb_name = """name of a population space."""

nb_path = """name of a population space."""

population_space_directory = """directory in which to save the
population space."""

population_space = vb

population_map = """instance of a population map."""

stimulus = """stimulus file.  File name containing the stimulus
instance for this session.  An stimulus instance contains all
information of the paradigm of this session, i.e., whether the paradigm
follows a block design, the names of the respective blocks, their
respective onsets and respective durations."""

template = """template file used for the population space. The file
should contain a 3D-image of a brain in any file format understood by
the NiBabel project, e.g, any of ANALYZE (plain, SPM99, SPM2 and later),
GIFTI, NIfTI1, NIfTI2, MINC1, MINC2, MGH and ECAT as well as Philips
PAR/REC (For more details see http://nipy.org/nibabel/). On the other
hand, it should also be understood by FSL."""

template_mask = """template mask file used for the population space. The
file should contain a 3D-image of a brain in any file format understood
by the NiBabel project, e.g, any of ANALYZE (plain, SPM99, SPM2 and
later), GIFTI, NIfTI1, NIfTI2, MINC1, MINC2, MGH and ECAT as well as
Philips PAR/REC (For more details see http://nipy.org/nibabel/). On the
other hand, it should also be understood by FSL."""

reference_maps = """reference maps."""

########################################################################
# Miscellaneous
########################################################################

verbose = """increase output verbosity."""

quite = """do not be verbose."""

skip = """do not create or fit anything"""

force = """by default {} files will not be overwritten but parsed
instead.  This will force to recreate any files."""

ignore_lock = """ignore whether a file has been locked by a different
program. Only use this if orphans exists, and you now for sure no other
program is running that may have looked the respected files. Together
with --skip this may be used to delete orphan lock files. Short -is"""

remove_lock = """removes the lock."""

# are these old?

delete_invalid = """delete all {} files which entries in the protocol
which are marked as invalid."""

conditional_force = """by default {} files will not be overwritten but
parsed instead.  This will force to re-create any of files on which
creation this file depends has a newer time stamp."""

fix = """by default {} files will not be overwritten but parsed
instead.  This will force to re-create any of these files."""

remove = """by default {} files will not be overwritten but parsed
instead.  This will force to re-create any of these files."""


########################################################################
# Arguments specific for using the protocol API
########################################################################

protocol = """
protocol file.
"""

protocol_log = """
Log file. Entries that could not be processed are marked as invalid
"""

cohort = """only use entries instances which belong to this cohort. If
no protocol file is used, then this is the name of the cohort to which
this subject belongs."""

j = """only use entries which belong to this subject id. If no protocol
file is used, then this is the id of the subject in its cohort."""

datetime = """only use entries which are stamped with this date and
time."""

paradigm = """only use entries which belong to this stimulus design.
The name of the stimulus will appear as part of the default file names
produced by various fmristats command line tools.  Unless you are
planing to provide custom file names to all of fmristats' command line
tools, it is recommend not to use any special characters in this
name."""

strftime = """convert time to string according to this format specification."""

########################################################################
# Arguments specific for the setup of a session instances
########################################################################

nii = """
the output of a fMRI session.  The file should contain a 4D-image of a
fMRI session in any file format understood by the NiBabel project, e.g,
any of ANALYZE (plain, SPM99, SPM2 and later), GIFTI, NIfTI1, NIfTI2,
MINC1, MINC2, MGH and ECAT as well as Philips PAR/REC.  For more details
see http://nipy.org/nibabel/.  Please note that fmristats has only been
tested with Nifti1 files.
"""

epi_code = """
code for the direction of the normal vector of the EPIs assuming the
data is in RAS+. For example, if the EPIs have been measured inferior to
superior, then set EP to 2, if they have been measured superior to
inferior, set EP to -2, if they have been measured posterior to
anterior, i.e., they are parallel to the left--right-inferior--superior
plain, then set EP to 1.  Values in a protocol file will take precedence
"""

onsetx = """times of onsets of condition x, which will be interpreted
as the control stimulus of the subject.  Should, e.g., be given in the
form:

                --onsetx 2 84 166 248 330

to be parsed correctly.  If onsets are not given in order, they will be
sorted automatically.  Must be specified if --stimulus is not.
"""

onsety = """times of onset of condition y, which will be interpreted as
the stimulus stimulus of the subject.  Should, e.g., be given in the
form:

                --onsety 37 119 201 283

to be parsed correctly.  If onsets are not given in order, they will be
sorted automatically.  Must be specified if --stimulus is not.  """

durationx = """duration of the stimulus of the subject under condition
x.  Must be specified if --stimulus is not."""

durationy = """duration of the stimulus of the subject under condition
y.  Must be specified if --stimulus is not.  """

detect_foreground = """detect foreground."""

set_foreground = """set foreground."""

grubbs = """an outlier detection is performed to identify scans which
may have been acquired during severe head movements. More precisely, a
Grubbs' outlying test will be performed on the set of estimated
principle semi axis for each full scan cycle on the given level of
significance.

When using fmririgid to create ReferenceMaps, the default is 0.1, and
the information of outlying scans is saved to disk together with the
estimated rigid body transformations. Then, when running fmrifit,
this information is used. When setting --grubbs in fmrifit, outlier
estimation is performed again."""

no_tracking = """skip tracking of the rigid body, i.e., skip fitting
head movements to the data.  This will have the effect that head
movements are set to the identity on the reference space.  This is
useful if your data consists of phantom data.  """

########################################################################
# Arguments specific for the RSM Signal Model: design matrix
########################################################################

burn_in = """acquisition burn in. The *first* scan cycle that should
enter the model."""

time_order = """degree of polynomial that models temporal fluctuation.
Default is 1."""

drop_stimulus = """if stimulus information is present in the model,
here, you may drop it."""

########################################################################
# Arguments specific for the RSM Signal Model
########################################################################

scale = """
standard deviation of a Gaussian kernel that defines the weighting
scheme of the underlying WLS regression.  If not given explicitly, scale
will be set to one half of the length of the diagonal of the
orthorhombic measure lattice of the session (and this is recommended).
This default behaviour corresponds to the SCALE_TYPE ``diagonal`` and
can be overwritten by setting SCALE_TYPE to a different value. The
parameter SCALE will determine the final curvature of the fitted effect
field. The larger SCALE, the flatter the fitted effect field will
appear.
"""

scale_type = """
should be either one of ``diagonal``, ``max``, or ``min``.

If SCALE_TYPE is ``diagonal``, then scale is equal to one half of the
length of the diagonal of the orthorhombic measure lattice.  If
SCALE_TYPE is ``min`` (or ``max``), then SCALE is one half of the
minimal (or maximal) edge length of the orthorhombic measure lattice.
"""

stimulus_block = """
name of the stimulus block.
"""

control_block = """
name of the control block.
"""

population_mask = """Fit the effect field **only** at points which are
within the template brain of the population space.  You should only use
this option, if you are sure that your diffeomorphism from population to
reference space of this sample is dead on. If this is the case -- and it
should be, right? -- this will considerably speed up the computation
time.  It does, however, make visual control of the fit more difficult,
and should be used with care."""

factor = """
Only observations within an area of FACTOR standard deviations under
the Gaussian distribution defined by SCALE will find their way into the
WLS regression.
"""

mass = """
Only observations within an area of mass `mass` under the Gaussian
distribution defined by SCALE will find their way into the WLS
regression.  If specified, must be a value between 0 and 1, and should
be a value close to 1.
"""

offset_beginning = """
The haemodynamic response to stimulus is not immediate.  It is usually
assumed that the HR-function spikes approximately five seconds after the
first stimulus.  The value of OFFSET_BEGINNING is added to the onset
times of the stimulus phases to allow you to wait till the subject's
brain is in the respected stimulus modus.
"""

offset_end = """
Similar to OFFSET_BEGINNING the value OFFSET_END is removed from the end
of each stimulus phase and not considered in the fitting.
"""


########################################################################
# Arguments specific for the RSM Population Model
########################################################################

########################################################################
# Epilog
########################################################################

epilog = """
For more information visit https://fmristats.github.io/
"""

cores = """
Number of cores to use. Default is the number of cores on the machine
"""
