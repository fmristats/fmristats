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

########################################################################
# Miscellaneous
########################################################################

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
