Rigorous statistical modelling of functional MRI data of the brain
==================================================================

The current state-of-the-art approach to the statistical analysis of
functional MR-images involves a variety of pre-processing steps, which
alter the signal to noise ratio of the original data.

This is a new and original approach for the statistical analysis of
functional MR-imaging data of brain scans. The method essentially fits a
weighted least squares model to arbitrary points of a 3D-random field.
Without prior spacial smoothing, i.e., without altering the original
4D-image, the method nevertheless results in a smooth fit of the
underlying activation pattern. More importantly, though, the method
yields a trustworthy estimate of the uncertainty of the estimated
activation field for each subject in a study. The availability of this
uncertainty field allows for the first time to model group studies and
group-wise comparisons using random effects meta regression models,
acknowledging the fact that (i) individual subjects are random entities
in group studies, and that (ii) the variability in the estimated
individual activation patterns varies across the brain and between
subjects.
