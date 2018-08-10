Data oriented method for the fitting of FMRI models
===================================================

Most current approaches to the statistical analysis of functional
magnetic resonance imaging (FMRI) data involve varieties of
preprocessing steps which alter the signal to noise ratio of the
original data.

Enhancing the SNR prior to a formal analysis, though, shakes at primary
principles of statistical decision making and it will generally inflate
the type I error of the analysis.

This is the first statistical software tool which implements the *data
oriented method (DOM)* estimator for FMRI data models, a new and
original method for the statistical analysis of FMRI data of brain
scans. The method fits a weighted least squares model to points of a
random vector field. Without prior spacial smoothings, i.e. without
altering the original 4D-image, the method nevertheless results in
smooth fits of the underlying activation parameter fields. More
importantly, though, the method yields a trustworthy estimate of the
uncertainty of the estimated activation field for each subject in a
study. The availability of these uncertainty fields allows to model FMRI
studies by random effects meta regression models acknowledging that
individual subjects are random entities, and that the variability in the
estimated individual activation patterns vary across the brain and
between subjects.
