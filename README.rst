Modelling the data and not the images in FMRI
=============================================

Current approaches to the analysis of functional magnetic resonance
imaging (FMRI) data apply various preprocessing steps to the original
FMRI. These preprocessings lead to a general underestimation of residual
variance in the downstream analysis. This negatively impacts the type I
error of statistical tests and increases the risk for reporting false
positive results.

This is the first statistical software tool which implements the *model
based* (MB) estimator for FMRI data models. It is a new and original
method for the statistical analysis of FMRI of brain scans. MB
estimation combines all preprocessing steps of the standard approaches
into one single modelling step. Without altering the original 4D-image,
the method results in smooth fits of the underlying parameter fields.
More importantly, the method yields a trustworthy estimate of the
uncertainty in BOLD effect estimation.

The availability of these uncertainty fields allows to model FMRI
studies by random effects meta regression models, acknowledging that
individual subjects are random entities, and that the certainty at which
the actual BOLD effect in an individual can be estimated from an FMRI
varies across the brain and between the subjects.

MB estimation allows to process and report BOLD effects in ati units. In
particular multicentre studies gain power by its use: if an effect is
present in your data, you will be more likely to find it.

Citing the MB estimator and this software:

    Thomas W. D. Möbius (2018) Modelling the data and not the images in
    FMRI, ArXiv e-prints, arXiv:1809.07232

    Thomas W. D. Möbius (2018) fmristats: Modelling the data and not the
    images in FMRI (Version 0.1.0) [Computer program]. Available at
    http://fmristats.github.io/

Thank you for citing this project.
