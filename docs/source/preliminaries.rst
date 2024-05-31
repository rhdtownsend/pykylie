*************
Preliminaries
*************

BRUCE / KYLIE
=============

:mad-star:`BRUCE / KYLIE <bruce>` are a pair of open-source Fortran 77
codes for synthesizing spectra of pulsating stars. BRUCE constructs a
point-sampled model for the surface of a rotating, gravity-darkened
star, and then subjects this model to perturbations arising from one
or more non-radial pulsation modes. Departures from adiabaticity can
be taken into account, as can the Coriolis force through adoption of
the traditional approximation of rotation.

BRUCE writes out a time-series of perturbed surface models. This
series is read in by KYLIE, which synthesizes disk-integrated spectra
for the models by co-adding the specific intensity emitted by each
visible surface point. The specific intensities are evaluated via
interpolation in a pre-calculated grid. More-detailed descriptions of
both BRUCE and KYLIE are provided in :ads_citet:`townsend:1997a` and
:ads_citet:`townsend:1997b`.

PyKYLIE
=======

Historically, it was not possible to distribute a specific intensity
grid with BRUCE/KYLIE. However, this situation changed with the
development of `Multidimensional Specgral Grids (MSG) <msg_>`__, an
open-source package for synthesizing stellar spectra and photometric
colors via interpolation in pre-calculated grids. PyKYLIE is a Python
re-implementation of the original KYLIE code, built on top of MSG.

Obtaining PyKYLIE
=================

The source code for PyKYLIE is hosted in the
:git:`rhdtownsend/pykylie` git repository on :git:`GitHub
<>`. Instructions for installing the software can be found in the
:ref:`install` chapter.

Development Team
================

PyKYLIE remains under active development by the following team:

* `Rich Townsend <http://www.astro.wisc.edu/~townsend>`__ (University of Wisconsin-Madison); project leader

If you are interested in contributing toward further development of
PyKYLIE, please see the :ref:`contributing` chapter.

Acknowledgments
================

MSG has been developed with financial support from the following grants:

* NASA award 80NSSC20K0515.

