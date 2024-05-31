.. _install:

************
Installation
************

This chapter reviews a couple of different ways of installing PyKYLIE
on your system.

Pre-requisites
==============

PyKYLIE requires that the following Python packages are already present:

* `NumPy <https://numpy.org/>`__
* `SciPy <https://scipy.org/>`__
* `Astropy <https://www.astropy.org/>`__
* `PyMSG <https://msg.readthedocs.io/>`__

If you opt to install from PyPI (below), then these pre-requisites
should be taken care of automatically, with the exception of PyMSG
(which has to be installed manually).

Installing from PyPI
====================

To install PyKYLIE from the `Python Package Index (PyPI)
<https://pypi.org/>`__, use the :command:`pip` command:

.. prompt:: bash

   pip install pykylie

If PyKYLIE is already installed, you can upgrade to a more-recent
version via

.. prompt:: bash

   pip install --upgrade pykylie

Installing from Source
======================

To install PyKYLIE from source, download the `source code
<tarball_>`__ and unpack it from the command line using the
:command:`tar` utility:

.. prompt:: bash
   :substitutions:

   tar xf pykylie-|version|.tar.gz

Then, use :command:`pip` to install it:

.. prompt:: bash
   :substitutions:

   pip install ./pykylie-|version|


	 



