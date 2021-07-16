.. _install:

===============
1. Installation
===============

1.1 Download SPINNER
====================
Please download the file from https://github.com/MDIL-SNU/SPINNER

1.2 Python requirements
=======================
SPINNER supports Python :code:`3.6` or higher version. SPINNER utilizes Python modules in the following:

  - mpi4py (https://bitbucket.org/mpi4py/mpi4py/)
  - PyYAML (https://pypi.org/project/PyYAML/)
  - numpy (https://numpy.org/)
  - pybind11 (https://github.com/pybind/pybind11/)

You can download these modules by

::

  pip3 install mpi4py PyYAML numpy pybind11

or directly downloading the source codes from the webpages.

1.3 C++ compilers
=================
SPINNER supports CMAKE 3.20 or higher version and gcc version 7.3.1 or higer version.

1.4 Python binding of randSpg and LAMMPS (important)
====================================================
SPINNER utilizes python binding of randSpg and LAMMPS code. However, those codes are slightly modified to be incorporated with SPINNER, so you need to comiple with source codes provided with SPINNER, not with the original ones.

1.4.1 Install randSpg
---------------------
To install randSpg, do follows

::

  cd /SPINNER-directory/randSpg-vspinner/
  mkdir build
  cd build
  cmake ..
  make -j3

to bind randSpg with python, do follows

::

  cd /SPINNER-directory/randSpg-vspinner/python
  mkdir build
  cd build
  cmake ..
  make –j3
  cp pyrandspg.cpython* /directory-where-your-python-is/lib/python3/site-packages/
 
1.4.2 Install LAMMPS
--------------------
This is the most tricky part of the installation steps. You can install python-lammps by follows but it may not work depending on your system. Please look into LAMMPS forum (https://www.lammps.org/forum.html) or LAMMPS manual Ch. 2 (https://docs.lammps.org/Python_run.html) for detail discussion.

::

  cd /SPINNER-directory/lammps-vspinner
  cd src
  make yes-python
  make XXX mode=shlib
  make install-python

Here, XXX is the name of the make file (Makefile.XXX in src/MAKE directory). Note that LAMMPS have to be installed with serial version (not mpi version). The optimization is recommended for your system but default serial option is also possible.

