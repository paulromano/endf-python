.. _getting_started:

Getting Started
---------------

What is ``endf``?
+++++++++++++++++

:mod:`endf` is a Python package for reading and interpreting `ENDF-6
<https://doi.org/10.2172/1425114>`_ and `ACE
<https://github.com/NuclearData/ACEFormat>`_ format nuclear data files. Compared
to other packages that provide functionality for working with ENDF and ACE
files, this package has numerous advantages:

- Easily installable through ``pip``
- Thoughtful API design targeting Python first
- Offers both a low-level interface for working with the raw data in an ENDF
  file as well as a more intuitive high-level interface
- Fast file-loading performance thanks to optimized read routines
- Fully documented, tested, and type-hinted

Installation
++++++++++++

The ``endf`` package can be installed by running:

.. code-block:: sh

    python -m pip install endf

