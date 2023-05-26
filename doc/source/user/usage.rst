.. _usage:

Detailed Usage
--------------

The :mod:`endf` package provides both low-level and high-level interfaces to
ENDF files. The low-level interface can be used when you simply need access to
the raw data in an ENDF file and you know exactly what you are looking for
(usually with a copy of the ENDF-6 format manual open next to you). The
high-level interface can be used if you are a normal human and you just want to,
say, look at a cross section for a specific reaction.

Low-level Interface
+++++++++++++++++++

To read an ENDF file with a single material, simply pass the filename to the
:class:`~endf.Material` class:


.. code-block::

    import endf

    mat = endf.Material('n-092_U_235.endf')

The ``section_text`` attribute holds a dictionary mapping (MF, MT) to the section
of text from the ENDF file:

.. code-block:: pycon

    >>> print(mat.section_text[12, 75])
     9.223500+4 2.330248+2          2          2         25          0922812 75
     4.386000+5 0.000000+0          0          0          3          1922812 75
     0.000000+0 1.000000+0 1.000000+0                                 922812 75

The ``section_data`` attribute also holds a dictionary with keys that are (MF, MT)
pairs, but the values are dictionaries that hold the individual pieces of data
from the ENDF file:

.. code-block:: pycon

    >>> mat.section_data[3, 16]
    {'ZA': 92235,
    'AWR': 233.0248,
    'QM': -5298000.0,
    'QI': -5298000.0,
    'LR': 0,
    'sigma': <Tabulated1D: 39 points, 1 regions>}

High-level Interface
++++++++++++++++++++

While this form is more useful, it still may be a little too "raw". The
:class:`~endf.Material` class has an :meth:`~endf.Material.interpret` method
that returns a class based on the ENDF sublibrary type (for example,
incident-neutron data will result in an instance of the ``IncidentNeutron``
class).
