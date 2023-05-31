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

The :attr:`~endf.Material.section_text` attribute holds a dictionary mapping
(MF, MT) to the section of text from the ENDF file:

.. code-block:: pycon

    >>> print(mat.section_text[12, 75])
     9.223500+4 2.330248+2          2          2         25          0922812 75
     4.386000+5 0.000000+0          0          0          3          1922812 75
     0.000000+0 1.000000+0 1.000000+0                                 922812 75

The :attr:`~endf.Material.section_data` attribute also holds a dictionary with
keys that are (MF, MT) pairs, but the values are dictionaries that hold the
individual pieces of data from the ENDF file:

.. code-block:: pycon

    >>> mat.section_data[3, 16]
    {'ZA': 92235,
    'AWR': 233.0248,
    'QM': -5298000.0,
    'QI': -5298000.0,
    'LR': 0,
    'sigma': <Tabulated1D: 39 points, 1 regions>}

Note that indexing the :class:`~endf.Material` object directly is equivalent to
indexing the :attr:`~endf.Material.section_data` attribute:

.. code-block:: pycon

    >>> mat[3, 16]
    {'ZA': 92235,
    'AWR': 233.0248,
    'QM': -5298000.0,
    'QI': -5298000.0,
    'LR': 0,
    'sigma': <Tabulated1D: 39 points, 1 regions>}


Tabulated1D Objects
~~~~~~~~~~~~~~~~~~~

As the above example shows, whenever a TAB1 record is encountered in an ENDF
file, it is represented as a :class:`~endf.Tabulated1D` object. The raw `x` and
`y` values can be accessed through the :attr:`~endf.Tabulated1D.x` and
:attr:`~endf.Tabulated1D.y` attributes. For example, if you wanted to make a
log-log plot of a cross section::

    import matplotlib.pyplot as plt

    xs = mat.section_data[3, 16]['sigma']
    plt.loglog(xs.x, xs.y, marker='.', label='(n,2n)')
    plt.show()

The :class:`~endf.Tabulated1D` acts like a function; you can pass it a given
value of the independent variable and it will return the corresponding value of
the dependent variable, interpolating between tabulated points where necessary.
For example, if you wanted to know the cross section at 9.5 MeV:

.. code-block:: pycon

    >>> xs(9.5e6)
    0.7428240476190476

Passing a list, array, or other iterable of values will return a :mod:`numpy`
array of corresponding values:

.. code-block:: pycon

    >>> energies = numpy.linspace(6.0e6, 8.0e6, 10)
    >>> xs(energies)
    array([0.183841  , 0.25245628, 0.3152985 , 0.3550485 , 0.3947985 ,
           0.43081767, 0.463106  , 0.49255778, 0.52108125, 0.55290937])

.. _high_level_interface:

High-level Interface
++++++++++++++++++++

While this form is more useful, it still may be a little too "raw". The
:class:`~endf.Material` class has an :meth:`~endf.Material.interpret` method
that returns a class based on the ENDF sublibrary type. For example,
incident-neutron data will result in an instance of the
:class:`endf.IncidentNeutron` class. These "interpreted" classes provide a much
more intuitive interface to data within an ENDF (or ACE) file.

Incident Neutron Data
~~~~~~~~~~~~~~~~~~~~~

The :class:`~endf.IncidentNeutron` class provides the high-level interface to
ENDF incident neutron sublibrary files. You can get an instance of this class
either by calling the :meth:`endf.Material.interpret` method::

    mat = endf.Material('n-092_U_235.endf')
    u235 = mat.interpret()

or by directly passing a filename to the :meth:`endf.IncidentNeutron.from_endf`
method:

.. code-block:: pycon

    >>> u235 = endf.IncidentNeutron.from_endf('n-092_U_235.endf')
    >>> u235
    <IncidentNeutron: U235, 85 reactions>

Most incident neutron data is collected into the
:attr:`~endf.IncidentNeutron.reactions` attribute, which is a dictionary that
maps the MT value to an instance of the :class:`~endf.Reaction` class:

.. code-block:: pycon

    >>> u235.reactions
    {1: <Reaction: MT=1 (n,total)>,
     2: <Reaction: MT=2 (n,elastic)>,
     4: <Reaction: MT=4 (n,level)>,
     5: <Reaction: MT=5 (n,misc)>,
     16: <Reaction: MT=16 (n,2n)>,
     17: <Reaction: MT=17 (n,3n)>,
     18: <Reaction: MT=18 (n,fission)>,
     ...
     835: <Reaction: MT=835 (n,a35)>}

To look at data for a specific reaction then, you can index the
:attr:`~endf.IncidentNeutron.reactions` attribute:

.. code-block:: pycon

    >>> u235.reactions[16]
    <Reaction: MT=16 (n,2n)>

Alteratively, you can pass the MT value as an index to
:class:`~endf.IncidentNeutron` directly:

.. code-block:: pycon

    >>> u235[16]
    <Reaction: MT=16 (n,2n)>

or even use the name of the reaction as an index:

.. code-block:: pycon

    >>> u235['n,2n']
    <Reaction: MT=16 (n,2n)>

The :class:`~endf.Reaction` class has several attributes, including
:attr:`~endf.Reaction.xs` (cross section), :attr:`~endf.Reaction.products`
(reaction products), :attr:`~endf.Reaction.q_reaction` (reaction Q-value) and
:attr:`~endf.Reaction.q_massdiff` (mass-difference Q value). The
:attr:`~endf.Reaction.xs` attribute is a dictionary mapping a temperature to the
integral cross section:

.. code-block:: pycon

    >>> n2n = u235['n,2n']
    >>> n2n.xs
    {'0K': <Tabulated1D: 39 points, 1 regions>}

For data that originates from an ENDF file, the cross sections are always
present at 0 K. However, the :class:`~endf.IncidentNeutron` class can also be
used for ACE data. In that case, cross sections can be present at temperatures
other than 0 K.

The :attr:`~endf.Reaction.products` attribute gives a list of reaction products
as :class:`~endf.Product` objects:

.. code-block:: pycon

    >>> n2n.products
    [<Product: neutron, emission=prompt, yield=2.0>,
     <Product: U234, emission=prompt, yield=1.0>,
     <Product: photon, emission=prompt, yield=tabulated>]

The yield of a given product is accessed through the
:attr:`~endf.Product.yield_` attribute:

.. code-block:: pycon

    >>> photon = n2n.products[-1]
    >>> photon.yield_
    <Tabulated1D: 39 points, 1 regions>

ACE Files
+++++++++

Working with ACE files is conceptually similar to ENDF files. A low-level
interface provides access to the raw data within an ACE file (the `NXS`, `JXS`,
and `XSS` arrays) and the same high-level interface classes, e.g.,
:class:`endf.IncidentNeutron` can be used to more easily inspect data. If you
have an ACE file with a single table within it, you can load it with the
:func:`endf.ace.get_table` function, which returns a :class:`endf.ace.Table`
object:

.. code-block:: pycon

    >>> table = endf.ace.get_table('80198.710nc')
    >>> table
    <ACE Table: 80198.710nc at 293.6 K>

Raw access to the underlying arrays in the ACE file is provided via the
:attr:`~endf.ace.Table.nxs`, :attr:`~endf.ace.Table.jxs`, and
:attr:`~endf.ace.Table.xss` attributes:

.. code-block:: pycon

    >>> table.nxs
    array([     0, 172118,  80198,   3180,     34,     25,     97,      4,
                0,      0,     80,    198,      0,      0,      0,      0,
                0])

Note that each of these arrays is prepended with an extra zero at the beginning
so that the indexing follows the Fortran 1-based indexing that is referenced in
the ACE format manual.

As with the :class:`~endf.Material` class, the :class:`~endf.ace.Table` class
has a :meth:`~endf.ace.Table.interpret` that returns a corresponding high-level
interface class. For continuous-energy neutron data, this method will return an
:class:`~endf.IncidentNeutron` object:

.. code-block:: pycon

    >>> hg198 = table.interpret()
    >>> hg198
    <IncidentNeutron: Hg198, 38 reactions>

From this point, the interface follows exactly as is shown in
:ref:`high_level_interface`.
