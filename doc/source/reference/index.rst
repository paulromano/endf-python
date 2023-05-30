.. _api:

API Reference
=============

.. module:: endf

Low-level Interface
-------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   Material

.. autosummary::
   :toctree: generated
   :nosignatures:

   get_materials

High-level Interface
--------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   IncidentNeutron
   Tabulated1D
   Reaction
   Product

ACE File Interface
------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   ace.Table
   ace.TableType

.. autosummary::
   :toctree: generated
   :nosignatures:

   ace.get_table
   ace.get_tables
   ace.get_libraries_from_xsdir
