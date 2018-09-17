.. _topology ::

Molecular topology representations
==================================

This module provides pure-Python classes for representing molecules and molecular systems.
These classes offer several advantages over corresponding ``Topology`` objects in `OpenMM <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.topology.Topology.html#simtk.openmm.app.topology.Topology>`_ and `MDTraj <http://mdtraj.org/latest/api/generated/mdtraj.Topology.html#mdtraj.Topology>`_,
including offering serialization to a variety of standard formats (including `XML <https://www.w3.org/XML/>`_, `JSON <https://www.json.org/>`_, `YAML <http://yaml.org/>`_, `BSON <http://bsonspec.org/>`_, `TOML <https://github.com/toml-lang/toml>`_, and `MessagePack <https://msgpack.org/index.html>`_).


Primary objects
---------------

.. currentmodule:: openforcefield.topology
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    FrozenMolecule
    Molecule
    Topology

Secondary objects
-----------------

.. currentmodule:: openforcefield.topology
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    Particle
    Atom
    Bond
    VirtualSite
