(topology)=

# Molecular topology representations

This module provides pure-Python classes for representing molecules and molecular
systems. These classes offer several advantages over corresponding `Topology`
objects in [OpenMM](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.topology.Topology.html#simtk.openmm.app.topology.Topology)
and [MDTraj](https://mdtraj.org/1.9.4/api/generated/mdtraj.Topology.html),
including offering serialization to a variety of standard formats (including [XML](https://www.w3.org/XML/), [JSON](https://www.json.org/), [YAML](http://yaml.org/), [BSON](http://bsonspec.org/), [TOML](https://github.com/toml-lang/toml), and [MessagePack](https://msgpack.org/index.html)).


## Primary objects

```{eval-rst}
.. currentmodule:: openff.toolkit.topology
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    FrozenMolecule
    Molecule
    Topology
```

## Secondary objects

```{eval-rst}
.. currentmodule:: openff.toolkit.topology
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    Particle
    Atom
    Bond
    ValenceDict
    ImproperDict
    HierarchyScheme
    HierarchyElement
```
