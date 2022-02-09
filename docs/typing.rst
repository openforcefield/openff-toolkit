.. _typing :

Force field typing tools
========================

Chemical environments
---------------------

Tools for representing and operating on chemical environments

.. currentmodule:: openff.toolkit.typing.chemistry
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ChemicalEnvironment


Force field typing engines
--------------------------

Engines for applying parameters to chemical systems

The SMIRks-Native Open Force Field (SMIRNOFF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A reference implementation of the SMIRNOFF specification for parameterizing biomolecular systems

ForceField
~~~~~~~~~~

The ``ForceField`` class is a primary part of the top-level toolkit API.
``ForceField`` objects are initialized from SMIRNOFF data sources (e.g. an ``OFFXML`` file).
For a basic example of OpenMM ``System`` creation using a ``ForceField``, see ``examples/SMIRNOFF_simulation``.


.. currentmodule:: openff.toolkit.typing.engines.smirnoff.forcefield
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ForceField
    get_available_force_fields

Parameter Type
~~~~~~~~~~~~~~

``ParameterType`` objects are representations of individual SMIRKS-based SMIRNOFF parameters.
These are usually initialized during ``ForceField`` creation, and can be inspected and modified by users via the Python API.
For more information, see ``examples/forcefield_modification``.

.. currentmodule:: openff.toolkit.typing.engines.smirnoff.parameters
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ParameterType
    ConstraintHandler.ConstraintType
    BondHandler.BondType
    AngleHandler.AngleType
    ProperTorsionHandler.ProperTorsionType
    ImproperTorsionHandler.ImproperTorsionType
    vdWHandler.vdWType
    LibraryChargeHandler.LibraryChargeType
    GBSAHandler.GBSAType
    ChargeIncrementModelHandler.ChargeIncrementType
    VirtualSiteHandler.VirtualSiteBondChargeType
    VirtualSiteHandler.VirtualSiteMonovalentLonePairType
    VirtualSiteHandler.VirtualSiteDivalentLonePairType
    VirtualSiteHandler.VirtualSiteTrivalentLonePairType

Parameter Handlers
~~~~~~~~~~~~~~~~~~

Each ``ForceField`` primarily consists of several ``ParameterHandler`` objects, which each contain the machinery to add one energy component to an OpenMM ``System``.
During ``System`` creation, each ``ParameterHandler`` registered to a ``ForceField`` has its ``assign_parameters()`` function called.

.. currentmodule:: openff.toolkit.typing.engines.smirnoff.parameters
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ParameterList
    ParameterHandler
    ConstraintHandler
    BondHandler
    AngleHandler
    ProperTorsionHandler
    ImproperTorsionHandler
    vdWHandler
    ElectrostaticsHandler
    LibraryChargeHandler
    ToolkitAM1BCCHandler
    GBSAHandler
    ChargeIncrementModelHandler
    VirtualSiteHandler

Parameter I/O Handlers
~~~~~~~~~~~~~~~~~~~~~~

``ParameterIOHandler`` objects handle reading and writing of serialzied SMIRNOFF data sources.

.. currentmodule:: openff.toolkit.typing.engines.smirnoff.io
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ParameterIOHandler
    XMLParameterIOHandler



Parameter Attributes
~~~~~~~~~~~~~~~~~~~~

``ParameterAttribute`` and ``IndexedParameterAttribute`` provide a standard backend for ParameterHandler and Parameter attributes, while also enforcing validation of types and units.

.. currentmodule:: openff.toolkit.typing.engines.smirnoff.parameters
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ParameterAttribute
    IndexedParameterAttribute
    MappedParameterAttribute
    IndexedMappedParameterAttribute
