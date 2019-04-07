.. _typing ::

Forcefield typing tools
=======================

Chemical environments
---------------------

Tools for representing and operating on chemical environments

    .. warning :: This class is largely redundant with the same one in the Chemper package, and will likely be removed.


.. currentmodule:: openforcefield.typing.chemistry
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ChemicalEnvironment


Forcefield typing engines
-------------------------

Engines for applying parameters to chemical systems

The SMIRks-Native Open Force Field (SMIRNOFF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A reference implementation of the SMIRNOFF specification for parameterizing biomolecular systems

ForceField
~~~~~~~~~~

The ``ForceField`` class is a primary part of the top-level toolkit API.
``ForceField`` objects are initialized from SMIRNOFF data sources (e.g. an ``OFFXML`` file).
For a basic example of system creation using a ``ForceField``, see ``examples/SMIRNOFF_simulation``.


.. currentmodule:: openforcefield.typing.engines.smirnoff.forcefield
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ForceField

Parameter Type
~~~~~~~~~~~~~~

``ParameterType`` objects are representations of individual SMIRKS-based SMIRNOFF parameters.
These are usually initialized during ``ForceField`` creation, and can be inspected and modified by users via the Python API.
For more information, see ``examples/forcefield_modification``.

.. currentmodule:: openforcefield.typing.engines.smirnoff.parameters
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ParameterType
    BondHandler.BondType
    AngleHandler.AngleType
    ProperTorsionHandler.ProperTorsionType
    ImproperTorsionHandler.ImproperTorsionType
    vdWHandler.vdWType

Parameter Handlers
~~~~~~~~~~~~~~~~~~

Each ``ForceField`` primarily consists of several ``ParameterHandler`` objects, which each contain the machinery to add one energy component to a system.
During system creation, each ``ParameterHandler`` registered to a ``ForceField`` has its ``assign_parameters()`` function called..

.. currentmodule:: openforcefield.typing.engines.smirnoff.parameters
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ParameterList
    ParameterHandler
    BondHandler
    AngleHandler
    ProperTorsionHandler
    ImproperTorsionHandler
    vdWHandler
    ElectrostaticsHandler
    ToolkitAM1BCCHandler

Parameter I/O Handlers
~~~~~~~~~~~~~~~~~~~~~~

``ParameterIOHandler`` objects handle reading and writing of serialzied SMIRNOFF data sources.

.. currentmodule:: openforcefield.typing.engines.smirnoff.io
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ParameterIOHandler
    XMLParameterIOHandler