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

.. currentmodule:: openforcefield.typing.engines.smirnoff.forcefield
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ForceField

.. currentmodule:: openforcefield.typing.engines.smirnoff.parameters
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ParameterType
    ParameterList
    ParameterHandler
    BondHandler
    AngleHandler
    ProperTorsionHandler
    ImproperTorsionHandler
    vdWHandler
    ElectrostaticsHandler
    ToolkitAM1BCCHandler

.. currentmodule:: openforcefield.typing.engines.smirnoff.io
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ParameterIOHandler