.. _properties ::

Property API
==================================

This an api for....

Primary objects
---------------

.. currentmodule:: openforcefield.properties
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    PhysicalProperty

Built-in Properties
-------------------

.. currentmodule:: openforcefield.properties
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    Density
    DielectricConstant

Workflow Components
-------------------

Schema

.. currentmodule:: openforcefield.properties.workflow
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    CalculationSchema

Protocol API

.. currentmodule:: openforcefield.properties.workflow.protocols
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    BaseProtocol
    ProtocolSchema
    ProtocolPath

Built in Protocols

.. currentmodule:: openforcefield.properties.workflow.protocols
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    BuildCoordinatesPackmol
    BuildSmirnoffTopology
    RunEnergyMinimisation
    RunOpenMMSimulation
    AveragePropertyProtocol
    AverageTrajectoryProperty
    ExtractUncorrelatedData
    ExtractUncorrelatedTrajectoryData

Protocol Groups

.. currentmodule:: openforcefield.properties.workflow.groups
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ProtocolGroup
    ConditionalGroup