.. _property_estimator ::

Property Estimator
==================================

This module provides an api for multi-fidelity property calculations.

Client Side API
---------------

.. currentmodule:: openforcefield.properties.estimator.client
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    PropertyEstimator
    PropertyEstimatorOptions
    PropertyEstimatorDataModel

Server Side API
---------------

.. currentmodule:: openforcefield.properties.estimator.runner
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    PropertyCalculationRunner
    PropertyRunnerDataModel

Calculation Layers
------------------

.. currentmodule:: openforcefield.properties.estimator.layers
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    PropertyCalculationLayer
    register_calculation_layer
    SurrogateLayer
    ReweightingLayer
    SimulationLayer

Calculation Backends
--------------------

.. currentmodule:: openforcefield.properties.estimator.backends
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    PropertyEstimatorBackend
    DaskLocalClusterBackend