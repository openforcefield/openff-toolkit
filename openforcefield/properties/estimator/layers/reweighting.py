# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Direct simulation layer.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import pickle

from openforcefield.properties.estimator.layers.base import register_calculation_layer, PropertyCalculationLayer
from openforcefield.typing.engines.smirnoff import ForceField


# =============================================================================================
# Reweighting Layer
# =============================================================================================

@register_calculation_layer()
class ReweightingLayer(PropertyCalculationLayer):
    """A calculation layer which aims to calculate physical properties by
    reweighting the results of previous calculations.

    .. warning :: This class has not yet been implemented.
    """

    @staticmethod
    def schedule_calculation(calculation_backend, storage_backend, layer_directory,
                             data_model, callback, synchronous=False):

        parameter_set = storage_backend.retrieve_force_field(data_model.parameter_set_id)

        reweighting_futures = []

        for physical_property in data_model.queued_properties:

            existing_data = storage_backend.retrieve_simulation_data(str(physical_property.substance))

            reweighting_future = calculation_backend.submit_task(ReweightingLayer.perform_reweighting,
                                                                 physical_property,
                                                                 parameter_set,
                                                                 existing_data)

            reweighting_futures.append(reweighting_future)

        PropertyCalculationLayer._await_results(calculation_backend,
                                                storage_backend,
                                                layer_directory,
                                                data_model,
                                                callback,
                                                reweighting_futures,
                                                synchronous)

    @staticmethod
    def perform_reweighting(physical_property, parameter_set, existing_data, **kwargs):
        """A placeholder method that would be used to attempt
        to reweight previous calculations to yield the desired
        property.

        .. warning :: This method has not yet been implemented.

        Parameters
        ----------
        physical_property: PhysicalProperty
            The physical property to attempt to estimate by reweighting.
        parameter_set: ForceField
            The force field parameters to use when estimating the property.
        existing_data: list of StoredSimulationData
            Data which has been stored from previous calculations on systems
            of the same composition as the desired property.
        """

        # A return value indicates that the reweighting layer did not
        # have access to enough information to accurately estimate the property.
        return None
