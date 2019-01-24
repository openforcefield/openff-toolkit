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
# Surrogate Layer
# =============================================================================================

@register_calculation_layer()
class SurrogateLayer(PropertyCalculationLayer):
    """A calculation layer which aims to calculate physical properties from
    a surrogate model, such as a Gaussian mixture model.

    .. warning :: This class has not yet been implemented.
    """

    @staticmethod
    def schedule_calculation(calculation_backend, storage_backend, layer_directory,
                             data_model, callback, synchronous=False):

        parameter_set = storage_backend.retrieve_force_field(data_model.parameter_set_id)

        surrogate_futures = []

        for physical_property in data_model.queued_properties:

            surrogate_future = calculation_backend.submit_task(SurrogateLayer.perform_surrogate_extrapolation,
                                                               physical_property,
                                                               parameter_set)

            surrogate_futures.append(surrogate_future)

        PropertyCalculationLayer._await_results(calculation_backend,
                                                storage_backend,
                                                layer_directory,
                                                data_model,
                                                callback,
                                                surrogate_futures,
                                                synchronous)

    @staticmethod
    def perform_surrogate_extrapolation(physical_property, parameter_set):
        """A placeholder method that would be used to spawn the surrogate
        model backend.

        .. warning :: This method has not yet been implemented.
        """

        # A return value indicates that the surrogate layer did not
        # have access to enough information to accurately estimate the property.
        return None
