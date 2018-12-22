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

from openforcefield.properties.estimator.layers.base import register_calculation_layer, PropertyCalculationLayer


# =============================================================================================
# Surrogate Layer
# =============================================================================================

@register_calculation_layer()
class SurrogateLayer(PropertyCalculationLayer):
    """Attempts to calculate properties from a precomputed
    surrogate model.
    """

    @staticmethod
    def perform_calculation(backend, data_model, existing_data, callback):
        """Attempt to calculate properties by surrogate modelling.
        """
        surrogate_futures = []

        for physical_property in data_model.queued_properties:

            surrogate_future = backend.submit_task(SurrogateLayer.perform_surrogate_extrapolation,
                                                   physical_property,
                                                   data_model.parameter_set)

            surrogate_futures.append(surrogate_future)

        PropertyCalculationLayer._await_results(backend, data_model, callback, surrogate_futures)

    @staticmethod
    def perform_surrogate_extrapolation(physical_property, parameter_set):
        """A placeholder method that would be used to spawn the surrogate
        model backend.
        """

        # For now the return tuple indicates that the surrogate modelling
        # was not sufficiently accurate to estimate the property (False)
        # and simply returns the property back to be passed to the next layer.
        return False, physical_property
