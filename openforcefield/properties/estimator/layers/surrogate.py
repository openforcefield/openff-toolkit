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
    def schedule_calculation(backend, data_model, existing_data, callback, synchronous=False):

        parameter_set = ForceField([])

        with open(data_model.parameter_set_path, 'rb') as file:
            parameter_set.__setstate__(pickle.load(file))

        surrogate_futures = []

        for physical_property in data_model.queued_properties:

            surrogate_future = backend.submit_task(SurrogateLayer.perform_surrogate_extrapolation,
                                                   physical_property,
                                                   parameter_set)

            surrogate_futures.append(surrogate_future)

        PropertyCalculationLayer._await_results(backend, data_model, callback, surrogate_futures, synchronous)

    @staticmethod
    def perform_surrogate_extrapolation(physical_property, parameter_set):
        """A placeholder method that would be used to spawn the surrogate
        model backend.

        .. warning :: This method has not yet been implemented.
        """

        # For now the return tuple indicates that the surrogate modelling
        # was not sufficiently accurate to estimate the property (False)
        # and simply returns the property back to be passed to the next layer.
        return False, physical_property
