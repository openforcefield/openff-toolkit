# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Base calculation layer API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import logging

from openforcefield.properties import PhysicalProperty
from openforcefield.properties.estimator.storage import StoredSimulationData
from openforcefield.properties.estimator.workflow.protocols import PropertyCalculatorException

# =============================================================================================
# Registration Decorators
# =============================================================================================

available_layers = {}


def register_calculation_layer():
    """A decorator which registers a class as being a calculation layer
    which may be used in property calculations.
    """

    def decorator(cls):

        if cls.__name__ in available_layers:
            raise ValueError('The {} layer is already registered.'.format(cls.__name__))

        available_layers[cls.__name__] = cls
        return cls

    return decorator


def return_args(*args):
    return args


# =============================================================================================
# Base Layer
# =============================================================================================

class CalculationLayerResult:
    """The output returned from attempting to calculate a property on
     a PropertyCalculationLayer."""

    def __init__(self):
        """Constructs a new CalculationLayerResult object.
        """
        self.property_id: str = None

        self.calculated_property: PhysicalProperty = None
        self.calculation_error: PropertyCalculatorException = None

        self.data_to_store: StoredSimulationData = None


class PropertyCalculationLayer:
    """An abstract representation of a calculation layer in the property calculation stack.

    Notes
    -----
    Calculation layers must inherit from this class, and must override the
    `schedule_calculation` method.
    """

    @staticmethod
    def _await_results(calculation_backend, storage_backend, layer_directory, data_model,
                       callback, submitted_futures, synchronous=False):
        """A helper method to handle passing the results of this layer back to
        the main thread.

        Parameters
        ----------
        calculation_backend: openforcefield.properties.estimator.backends.PropertyEstimatorBackend
            The backend to the submit the calculations to.
        storage_backend: openforcefield.properties.estimator.storage.PropertyEstimatorStorage
            The backend used to store / retrieve data from previous calculations.
        layer_directory: str
            The local directory in which to store all local, temporary calculation data from this layer.
        data_model: openforcefield.properties.estimator.runner.PropertyRunnerDataModel
            The data model encoding the awaited calculation.
        callback: function
            The function to call when the backend returns the results (or an error).
        submitted_futures: list(dask.distributed.Future)
            A list of the futures returned by the backed when submitting the calculation.
        synchronous: bool
            If true, this function will block until the calculation has completed.
        """

        callback_future = calculation_backend.submit_task(return_args, data_model, *submitted_futures)

        def callback_wrapper(future_object):

            results = list(future_object.result())
            returned_data_model = results.pop(0)

            for returned_output in results:

                if returned_output is None:
                    # Indicates the layer could not calculate this
                    # particular property.
                    continue

                return_object = None

                if returned_output.calculation_error is not None:
                    return_object = returned_output.calculation_error
                else:
                    return_object = returned_output.calculated_property

                # Make sure to store any important calculation data.
                if returned_output.data_to_store is not None and returned_output.calculated_property is not None:

                    if returned_output.data_to_store.parameter_set_id is None:
                        returned_output.data_to_store.parameter_set_id = data_model.parameter_set_id

                    substance_id = str(returned_output.calculated_property.substance)

                    storage_backend.store_simulation_data(substance_id, returned_output.data_to_store)

                matches = [x for x in returned_data_model.queued_properties if x.id == returned_output.property_id]

                if len(matches) != 1:
                    logging.info('An id conflict occurred... unexpected results may ensue.')

                for match in matches:
                    returned_data_model.queued_properties.remove(match)

                returned_data_model.calculated_properties.append(return_object)

            callback(returned_data_model)

        if synchronous:
            callback_wrapper(callback_future)
        else:
            callback_future.add_done_callback(callback_wrapper)

    @staticmethod
    def schedule_calculation(calculation_backend, storage_backend, layer_directory,
                             data_model, callback, synchronous=False):
        """Submit the proposed calculation to the backend of choice.

        Parameters
        ----------
        calculation_backend: openforcefield.properties.estimator.backends.PropertyEstimatorBackend
            The backend to the submit the calculations to.
        storage_backend: openforcefield.properties.estimator.storage.PropertyEstimatorStorage
            The backend used to store / retrieve data from previous calculations.
        layer_directory: str
            The local directory in which to store all local, temporary calculation data from this layer.
        data_model: openforcefield.properties.estimator.runner.PropertyRunnerDataModel
            The data model encoding the proposed calculation.
        callback: function
            The function to call when the backend returns the results (or an error).
        synchronous: bool
            If true, this function will block until the calculation has completed.
            This is mainly intended for debugging purposes.
        """
        raise NotImplementedError()
