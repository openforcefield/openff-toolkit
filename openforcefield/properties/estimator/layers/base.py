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

from pydantic import BaseModel

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

class PropertyCalculationLayer:
    """Represents a layer in the property calculation stack.
    """

    @staticmethod
    def _await_results(backend, data_model, callback, submitted_futures, synchronous=False):
        """A method to handle passing the results of this layer back to
        the main thread.

        Parameters
        ----------
        backend: openforcefield.properties.estimator.backends.base.PropertyEstimatorBackend
            The backend to the submit the calculations to.
        data_model: openforcefield.properties.estimator.runner.PropertyRunnerDataModel
            The data model encoding the awaited calculation.
        callback: function
            The function to call when the backend returns the results (or an error).
        submitted_futures: list(dask.distributed.Future)
            A list of the futures returned by the backed when submitting the calculation.
        synchronous: bool
            If true, this function will block until the calculation has completed.
        """

        callback_future = backend.submit_task(return_args, data_model, *submitted_futures)

        def callback_wrapper(future_object):

            results = list(future_object.result())
            returned_data_model = results.pop(0)

            for succeeded, returned_output in results:

                if not succeeded:
                    continue

                matches = [x for x in returned_data_model.queued_properties if x.id == returned_output.id]

                if len(matches) != 1:
                    logging.info('An id conflict occurred... unexpected results may ensue.')

                for match in matches:
                    returned_data_model.queued_properties.remove(match)

                returned_data_model.calculated_properties.append(returned_output)

            callback(returned_data_model)

        if synchronous:
            callback_wrapper(callback_future)
        else:
            callback_future.add_done_callback(callback_wrapper)

    @staticmethod
    def schedule_calculation(backend, data_model, existing_data, callback, synchronous=False):
        """Submit the proposed calculation to the backend of choice.

        Parameters
        ----------
        backend: openforcefield.properties.estimator.backends.base.PropertyEstimatorBackend
            The backend to the submit the calculations to.
        data_model: openforcefield.properties.estimator.runner.PropertyRunnerDataModel
            The data model encoding the proposed calculation.
        existing_data: dict of str and Any
            Data which has already been calculated by a previous layer.
        callback: function
            The function to call when the backend returns the results (or an error).
        synchronous: bool
            If true, this function will block until the calculation has completed.
            This is mainly intended for debugging purposes.
        """
        raise NotImplementedError()
