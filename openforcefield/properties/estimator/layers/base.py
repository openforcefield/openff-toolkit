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
    def _await_results(backend, data_model, callback, submitted_futures):
        """A method to handle passing the results of this layer back to
        the main thread."""

        callback_future = backend.submit_task(return_args, data_model, *submitted_futures)

        def callback_wrapper(future_object):

            results = list(future_object.result())
            returned_data_model = results.pop(0)

            for succeeded, returned_property in results:

                if not succeeded:
                    continue

                returned_data_model.queued_properties.remove(returned_property)
                returned_data_model.calculated_properties.append(returned_property)

            callback(returned_data_model)

        callback_future.add_done_callback(callback_wrapper)

    @staticmethod
    def perform_calculation(backend, data_model, existing_data, callback):
        raise NotImplementedError()
