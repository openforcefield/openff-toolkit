#!/usr/bin/env python

"""
Utility classes and functions to create classes supporting callbacks.

"""

__all__ = [
    'Callbackable',
    'callback_method'
    'CallbackRegistrationError'
]


# =====================================================================
# GLOBAL IMPORTS
# =====================================================================

import functools
import inspect


# =====================================================================
# CALLBACKABLE
# =====================================================================

class CallbackRegistrationError(TypeError):
    """Error raised during callback registration."""
    pass


def callback_method(func=None, event_groups=()):

    if func is None and len(event_groups) > 0:
        # Decorator called with @callback_method(event_groups=[...]).
        return functools.partial(callback_method, event_groups=event_groups)

    @functools.wraps(func)
    def callbacked_func(self, *args, **kwargs):
        # Invoke callbacks that must be run before the function.
        self._raise_callback_events(func.__name__, *args, before=True, **kwargs)
        # Store the returned value.
        returned_value = func(self, *args, **kwargs)
        # Invoke callbacks that must be run after the function is executed.
        self._raise_callback_events(func.__name__, *args, before=False, **kwargs)
        return returned_value

    # Assign an attribute so that we can distinguish between
    # callbackable and not callbackable methods.
    callbacked_func._callback_event_groups = set(event_groups)

    return callbacked_func


class Callbackable:
    """A base class for registering and handling callbacks."""

    def __init__(self, *args, **kwargs):
        # Map event_name -> list of callbacks.
        self._before_callbacks = {}
        self._after_callbacks = {}
        # Forward arguments to the next class of the MRO to support multiple inheritance.
        super().__init__(*args, **kwargs)

    def register_callback(self, event_name, callback, before=False):
        # Check if a callback can be registered for this function. Methods
        # that are not decorated by @callback_method raise an error.
        try:
            attr = getattr(self, event_name)
        except AttributeError:
            # event_name may be a group name.
            self._check_event_group_exist(event_name)
        else:
            # Check if the attribute has been decorated.
            try:
                attr._callback_event_groups
            except AttributeError:
                raise CallbackRegistrationError(f'{self.__class__}.{event_name} is not '
                                                f'tagged with the @callback_method decorator')

        # Determine if this callback must be called before or after calling the method.
        callbacks = self._get_callbacks(before)

        # Update the instance callbacks dictionary.
        try:
            callbacks[event_name].append(callback)
        except KeyError:
            callbacks[event_name] = [callback]

    def _check_event_group_exist(self, event_group_name):
        for member in inspect.getmembers(self):
            # If this is not a callback_method, skip.
            try:
                event_groups = member[1]._callback_event_groups
            except AttributeError:
                continue

            # If the group exist, return.
            if event_group_name in event_groups:
                return

        # The group wasn't found.
        raise CallbackRegistrationError(f'No method of {self.__class__} is associated '
                                        f'to the callback event group "{event_group_name}".')

    def _get_callbacks(self, before):
        if before:
            return self._before_callbacks
        return self._after_callbacks

    def _raise_callback_events(self, func_name, *args, before=False, **kwargs):
        event_groups = getattr(self, func_name)._callback_event_groups

        # Retrieve all the callbacks associated to the function and its groups.
        for event_name in [func_name, *event_groups]:
            event_callbacks = self._get_callbacks(before).get(event_name, [])
            for callback in event_callbacks:
                callback(self, func_name, *args, **kwargs)
