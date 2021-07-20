#!/usr/bin/env python

"""
Utility classes and functions to create objects supporting callback registration.

"""

__all__ = [
    "callback_method",
    "Callbackable",
]


import functools
import inspect

# =====================================================================
# GLOBAL IMPORTS
# =====================================================================
from openff.toolkit.utils.exceptions import CallbackRegistrationError

# =====================================================================
# CALLBACKABLE CLASSES
# =====================================================================


def callback_method(func=None, events=()):
    """Decorator used to mark a method as callbackable.

    This decorator is designed to be used by ``Callbackable`` children
    classes. Usage examples can be found in the documentation of ``Callbackable``.

    Parameters
    ----------
    func : MethodType
        The instance method to decorate. Static and class methods are
        currently not supported.
    events : Iterable[str]
        The name of the callback events associated to this method.

    See Also
    --------
    Callbackable

    """

    if func is None and len(events) > 0:
        # Decorator called with @callback_method(events=[...]).
        return functools.partial(callback_method, events=events)

    @functools.wraps(func)
    def callbacked_func(self, *args, **kwargs):
        # Store the returned value.
        returned_value = func(self, *args, **kwargs)
        # Invoke callbacks that must be run after the function is executed.
        self._raise_callback_events(func.__name__, *args, **kwargs)
        return returned_value

    # Assign an attribute so that we can distinguish between
    # callbackable and not callbackable methods.
    callbacked_func._callback_events = set(events)

    return callbacked_func


class Callbackable:
    """A base class for registering and handling callbacks.

    The class allows easy implementation of classes that support the
    registration of callback function when some instance methods are
    called. Callback registered to static and class methods are currently
    not supported.

    The class allows callbacks to be registered to single methods or to
    an event, which essentially provides a way to register the same callback
    to multiple methods with a single line.

    Callbacks are called after the function has been executed, and they
    must have the signature

        callback(callbackable, function_name, *args, **kwargs)

    where ``callbackable`` is the ``Callbackable`` instance that raised
    the event, ``function_name`` is the name of the function of ``callbackable``
    that was called, and ``*args`` and ``**kwargs`` are the arguments that
    were passed to the function.

    Examples
    --------

    We define a callbackable list that allows registering callback
    functions to methods that change the list. The first four methods
    summarize the different ways in which the ``callback_method``
    decorator can be used to mark a method for callback registration.

    >>> class CallbackableList(Callbackable, list):
    ...
    ...     @callback_method
    ...     def pop(self, *args, **kwargs):
    ...         super().pop(*args, **kwargs)
    ...
    ...     @callback_method(events=['new_element'])
    ...     def extend(self, *args, **kwargs):
    ...         super().extend(*args, **kwargs)
    ...
    ...     remove = callback_method(list.remove)
    ...     append = callback_method(list.append, events=['new_element'])
    ...
    ...     insert = callback_method(list.insert, events=['new_element'])
    ...     __iadd__ = callback_method(list.__iadd__, events=['new_element'])
    ...     __setitem__ = callback_method(list.__setitem__, events=['new_element'])
    ...     __imul__ = callback_method(list.__imul__)
    ...     __delitem__ = callback_method(list.__delitem__)
    ...

    We can now register a callback to any of the callback methods.

    >>> def callback(callbackable, func_name, *args, **kwargs):
    ...     print(f'callback: {func_name}{args}')
    ...
    >>> l = CallbackableList([1, 2, 3])
    >>> l.register_callback('__delitem__', callback=callback)
    >>> del l[2]
    callback: __delitem__(2,)

    We can also register a callback to a generic event associated to one
    or more methods.

    >>> l.register_callback('new_element', callback=callback)
    >>> l.append(5)
    callback: append(5,)
    >>> l.insert(1, 5)
    callback: insert(1, 5)

    """

    def __init__(self, *args, **kwargs):
        # Map event_name -> list of callbacks.
        self._callbacks = {}
        # Forward arguments to the next class of the MRO to support multiple inheritance.
        super().__init__(*args, **kwargs)

    def register_callback(self, event_name, callback):
        # Check if a callback can be registered for this function. Methods
        # that are not decorated by @callback_method raise an error.
        try:
            attr = getattr(self, event_name)
        except AttributeError:
            # event_name may be an event name. This raises an exception
            # if no method associated to this event is found.
            self._check_event_exist(event_name)
        else:
            # Check if the attribute has been decorated.
            try:
                attr._callback_events
            except AttributeError:
                raise CallbackRegistrationError(
                    f"{self.__class__}.{event_name} is not "
                    f"tagged with the @callback_method decorator"
                )

        # Update the instance callbacks dictionary.
        try:
            self._callbacks[event_name].append(callback)
        except KeyError:
            self._callbacks[event_name] = [callback]

    def _check_event_exist(self, event_name):
        for member in inspect.getmembers(self):
            # If this is not a callback_method, skip.
            try:
                events = member[1]._callback_events
            except AttributeError:
                continue

            # If the event exist, return.
            if event_name in events:
                return

        # The event wasn't found.
        raise CallbackRegistrationError(
            f"No method of {self.__class__} is associated "
            f'to the callback event "{event_name}".'
        )

    def _raise_callback_events(self, func_name, *args, **kwargs):
        events = getattr(self, func_name)._callback_events

        # Retrieve all the callbacks associated to the function and its events.
        for event_name in [func_name, *events]:
            for callback in self._callbacks.get(event_name, ()):
                callback(self, func_name, *args, **kwargs)


if __name__ == "__main__":
    import doctest

    doctest.run_docstring_examples(Callbackable, globals())
