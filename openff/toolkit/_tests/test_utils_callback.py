"""
Tests for callback utility classes and functions.

"""

import pytest

from openff.toolkit.utils.callback import Callbackable, callback_method
from openff.toolkit.utils.exceptions import CallbackRegistrationError


class CallHistory:
    """Used to keep track of the order in which callbacks and methods are called by Callbackable."""

    history = None

    def reset_history(self):
        CallHistory.history = []

    @classmethod
    def add_history_entry(cls, name, *args, **kwargs):
        # Store args and kwargs in the history only if they are given.
        history_entry = []
        if len(args) != 0:
            history_entry.append(args)
        if len(kwargs) != 0:
            history_entry.append(kwargs)

        if len(history_entry) == 0:
            cls.history.append(name)
        else:
            cls.history.append([name, *history_entry])

    def instance_callback(self, callbackable, event_name, *args, **kwargs):
        assert isinstance(self, object)
        CallHistory.add_history_entry("callback_" + event_name, *args, **kwargs)

    @classmethod
    def class_callback(cls, callbackable, event_name, *args, **kwargs):
        assert isinstance(cls, type)
        CallHistory.add_history_entry("callback_" + event_name, *args, **kwargs)

    @staticmethod
    def static_callback(callbackable, event_name, *args, **kwargs):
        CallHistory.add_history_entry("callback_" + event_name, *args, **kwargs)


call_history = CallHistory()


class TestCallbackable:
    """Test suite for the Callbackable base class."""

    class MyCallbackable(Callbackable):
        @callback_method
        def instance_method(self, *args, **kwargs):
            CallHistory.add_history_entry("instance_method", *args, **kwargs)

        @callback_method
        def __iadd__(self, other):
            CallHistory.add_history_entry("__iadd__", other)

        @callback_method(events=["event1"])
        def event_method1(self, *args, **kwargs):
            CallHistory.add_history_entry("event_method1", *args, **kwargs)

        @callback_method(events=["event1", "event2"])
        def event_method2(self, *args, **kwargs):
            CallHistory.add_history_entry("event_method2", *args, **kwargs)

    def check_method_call_order(
        self, callbackable, event_name, event_sequence, *args, **kwargs
    ):
        """Check that callback and methods/attributes are invoked in the correct order.

        This also formats the history correctly if args and kwargs are given.
        """
        # Modify the expected history if args and kwargs are given.
        if len(args) == 0 and len(kwargs) == 0:
            expected_history = event_sequence
        else:
            expected_history = [[event_name] for event_name in event_sequence]
        if len(args) != 0:
            for event in expected_history:
                event.append(args)
        if len(kwargs) != 0:
            for event in expected_history:
                event.append(kwargs)

        # Reset history and verify that the calls are in the correct order.
        call_history.reset_history()
        getattr(callbackable, event_name)(*args, **kwargs)
        assert call_history.history == expected_history

    @pytest.mark.parametrize("event_name", ["instance_method"])
    @pytest.mark.parametrize(
        "callback",
        [
            call_history.instance_callback,
            CallHistory.class_callback,
            CallHistory.static_callback,
        ],
    )
    @pytest.mark.parametrize(
        "args,kwargs", [([], {}), ([1, 2.0], {"kwarg1": 0, "kwarg2": None})]
    )
    def test_register_method_callback(self, event_name, callback, args, kwargs):
        """Methods' callbacks are invoked in the correct order and with the correct arguments."""
        callbackable = TestCallbackable.MyCallbackable()

        # No callback is called before registration.
        event_sequence = [event_name]
        self.check_method_call_order(
            callbackable, event_name, event_sequence, *args, **kwargs
        )

        # Register the callback.
        callbackable.register_callback(event_name, callback)

        # After the registration, the callback is invoked correctly.
        event_sequence = [event_name, "callback_" + event_name]
        self.check_method_call_order(
            callbackable, event_name, event_sequence, *args, **kwargs
        )

    def test_register_magic_method_callback(self):
        """Callbacks registered to magic methods are invoked correctly."""
        callbackable = TestCallbackable.MyCallbackable()
        callbackable.register_callback("__iadd__", call_history.instance_callback)

        extension = [1, 2]
        call_history.reset_history()
        callbackable += extension
        assert call_history.history == [
            ["__iadd__", (extension,)],
            ["callback___iadd__", (extension,)],
        ]

    def test_register_event_callback(self):
        """Callbacks registered to a event are handled corectly."""
        callbackable = TestCallbackable.MyCallbackable()

        # Register the callbacks to event1 (event_method1 and event_method2).
        callbackable.register_callback("event1", call_history.instance_callback)
        callbackable.register_callback("event1", CallHistory.class_callback)
        # Register one callback to event2 (only event_method2).
        callbackable.register_callback("event2", CallHistory.static_callback)

        # Check the event sequence for both methods belong to the two events.
        event_sequence = [
            "event_method1",
            "callback_event_method1",
            "callback_event_method1",
        ]
        self.check_method_call_order(callbackable, "event_method1", event_sequence)

        event_sequence = [
            "event_method2",
            "callback_event_method2",
            "callback_event_method2",
            "callback_event_method2",
        ]
        self.check_method_call_order(callbackable, "event_method2", event_sequence)

    def test_not_callback_method_raise_exception(self):
        """An exception is raised if a callback is registered for a method not tagged with callback_method."""

        class TempCallbackable(Callbackable):
            def not_callback_method(self):
                pass

        callbackable = TempCallbackable()
        with pytest.raises(
            CallbackRegistrationError,
            match="is not tagged with the @callback_method decorator",
        ):
            callbackable.register_callback(
                "not_callback_method", call_history.instance_callback
            )

    def test_unknown_event_raise_exception(self):
        """An exception is raised if a callback is registered for an unknown callback event."""
        callbackable = TestCallbackable.MyCallbackable()
        with pytest.raises(
            CallbackRegistrationError,
            match='is associated to the callback event "unknown"',
        ):
            callbackable.register_callback("unknown", call_history.instance_callback)
