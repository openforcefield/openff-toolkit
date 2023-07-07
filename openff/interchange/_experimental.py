import os
from functools import wraps
from inspect import cleandoc

from openff.interchange.exceptions import ExperimentalFeatureException


def experimental(func):
    """
    Decorate a function as experimental, requiring environment variable opt-in.

    To use the wrapped function, set the environment variable INTERCHANGE_EXPERIMENTAL=1.

    """
    firstline, _, remainder = cleandoc(func.__doc__).partition("\n")

    func.__doc__ = (
        "🧪 "
        f"{firstline}"
        "\n\n.. admonition:: Experimental\n\n    This object is experimental and should not be used in production.\n"
        f"{remainder}"
    )

    @wraps(func)
    def wrapper(*args, **kwargs):
        if os.environ.get("INTERCHANGE_EXPERIMENTAL", "0") != "1":
            raise ExperimentalFeatureException(
                f"\n\tFunction or method {func.__name__} is experimental. This feature is not "
                "complete, not yet reliable, and/or needs more testing to be considered suitable "
                "for production.\n"
                "\tTo use this feature on a provisional basis, set the environment variable "
                "INTERCHANGE_EXPERIMENTAL=1.",
            )

        return func(*args, **kwargs)

    return wrapper
