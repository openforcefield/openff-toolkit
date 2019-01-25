# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
An API to define the types of messages which may be passed between a client and server.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from enum import IntEnum


# =============================================================================================
# PropertyEstimatorMessageTypes
# =============================================================================================
from pydantic import BaseModel


class PropertyEstimatorMessageTypes(IntEnum):

    Undefined = 0
    Submission = 1
    Query = 2


class PropertyEstimatorException(BaseModel):
    """A json serializable object wrapper containing information about
    a failed property calculation.

    .. todo:: Flesh out more fully.
    """
    directory: str
    message: str
