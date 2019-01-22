# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
File based storage API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from .base import PropertyEstimatorStorage


# =============================================================================================
# Base Backend Definition
# =============================================================================================


class MongoDBStorageBackend(PropertyEstimatorStorage):
    """A storage backend which stores files within a mongo
    database.

    .. warning :: This class has not yet been implemented.
    """

    def store_object(self, storage_key, object_to_store):
        raise NotImplementedError()

    def retrieve_object(self, storage_key):
        raise NotImplementedError()

    def has_object(self, storage_key):
        return NotImplementedError()
