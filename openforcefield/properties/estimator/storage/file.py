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

import logging
import os
import pickle
from os import path
from pickle import PicklingError, UnpicklingError

from .base import PropertyEstimatorStorage


# =============================================================================================
# Base Backend Definition
# =============================================================================================


class LocalFileStorage(PropertyEstimatorStorage):
    """A storage backend which stores files normally on the local
    disk.
    """

    def __init__(self, root_key='property_estimator'):
        super().__init__(root_key)

        if not path.isdir(root_key):
            os.makedirs(root_key)

    def store_object(self, storage_key, object_to_store):
        """Store an object in the estimators storage system.

        Parameters
        ----------
        storage_key: str
            A unique key that describes where (/ how) in the storage system an
            object gets stored. This may be a file path for example.
        object_to_store: Any
            The object to store. The object must be pickle serializable.
        """

        file_path = storage_key

        if file_path.find(self._root_key) < 0:
            file_path = path.join(self._root_key, file_path)

        try:

            with open(file_path, 'wb') as file:
                pickle.dump(object_to_store, file)

        except PicklingError:
            logging.warning('Unable to pickle an object to {}'.format(storage_key))

        super(LocalFileStorage, self).store_object(file_path, object_to_store)

    def retrieve_object(self, storage_key):
        """Retrieves a stored object for the estimators storage system.

        Parameters
        ----------
        storage_key: str
            A unique key that describes where (/ how) in the storage system an
            object gets stored. This may be a file path for example.

        Returns
        -------
        Any, None
            The stored object if present, otherwise None.
        """
        file_path = storage_key

        if file_path.find(self._root_key) < 0:
            file_path = path.join(self._root_key, file_path)

        if not self.has_object(file_path):
            return None

        loaded_object = None

        try:

            with open(file_path, 'rb') as file:
                loaded_object = pickle.load(file)

        except UnpicklingError:
            logging.warning('Unable to unpickle the object at {}'.format(storage_key))

        return loaded_object

    def has_object(self, storage_key):
        """Check whether an object with the specified key exists in the
        storage system.

        Parameters
        ----------
        storage_key: str
            A unique key that describes where (/ how) in the storage system an
            object gets stored. This may be a file path for example.

        Returns
        -------
        True if the object is within the storage system.
        """

        file_path = storage_key

        if file_path.find(self._root_key) < 0:
            file_path = path.join(self._root_key, file_path)

        if not path.isfile(file_path):
            return False

        return True
