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

    @property
    def root_directory(self):
        """str: Returns the directory in which all stored objects are located."""
        return self.root_directory

    def __init__(self, root_directory='stored_data'):

        self._root_directory = root_directory

        if not path.isdir(root_directory):
            os.makedirs(root_directory)

        super().__init__()

    def store_object(self, storage_key, object_to_store):

        file_path = path.join(self._root_directory, storage_key)

        try:

            with open(file_path, 'wb') as file:
                pickle.dump(object_to_store, file)

        except PicklingError:
            logging.warning('Unable to pickle an object to {}'.format(storage_key))

        super(LocalFileStorage, self).store_object(storage_key, object_to_store)

    def retrieve_object(self, storage_key):

        if not self.has_object(storage_key):
            return None

        file_path = path.join(self._root_directory, storage_key)

        loaded_object = None

        try:

            with open(file_path, 'rb') as file:
                loaded_object = pickle.load(file)

        except UnpicklingError:
            logging.warning('Unable to unpickle the object at {}'.format(storage_key))

        return loaded_object

    def has_object(self, storage_key):

        file_path = path.join(self._root_directory, storage_key)

        if not path.isfile(file_path):
            return False

        return True
