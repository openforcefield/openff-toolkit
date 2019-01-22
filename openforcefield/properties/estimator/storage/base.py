# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Base Property Estimator Storage API

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import hashlib
import os
import pickle
from os import path


# =============================================================================================
# Base Backend Definition
# =============================================================================================

class PropertyEstimatorStorage:
    """An abstract base representation of how the property estimator will interact
    with and store simulation data.

    Notes
    -----
    Any inheriting class must provide an implementation for the
    `store_object`, `retrieve_object` and `has_object` methods
    """

    def __init__(self, root_key='property_estimator'):
        """Constructs a new PropertyEstimatorStorage object.

        Parameters
        ----------
        root_key: str
            The root key of the storage system. For a file based storage system for
            example, this will be the root folder in which all files are stored
        """

        self._root_key = root_key

        self._stored_object_keys = set()
        self._stored_object_keys_file = 'object_keys'

        self._force_field_hashes = {}

        self._force_field_root = path.join(self._root_key, 'cached_force_fields')
        self._force_field_hash_file = 'hash_keys'

        if not path.isdir(self._force_field_root):
            os.makedirs(self._force_field_root)

        self._load_stored_object_keys()
        self._load_force_field_hashes()

    def _load_stored_object_keys(self):
        """Load the unique path to each object stored in the storage system.
        """
        key_path = path.join(self._root_key, self._stored_object_keys_file)

        if not self.has_object(key_path):
            return

        stored_object_keys = self.retrieve_object(key_path)

        if stored_object_keys is None:
            stored_object_keys = {}

        for unique_key in stored_object_keys:

            if not self.has_object(unique_key):
                # The force field file does not exist, so skip the entry.
                continue

            if unique_key not in self._stored_object_keys:
                self._stored_object_keys.add(unique_key)

        # Re-write a fresh copy of the file so that only force fields that
        # exist are actually referenced.
        self._save_stored_object_keys()

    def _save_stored_object_keys(self):
        """Save the unique path of each of the objects stored in the storage system.
        """
        key_path = path.join(self._root_key, self._stored_object_keys_file)
        self.store_object(key_path, self._stored_object_keys)

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
        if storage_key not in self._stored_object_keys:
            self._stored_object_keys.add(storage_key)

            self._save_stored_object_keys()

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
        raise NotImplementedError()

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
        raise NotImplementedError()

    def _load_force_field_hashes(self):
        """Load the unique id and hash keys of each of the force fields which
         have been stored in the force field directory (``self._force_field_root``).
        """
        hash_file_path = path.join(self._force_field_root,
                                   self._force_field_hash_file)

        if not self.has_object(hash_file_path):
            return

        force_field_hashes = self.retrieve_object(hash_file_path)

        if force_field_hashes is None:
            force_field_hashes = {}

        for unique_id in force_field_hashes:

            force_field_path = path.join(self._force_field_root, unique_id)

            if not self.has_object(force_field_path):
                # The force field file does not exist, so skip the entry.
                continue

            self._force_field_hashes[unique_id] = force_field_hashes[unique_id]

        # Re-write a fresh copy of the file so that only force fields that
        # exist are actually referenced.
        self._save_force_field_hashes()

    def _save_force_field_hashes(self):
        """Save the unique id and force field hash key dictionary.
        """
        hash_file_path = path.join(self._force_field_root,
                                   self._force_field_hash_file)

        self.store_object(hash_file_path, self._force_field_hashes)

    @staticmethod
    def _force_field_to_hash(force_field):
        """Converts a ForceField object to a hash
        string.

        Parameters
        ----------
        force_field: ForceField
            The force field to hash.

        Returns
        -------
        str
            The hash key of the force field.
        """
        force_field_pickle = pickle.dumps(force_field.__getstate__())
        return hashlib.sha256(force_field_pickle).hexdigest()

    def has_force_field(self, force_field):
        """Checks whether the force field has been previously
        stored in the force field directory.

        Parameters
        ----------
        force_field: ForceField
            The force field to check for.

        Returns
        -------
        str, optional
            None if the force field has not been cached, otherwise
            the unique id of the cached force field.
        """

        hash_string = self._force_field_to_hash(force_field)

        for unique_id in self._force_field_hashes:

            existing_hash = self._force_field_hashes[unique_id]

            if hash_string != existing_hash:
                continue

            existing_path = self.get_force_field_path(unique_id)

            if not self.has_object(existing_path):
                # For some reason the force field got deleted..
                continue

            return unique_id

        return None

    def retrieve_force_field(self, unique_id):
        """Retrieves a force field from storage, if it exists.

        Parameters
        ----------
        unique_id: str
            The unique id of the force field to retrieve

        Returns
        -------
        ForceField, optional
            The force field if present in the storage system with the given key, otherwise None.
        """

        force_field_path = self.get_force_field_path(unique_id)
        return self.retrieve_object(force_field_path)

    def store_force_field(self, unique_id, force_field):
        """Store the force field in the cached force field
        directory.

        .. todo::
            In future will most likely need to also serialize the
            SMIRNOFF engine version for correct comparision

        Parameters
        ----------
        unique_id: str
            The unique id assigned to the force field.
        force_field: ForceField
            The force field to cache.
        """

        hash_string = self._force_field_to_hash(force_field)
        force_field_path = self.get_force_field_path(unique_id)

        self.store_object(force_field_path, force_field)
        self._force_field_hashes[unique_id] = hash_string

        self._save_force_field_hashes()

    def get_force_field_path(self, unique_id):
        """Returns the storage path to the force field with the provided id.

        Parameters
        ----------
        unique_id: str
            The unique id of the force field whose path should be retrieved

        Returns
        -------
        str
            The path to the force field.
        """
        return path.join(self._force_field_root, unique_id)
