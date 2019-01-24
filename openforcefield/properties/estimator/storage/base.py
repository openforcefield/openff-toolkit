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
import pickle
import uuid


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

    def __init__(self):
        """Constructs a new PropertyEstimatorStorage object.
        """

        self._stored_object_keys = set()
        self._stored_object_keys_file = 'internal_object_keys'

        # Store a map between the unique id of a force field,
        # and its hash value for easy comparision of force fields.
        self._force_field_id_map = {}
        self._force_field_id_map_file = 'internal_force_field_map'

        self._simulation_data_by_substance = {}
        self._simulation_data_by_substance_file = 'internal_simulation_data_map'

        self._load_stored_object_keys()
        self._load_force_field_hashes()
        self._load_simulation_data_map()

    def _load_stored_object_keys(self):
        """Load the unique key to each object stored in the storage system.
        """
        stored_object_keys = self.retrieve_object(self._stored_object_keys_file)

        if stored_object_keys is None:
            stored_object_keys = {}

        for unique_key in stored_object_keys:

            if not self.has_object(unique_key):
                # The stored entry key does not exist in the system, so skip the entry.
                continue

            if unique_key not in self._stored_object_keys:
                self._stored_object_keys.add(unique_key)

        # Store a fresh copy of the key dictionary so that only entries
        # that exist in the system actually referenced.
        self._save_stored_object_keys()

    def _save_stored_object_keys(self):
        """Save the unique key of each of the objects stored in the storage system.
        """
        self.store_object(self._stored_object_keys_file, self._stored_object_keys)

    def store_object(self, storage_key, object_to_store):
        """Store an object in the estimators storage system.

        Parameters
        ----------
        storage_key: str
            A unique key that describes where the object will be stored in
             the storage system.
        object_to_store: Any
            The object to store. The object must be pickle serializable.
        """
        if storage_key in self._stored_object_keys:
            return

        self._stored_object_keys.add(storage_key)
        self._save_stored_object_keys()

    def retrieve_object(self, storage_key):
        """Retrieves a stored object for the estimators storage system.

        Parameters
        ----------
        storage_key: str
            A unique key that describes where the stored object can be found
            within the storage system.

        Returns
        -------
        Any, optional
            The stored object if the object key is found, otherwise None.
        """
        raise NotImplementedError()

    def has_object(self, storage_key):
        """Check whether an object with the specified key exists in the
        storage system.

        Parameters
        ----------
        storage_key: str
            A unique key that describes where the stored object can be found
            within the storage system.

        Returns
        -------
        True if the object is within the storage system.
        """
        raise NotImplementedError()

    def _load_force_field_hashes(self):
        """Load the unique id and hash keys of each of the force fields which
         have been stored in the force field directory (``self._force_field_root``).
        """
        force_field_id_map = self.retrieve_object(self._force_field_id_map_file)

        if force_field_id_map is None:
            force_field_id_map = {}

        for unique_id in force_field_id_map:

            force_field_key = 'force_field_{}'.format(unique_id)

            if not self.has_object(force_field_key):
                # The force field file does not exist, so skip the entry.
                continue

            self._force_field_id_map[unique_id] = force_field_id_map[unique_id]

        # Store a fresh copy of the hashes so that only force fields that
        # exist are actually referenced.
        self._save_force_field_hashes()

    def _save_force_field_hashes(self):
        """Save the unique id and force field hash key dictionary.
        """
        self.store_object(self._force_field_id_map_file, self._force_field_id_map)

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

        for unique_id in self._force_field_id_map:

            existing_hash = self._force_field_id_map[unique_id]

            if hash_string != existing_hash:
                continue

            force_field_key = 'force_field_{}'.format(unique_id)

            if not self.has_object(force_field_key):
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
        force_field_key = 'force_field_{}'.format(unique_id)
        return self.retrieve_object(force_field_key)

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
        force_field_key = 'force_field_{}'.format(unique_id)

        self.store_object(force_field_key, force_field)

        if unique_id not in self._force_field_id_map or hash_string != self._force_field_id_map[unique_id]:

            self._force_field_id_map[unique_id] = hash_string
            self._save_force_field_hashes()

    def _load_simulation_data_map(self):
        """Load the dictionary which tracks which stored simulation data
        was calculated for a specific substance.
        """
        _simulation_data_by_substance = self.retrieve_object(self._simulation_data_by_substance_file)

        if _simulation_data_by_substance is None:
            _simulation_data_by_substance = {}

        for substance_id in _simulation_data_by_substance:

            self._simulation_data_by_substance[substance_id] = []

            for simulation_data_key in _simulation_data_by_substance[substance_id]:

                if not self.has_object(simulation_data_key):
                    # The force field file does not exist, so skip the entry.
                    continue

                self._simulation_data_by_substance[substance_id].append(simulation_data_key)

        # Store a fresh copy of the hashes so that only force fields that
        # exist are actually referenced.
        self._save_simulation_data_map()

    def _save_simulation_data_map(self):
        """Save the unique id and simulation data key by substance dictionary.
        """
        self.store_object(self._simulation_data_by_substance_file, self._simulation_data_by_substance)

    def retrieve_simulation_data(self, substance_id):
        """Retrieves any data that has been stored
        for a given substance.

        Parameters
        ----------
        substance_id: str
            The id of the substance to check for.

        Returns
        -------
        list of StoredSimulationData, optional
            A list of the stored data if present in the storage system, otherwise None.
        """

        if substance_id not in self._simulation_data_by_substance:
            return None

        returned_data = []

        for simulation_data_key in self._simulation_data_by_substance[substance_id]:
            returned_data.append(self.retrieve_object(simulation_data_key))

        return returned_data

    def store_simulation_data(self, substance_id, simulation_data):
        """Store the simulation data.

        Notes
        -----
        If the storage system already contains equivalent information (i.e data stored
        for the same substance, thermodynamic state and parameter set) then the
        data with the longest autocorrelation time will be retained.

        Parameters
        ----------
        substance_id: str
            The id of the substance to which the data belongs.
        simulation_data: StoredSimulationData
            The simulation data to store.
        """

        simulation_data_key = None
        data_to_store = None

        if substance_id in self._simulation_data_by_substance:

            for stored_data_key in self._simulation_data_by_substance[substance_id]:

                stored_data = self.retrieve_object(stored_data_key)

                if stored_data is None:
                    continue

                if simulation_data.thermodynamic_state != stored_data.thermodynamic_state:
                    continue

                if simulation_data.parameter_set_id != stored_data.parameter_set_id:
                    continue

                if stored_data.autocorrelation_time < simulation_data.autocorrelation_time:
                    continue

                if (simulation_data.autocorrelation_time == stored_data.autocorrelation_time and
                    stored_data.effective_samples < simulation_data.effective_samples):
                    continue

                data_to_store = stored_data
                simulation_data_key = stored_data_key

        if simulation_data_key is None:

            simulation_data_key = "{}_{}".format(substance_id, uuid.uuid4())
            data_to_store = simulation_data

        self.store_object(simulation_data_key, data_to_store)

        if (substance_id not in self._simulation_data_by_substance or
                simulation_data_key not in self._simulation_data_by_substance[substance_id]):

            if substance_id not in self._simulation_data_by_substance:
                self._simulation_data_by_substance[substance_id] = []

            self._simulation_data_by_substance[substance_id].append(simulation_data_key)
            self._save_simulation_data_map()
