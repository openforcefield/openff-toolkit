.. _utils ::

Utilities
=========

Toolkit wrappers
----------------

Tools for representing and operating on chemical environments

.. currentmodule:: openforcefield.utils.toolkits
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ToolkitRegistry
    ToolkitWrapper
    OpenEyeToolkitWrapper
    RDKitToolkitWrapper
    AmberToolsWrapper

Serialization support
---------------------

.. currentmodule:: openforcefield.utils.serialization
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    Serializable

Structure tools
---------------

Tools for manipulating molecules and structures

.. todo :: These methods are deprecated and will be removed.

.. currentmodule:: openforcefield.utils.structure
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    generateSMIRNOFFStructure
    generateProteinStructure
    combinePositions
    mergeStructure
    generateTopologyFromOEMol
    get_testdata_filename
    normalize_molecules
    read_molecules
    setPositionsInOEMol
    extractPositionsFromOEMol
    read_typelist
    positions_from_oemol
    check_energy_is_finite
    get_energy
    get_molecule_parameterIDs
    getMolParamIDToAtomIndex
    merge_system
    save_system_to_amber
    save_system_to_gromacs

Miscellaneous utilities
-----------------------

.. currentmodule:: openforcefield.utils.utils
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    all_subclasses
    temporary_cd
    temporary_directory
    get_data_filename
