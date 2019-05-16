.. _utils ::

Utilities
=========

Toolkit wrappers
----------------

The toolkit wrappers provide a simple uniform API for accessing minimal functionality of cheminformatics toolkits.

These toolkit wrappers are generally used through a :class:`ToolkitRegistry`, which can be constructed with a desired precedence of toolkits:

.. code-block:: python

   from openforcefield.utils.toolkits import ToolkitRegistry
   toolkit_registry = ToolkitRegistry()
   toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]
   [ toolkit_registry.register(toolkit) for toolkit in toolkit_precedence if toolkit.is_available() ]

Alternatively, the global toolkit registry (which will attempt to register any available toolkits) can be used:

.. code-block:: python

   from openforcefield.utils.toolkits import DEFAULT_TOOLKIT_REGISTRY as toolkit_registry

The toolkit wrappers can then be accessed through the registry:

.. code-block:: python

   molecule = Molecule.from_smiles('Cc1ccccc1')
   smiles = toolkit_registry.call('to_smiles', molecule)

.. currentmodule:: openforcefield.utils.toolkits
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ToolkitRegistry
    ToolkitWrapper
    OpenEyeToolkitWrapper
    RDKitToolkitWrapper
    AmberToolsToolkitWrapper

Serialization support
---------------------

.. currentmodule:: openforcefield.utils.serialization
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    Serializable

Miscellaneous utilities
-----------------------

Miscellaneous utility functions.

.. currentmodule:: openforcefield.utils.utils
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    inherit_docstrings
    all_subclasses
    temporary_cd
    temporary_directory
    get_data_file_path

Structure tools
---------------

Tools for manipulating molecules and structures

.. warning ::

   These methods are deprecated and will be removed and replaced with integrated API methods.
   We recommend that no new code makes use of these functions.

.. currentmodule:: openforcefield.utils.structure
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    get_molecule_parameterIDs
