.. _utils ::

Utilities
=========

Toolkit wrappers
----------------

The toolkit wrappers provide a simple uniform API for accessing minimal functionality of cheminformatics toolkits.

These toolkit wrappers are generally used through a :class:`ToolkitRegistry`, which can be constructed with a desired precedence of toolkits:

.. code-block:: python

    from openforcefield.utils.toolkits import ToolkitRegistry, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper
    toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]
    [toolkit.is_available() for toolkit in toolkit_precedence]
    toolkit_registry = ToolkitRegistry(toolkit_precedence=toolkit_precedence)
    toolkit_registry.registered_toolkits

The toolkit wrappers can then be accessed through the registry:

.. code-block:: python

    molecule = Molecule.from_smiles('Cc1ccccc1')
    smiles = toolkit_registry.call('to_smiles', molecule)

The order of toolkits, as specified in ``toolkit_precedence`` above, determines the order in which
the called method is resolved, i.e. if the toolkit with highest precedence has a ``to_smiles``
method, that is the toolkit that will be called. If the toolkit with highest precedence does not
have such a method, it is attempted with other toolkits until one iss found. If no registered
toolkits have the method, an informative error message is presented.

Alternatively, the global toolkit registry (which will attempt to register any available toolkits) can be used:

.. code-block:: python

    from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    toolkit_registry.registered_toolkits

Individual toolkits can be registered or deregistered to control which toolkits are called. This can
be useful for debugging and exploring subtley different behavior between toolkit wrappers.

.. code-block:: python

    from openforcefield.utils.toolkits import OpenEyeToolkitWrapper, BuiltInToolkitWrapper
    from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    toolkit_registry.deregister_toolkit(OpenEyeToolkitWrapper)
    toolkit_registry.register_toolkit(BuiltInToolkitWrapper)
    toolkit_registry.registered_toolkits

For example, differences in ``to_smiles`` functionality between OpenEye toolkits and The RDKit can
be explored by selecting which toolkit(s) are and are not registered.

.. code-block::python

    from openforcefield.utils.toolkits import ToolkitRegistry, OpenEyeToolkitWrapper, RDKitToolkitWrapper
    toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper, RDKitToolkitWrapper])
    smiles_via_openeye = toolkit_registry.call('to_smiles', molecule)
    toolkit_registry.deregister_toolkit(OpenEyeToolkitWrapper)
    smiles_via_rdkit = toolkit_registry.call('to_smiles', molecule)


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
    get_data_file_path
    convert_0_1_smirnoff_to_0_2
    convert_0_2_smirnoff_to_0_3
    get_molecule_parameterIDs
