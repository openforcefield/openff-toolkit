.. _utils :

Utilities
=========

Toolkit wrappers
----------------

The toolkit wrappers provide a simple uniform API for accessing minimal functionality of cheminformatics toolkits.

These toolkit wrappers are generally used through a :class:`ToolkitRegistry`, which can be constructed with a desired precedence of toolkits:

.. code-block:: python

    >>> from openff.toolkit.utils.toolkits import ToolkitRegistry, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper
    >>> toolkit_registry = ToolkitRegistry()
    >>> toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]
    >>> [ toolkit_registry.register_toolkit(toolkit) for toolkit in toolkit_precedence if toolkit.is_available() ]
    [None, None, None]

The toolkit wrappers can then be accessed through the registry:

.. code-block:: python

    >>> from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    >>> from openff.toolkit import Molecule
    >>> molecule = Molecule.from_smiles('Cc1ccccc1')
    >>> smiles = toolkit_registry.call('to_smiles', molecule)

The order of toolkits, as specified in ``toolkit_precedence`` above, determines the order in which
the called method is resolved, i.e. if the toolkit with highest precedence has a ``to_smiles``
method, that is the toolkit that will be called. If the toolkit with highest precedence does not
have such a method, it is attempted with other toolkits until one is found. By default, if a toolkit with an appropriately-named method raises an exception of any type, then iteration over the registered toolkits stops and that exception is raised. To continue iteration if specific exceptions are encountered, customize this behavior using the optional ``raise_exception_types`` keyword argument to ``ToolkitRegistry.call``. If no registered
toolkits have the method, a ValueError is raised, containing a message listing the registered toolkits and exceptions (if any) that were ignored. 

Alternatively, the global toolkit registry (which will attempt to register any available toolkits) can be used:

.. code-block:: python

    >>> from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    >>> len(toolkit_registry.registered_toolkits)
    4

Individual toolkits can be registered or deregistered to control the backend that ToolkitRegistry calls resolve to. This can
be useful for debugging and exploring subtley different behavior between toolkit wrappers.

To temporarily change the state of ``GLOBAL_TOOLKIT_REGISTRY``, we provide the ``toolkit_registry_manager``
context manager.

.. code-block:: python

    >>> from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, AmberToolsToolkitWrapper, GLOBAL_TOOLKIT_REGISTRY
    >>> from openff.toolkit.utils import toolkit_registry_manager
    >>> print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
    4
    >>> with toolkit_registry_manager(ToolkitRegistry([RDKitToolkitWrapper(), AmberToolsToolkitWrapper()])):
    ...     print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
    2

To remove ``ToolkitWrappers`` permanently from a ``ToolkitRegistry``, the ``deregister_toolkit`` method can be used:

.. code-block:: python

    >>> from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper, BuiltInToolkitWrapper
    >>> from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    >>> print(len(toolkit_registry.registered_toolkits))
    4
    >>> toolkit_registry.deregister_toolkit(RDKitToolkitWrapper)
    >>> print(len(toolkit_registry.registered_toolkits))
    3
    >>> toolkit_registry.register_toolkit(RDKitToolkitWrapper)
    >>> print(len(toolkit_registry.registered_toolkits))
    4

For example, differences in ``to_smiles`` functionality between OpenEye toolkits and The RDKit can
be explored by selecting which toolkit(s) are and are not registered.

.. code-block:: python

    >>> from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper, GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    >>> from openff.toolkit import Molecule
    >>> molecule = Molecule.from_smiles('Cc1ccccc1')
    >>> smiles_via_openeye = toolkit_registry.call('to_smiles', molecule)
    >>> print(smiles_via_openeye)
    [H]c1c(c(c(c(c1[H])[H])C([H])([H])[H])[H])[H]

    >>> toolkit_registry.deregister_toolkit(OpenEyeToolkitWrapper)
    >>> smiles_via_rdkit = toolkit_registry.call('to_smiles', molecule)
    >>> print(smiles_via_rdkit)
    [H][c]1[c]([H])[c]([H])[c]([C]([H])([H])[H])[c]([H])[c]1[H]

.. currentmodule:: openff.toolkit.utils.toolkits
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ToolkitRegistry
    ToolkitWrapper
    OpenEyeToolkitWrapper
    RDKitToolkitWrapper
    AmberToolsToolkitWrapper
    NAGLToolkitWrapper
    BuiltInToolkitWrapper

.. currentmodule:: openff.toolkit.utils.toolkit_registry
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    toolkit_registry_manager

Serialization support
---------------------

.. currentmodule:: openff.toolkit.utils.serialization
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    Serializable

Collections
---------------------

Custom collections for the toolkit

.. currentmodule:: openff.toolkit.utils.collections
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    ValidatedList
    ValidatedDict

Miscellaneous utilities
-----------------------

Miscellaneous utility functions.

.. currentmodule:: openff.toolkit.utils.utils
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
    unit_to_string

Exceptions
----------

Exceptions raised by the Toolkit.

.. currentmodule:: openff.toolkit.utils
.. autosummary::
    :nosignatures:
    :toctree: api/generated/

    exceptions
