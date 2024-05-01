# Toolkit Wrappers and Registries

The toolkit wrappers provide a simple uniform API for accessing minimal functionality of cheminformatics toolkits.

These toolkit wrappers are generally used through a [`ToolkitRegistry`], which can be constructed with a list of toolkits ordered by precedence. 

```python
>>> from openff.toolkit import ToolkitRegistry, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper
>>> toolkit_registry = ToolkitRegistry([
...     OpenEyeToolkitWrapper, 
...     RDKitToolkitWrapper, 
...     AmberToolsToolkitWrapper,
... ])

```

The toolkit wrapper's functionality can then be accessed through the registry. The first toolkit in the list that provides a method with the given name will be used:

```python
>>> from openff.toolkit import Molecule
>>> molecule = Molecule.from_smiles('Cc1ccccc1')
>>> smiles = toolkit_registry.call('to_smiles', molecule)

```

For further details on how this search is performed and how it handles exceptions, see the [`ToolkitRegistry.call()` API docs].

Many functions in the OpenFF Toolkit API include a `toolkit_registry` argument that can be used to specify the toolkit wrappers used by a call to that function. The value of this argument can be either a single toolkit wrapper instance, or an entire toolkit registry:

```python
>>> smiles = molecule.to_smiles(toolkit_registry=RDKitToolkitWrapper())
>>> smiles = molecule.to_smiles(toolkit_registry=toolkit_registry)

```

For example, differences in `to_smiles` functionality between OpenEye toolkits and The RDKit can be explored by specifying the desired toolkit wrapper:

```python
>>> molecule.to_smiles(toolkit_registry=RDKitToolkitWrapper())
[H]c1c(c(c(c(c1[H])[H])C([H])([H])[H])[H])[H]
>>> molecule.to_smiles(toolkit_registry=OpenEyeToolkitWrapper())
[H][c]1[c]([H])[c]([H])[c]([C]([H])([H])[H])[c]([H])[c]1[H]

```

The default value of this argument is the [`GLOBAL_TOOLKIT_REGISTRY`], which by default includes all the toolkits that OpenFF recommends for everyday use:

```python
>>> from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY
>>> len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits)
4

```

To temporarily change the state of `GLOBAL_TOOLKIT_REGISTRY`, we provide the `toolkit_registry_manager` context manager.

```python
>>> from openff.toolkit.utils import toolkit_registry_manager
>>> print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
4
>>> with toolkit_registry_manager(ToolkitRegistry([RDKitToolkitWrapper(), AmberToolsToolkitWrapper()])):
...     print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
2

```

To change the value of `GLOBAL_TOOLKIT_REGISTRY` (or any registry) for the remainder of the Python process, individual toolkits can be registered or deregistered. This can be useful for debugging and exploring subtly different behavior between toolkit wrappers. To remove `ToolkitWrappers` from a `ToolkitRegistry`, the `deregister_toolkit` method can be used:

```python
>>> print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
4
>>> GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(RDKitToolkitWrapper)
>>> print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
3
>>> GLOBAL_TOOLKIT_REGISTRY.register_toolkit(RDKitToolkitWrapper)
>>> print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
4

```

Note that as with other global attributes in Python, assigning a new toolkit registry to `GLOBAL_TOOLKIT_REGISTRY` is quite difficult to get right and very likely to fail silently - we recommend modifying the existing value instead in most cases.

[`ToolkitRegistry`]: openff.toolkit.utils.toolkits.ToolkitRegistry
[`ToolkitRegistry.call` API docs]: openff.toolkit.utils.toolkits.ToolkitRegistry.call
[`ToolkitRegistry.call`]: openff.toolkit.utils.toolkits.ToolkitRegistry.call
[`GLOBAL_TOOLKIT_REGISTRY`]: openff.toolkit.utils.toolkits.GLOBAL_TOOLKIT_REGISTRY

```{eval-rst}
.. currentmodule:: openff.toolkit.utils.toolkits
.. autosummary::
    :nosignatures:
    :toctree: generated/

    ToolkitRegistry
    ToolkitWrapper
    OpenEyeToolkitWrapper
    RDKitToolkitWrapper
    AmberToolsToolkitWrapper
    NAGLToolkitWrapper
    BuiltInToolkitWrapper
    GLOBAL_TOOLKIT_REGISTRY
    toolkit_registry_manager
```