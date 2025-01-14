# Toolkit Wrappers and Registries

The toolkit wrappers provide a simple uniform API for accessing some functionality of cheminformatics toolkits. They are used internally by the OpenFF Toolkit to avoid re-implementing existing scientific algorithms, and can be used by users to specify exactly what software is used for any given section of code.

These toolkit wrappers are generally used through a [`ToolkitRegistry`], which facilitates combining software with different capabilities and can be constructed with a list of toolkit wrappers ordered by precedence. 

```pycon
>>> from openff.toolkit import (
...     ToolkitRegistry,
...     OpenEyeToolkitWrapper,
...     RDKitToolkitWrapper,
...     AmberToolsToolkitWrapper,
... )
>>> toolkit_registry = ToolkitRegistry(
...     [
...         OpenEyeToolkitWrapper,
...         RDKitToolkitWrapper,
...         AmberToolsToolkitWrapper,
...     ]
... )

```

The toolkit wrappers' functionality can then be accessed through the registry. The first toolkit in the list that provides a method with the given name will be used:

```pycon
>>> from openff.toolkit import Molecule
>>> molecule = Molecule.from_smiles("Cc1ccccc1")
>>> smiles = toolkit_registry.call("to_smiles", molecule)

```

For further details on how this search is performed and how it handles exceptions, see the [`ToolkitRegistry.call()` API docs].

Many functions in the OpenFF Toolkit API include a `toolkit_registry` argument that can be used to specify the toolkit wrappers used by a call to that function. The value of this argument can be either a single toolkit wrapper instance, or an entire toolkit registry:

```pycon
>>> smiles = molecule.to_smiles(toolkit_registry=RDKitToolkitWrapper())
>>> smiles = molecule.to_smiles(toolkit_registry=toolkit_registry)

```

For example, differences in `to_smiles` functionality between the OpenEye and RDKit can be explored by specifying the desired toolkit wrapper:

```pycon
>>> molecule.to_smiles(toolkit_registry=RDKitToolkitWrapper())
'[H][c]1[c]([H])[c]([H])[c]([C]([H])([H])[H])[c]([H])[c]1[H]'
>>> molecule.to_smiles(toolkit_registry=OpenEyeToolkitWrapper())
'[H]c1c(c(c(c(c1[H])[H])C([H])([H])[H])[H])[H]'

```

The default value of this argument is the [`GLOBAL_TOOLKIT_REGISTRY`], which by default includes all the toolkits that OpenFF recommends for everyday use:

```pycon
>>> from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY
>>> len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits)
4

```

The [`toolkit_registry_manager`] context manager allows `GLOBAL_TOOLKIT_REGISTRY` to be changed temporarily:

```pycon
>>> from openff.toolkit.utils import toolkit_registry_manager
>>> print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
4
>>> with toolkit_registry_manager(
...     ToolkitRegistry([RDKitToolkitWrapper(), AmberToolsToolkitWrapper()])
... ):
...     print(len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits))
...
2

```

For more information about modifying `GLOBAL_TOOLKIT_REGISTRY`, see the [`GLOBAL_TOOLKIT_REGISTRY` API docs].

[`ToolkitRegistry`]: openff.toolkit.utils.toolkits.ToolkitRegistry
[`ToolkitRegistry.call()` API docs]: openff.toolkit.utils.toolkits.ToolkitRegistry.call
[`ToolkitRegistry.call`]: openff.toolkit.utils.toolkits.ToolkitRegistry.call
[`GLOBAL_TOOLKIT_REGISTRY`]: openff.toolkit.utils.toolkits.GLOBAL_TOOLKIT_REGISTRY
[`toolkit_registry_manager`]: openff.toolkit.utils.toolkits.toolkit_registry_manager
[`GLOBAL_TOOLKIT_REGISTRY` API docs]: openff.toolkit.utils.toolkits.GLOBAL_TOOLKIT_REGISTRY

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