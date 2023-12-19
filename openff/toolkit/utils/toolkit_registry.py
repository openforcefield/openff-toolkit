"Registry for ToolkitWrapper objects"

__all__ = ("ToolkitRegistry", "toolkit_registry_manager")

import inspect
import logging
from contextlib import contextmanager
from typing import Callable, Optional, Union

from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
from openff.toolkit.utils.base_wrapper import ToolkitWrapper
from openff.toolkit.utils.builtin_wrapper import BuiltInToolkitWrapper
from openff.toolkit.utils.exceptions import (
    InvalidToolkitError,
    LicenseError,
    ToolkitUnavailableException,
)
from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper
from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
from openff.toolkit.utils.utils import all_subclasses

logger = logging.getLogger(__name__)


class ToolkitRegistry:
    """
    Registry for ToolkitWrapper objects

    Examples
    --------

    Register toolkits in a specified order, skipping if unavailable

    >>> from openff.toolkit.utils.toolkits import ToolkitRegistry
    >>> toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]
    >>> toolkit_registry = ToolkitRegistry(toolkit_precedence)
    >>> toolkit_registry
    <ToolkitRegistry containing OpenEye Toolkit, The RDKit, AmberTools>

    Register all available toolkits (in the order OpenEye, RDKit, AmberTools, built-in)

    >>> toolkits = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, BuiltInToolkitWrapper]
    >>> toolkit_registry = ToolkitRegistry(toolkit_precedence=toolkits)
    >>> toolkit_registry
    <ToolkitRegistry containing OpenEye Toolkit, The RDKit, AmberTools, Built-in Toolkit>

    Retrieve the global singleton toolkit registry, which is created when this module is imported from all available
    toolkits:

    >>> from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY as toolkit_registry
    >>> toolkit_registry
    <ToolkitRegistry containing OpenEye Toolkit, The RDKit, AmberTools, Built-in Toolkit>

    Note that this will contain different ToolkitWrapper objects based on what toolkits
    are currently installed.

    .. warning :: This API is experimental and subject to change.
    """

    _toolkits: list[ToolkitWrapper]

    def __init__(
        self,
        toolkit_precedence: Optional[list[type[ToolkitWrapper]]] = None,
        exception_if_unavailable: bool = True,
        _register_imported_toolkit_wrappers: bool = False,
    ):
        """
        Create an empty toolkit registry.

        Parameters
        ----------
        toolkit_precedence
            List of toolkit wrapper classes, in order of desired precedence when performing molecule operations. If
            None, no toolkits will be registered.

        exception_if_unavailable
            If True, an exception will be raised if the toolkit is unavailable

        _register_imported_toolkit_wrappers
            If True, will attempt to register all imported ToolkitWrapper subclasses that can be
            found in the order of toolkit_precedence, if specified. If toolkit_precedence is not
            specified, the default order is [OpenEyeToolkitWrapper, RDKitToolkitWrapper,
            AmberToolsToolkitWrapper, BuiltInToolkitWrapper].

        """
        self._toolkits = list()

        toolkit_classes_to_register: list[type[ToolkitWrapper]] = list()

        if _register_imported_toolkit_wrappers:
            if toolkit_precedence is None:
                toolkit_precedence = [
                    OpenEyeToolkitWrapper,
                    RDKitToolkitWrapper,
                    AmberToolsToolkitWrapper,
                    BuiltInToolkitWrapper,
                ]
            all_importable_toolkit_wrappers = all_subclasses(ToolkitWrapper)
            for toolkit in toolkit_precedence:
                if toolkit in all_importable_toolkit_wrappers:
                    toolkit_classes_to_register.append(toolkit)
        else:
            if toolkit_precedence is not None:
                toolkit_classes_to_register = toolkit_precedence

        if toolkit_classes_to_register:
            for toolkit in toolkit_classes_to_register:
                self.register_toolkit(
                    toolkit_wrapper=toolkit,
                    exception_if_unavailable=exception_if_unavailable,
                )

    @property
    def registered_toolkits(self) -> list[ToolkitWrapper]:
        """
        List registered toolkits.

        .. warning :: This API is experimental and subject to change.

        .. todo :: Should this return a generator? Deep copies? Classes? Toolkit names?

        Returns
        -------
        toolkits
        """
        return list(self._toolkits)

    @property
    def registered_toolkit_versions(self) -> dict[str, str]:
        """
        Return a dict containing the version of each registered toolkit.

        .. warning :: This API is experimental and subject to change.

        Returns
        -------
        toolkit_versions
            A dictionary mapping names and versions of wrapped toolkits

        """
        return {tk.toolkit_name: tk.toolkit_version for tk in self.registered_toolkits}

    def register_toolkit(
        self,
        toolkit_wrapper: Union[ToolkitWrapper, type[ToolkitWrapper]],
        exception_if_unavailable: bool = True,
    ):
        """
        Register the provided toolkit wrapper class, instantiating an object of it.

        .. warning :: This API is experimental and subject to change.

        .. todo ::

           This method should raise an exception if the toolkit is unavailable, unless an optional argument
           is specified that silently avoids registration of toolkits that are unavailable.

        Parameters
        ----------
        toolkit_wrapper
            The toolkit wrapper to register or its class.
        exception_if_unavailable
            If True, an exception will be raised if the toolkit is unavailable

        """
        # Instantiate class if class, or just add if already instantiated.
        if isinstance(toolkit_wrapper, ToolkitWrapper):
            self._toolkits.append(toolkit_wrapper)

        elif issubclass(toolkit_wrapper, ToolkitWrapper):
            try:
                _toolkit_wrapper = toolkit_wrapper()

            # This exception can be raised by OpenEyeToolkitWrapper
            except LicenseError as license_exception:
                if exception_if_unavailable:
                    raise ToolkitUnavailableException(license_exception.msg)
                else:
                    logger.warning(license_exception)
                return
            except ToolkitUnavailableException:
                if exception_if_unavailable:
                    raise ToolkitUnavailableException(
                        "Unable to load toolkit '{_toolkit_wrapper._toolkit_name}'. "
                    )
                return

            self._toolkits.append(_toolkit_wrapper)

        else:
            raise ValueError(f"Given unexpected argument {type(toolkit_wrapper)=}")

    def deregister_toolkit(self, toolkit_wrapper: ToolkitWrapper):
        """
        Remove a ToolkitWrapper from the list of toolkits in this ToolkitRegistry

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        toolkit_wrapper
            The toolkit wrapper to remove from the registry

        Raises
        ------
        InvalidToolkitError
            If toolkit_wrapper is not a ToolkitWrapper or subclass
        ToolkitUnavailableException
            If toolkit_wrapper is not found in the registry
        """
        # If passed a class, instantiate it
        if inspect.isclass(toolkit_wrapper):
            toolkit_wrapper = toolkit_wrapper()

        if not isinstance(toolkit_wrapper, ToolkitWrapper):
            msg = (
                f"Argument {toolkit_wrapper} must an ToolkitWrapper "
                f"or subclass of it. Found type {type(toolkit_wrapper)}."
            )
            raise InvalidToolkitError(msg)

        toolkits_to_remove = []

        for toolkit in self._toolkits:
            if type(toolkit) is type(toolkit_wrapper):
                toolkits_to_remove.append(toolkit)

        if not toolkits_to_remove:
            msg = (
                f"Did not find {toolkit_wrapper} in registry. "
                f"Currently registered toolkits are {self._toolkits}"
            )
            raise ToolkitUnavailableException(msg)

        for toolkit_to_remove in toolkits_to_remove:
            self._toolkits.remove(toolkit_to_remove)

    def add_toolkit(self, toolkit_wrapper: ToolkitWrapper):
        """
        Append a ToolkitWrapper onto the list of toolkits in this ToolkitRegistry

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        toolkit_wrapper
            The ToolkitWrapper object to add to the list of registered toolkits

        Raises
        ------
        InvalidToolkitError
            If toolkit_wrapper is not a ToolkitWrapper or subclass
        """
        if not isinstance(toolkit_wrapper, ToolkitWrapper):
            raise InvalidToolkitError(
                "Something other than a ToolkitWrapper object was passed to ToolkitRegistry.add_toolkit()\n"
                f"Given object {toolkit_wrapper} of type {type(toolkit_wrapper)}"
            )

        self._toolkits.append(toolkit_wrapper)

    # TODO: Can we automatically resolve calls to methods that are not explicitly defined using some Python magic?

    def resolve(self, method_name: str) -> Callable:
        """
        Resolve the requested method name by checking all registered toolkits in
        order of precedence for one that provides the requested method.

        Parameters
        ----------
        method_name
            The name of the method to resolve

        Returns
        -------
        method
            The method of the first registered toolkit that provides the requested method name

        Raises
        ------
        NotImplementedError if the requested method cannot be found among the registered toolkits

        Examples
        --------

        Create a molecule, and call the toolkit ``to_smiles()`` method directly

        >>> from openff.toolkit import Molecule
        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry([OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper])
        >>> method = toolkit_registry.resolve('to_smiles')
        >>> smiles = method(molecule)

        .. todo :: Is there a better way to figure out which toolkits implement given methods by introspection?

        """
        for toolkit in self._toolkits:
            if hasattr(toolkit, method_name):
                method = getattr(toolkit, method_name)
                return method

        # No toolkit was found to provide the requested capability
        # TODO: Can we help developers by providing a check for typos in expected method names?
        raise NotImplementedError(
            f'No registered toolkits can provide the capability "{method_name}".\n'
            f"Available toolkits are: {self.registered_toolkits}\n"
        )

    # TODO: Can we instead register available methods directly with `ToolkitRegistry`,
    # so we can just use `ToolkitRegistry.method()`?
    def call(
        self,
        method_name: str,
        *args,
        raise_exception_types: Optional[list[type[Exception]]] = None,
        **kwargs,
    ):
        """
        Execute the requested method by attempting to use all registered toolkits in order of precedence.

        ``*args`` and ``**kwargs`` are passed to the desired method, and return values of the method are returned

        This is a convenient shorthand for ``toolkit_registry.resolve_method(method_name)(*args, **kwargs)``

        Parameters
        ----------
        method_name
            The name of the method to execute
        raise_exception_types
            A list of exception-derived types to catch and raise immediately. If None, this will be set to [Exception],
            which will raise an error immediately if the first ToolkitWrapper in the registry fails. To try each
            ToolkitWrapper that provides a suitably-named method, set this to the empty list ([]). If all
            ToolkitWrappers run without raising any exceptions in this list, a single ValueError will be raised
            containing the each ToolkitWrapper that was tried and the exception it raised.

        Raises
        ------
        NotImplementedError if the requested method cannot be found among the registered toolkits

        ValueError if no exceptions in the raise_exception_types list were raised by ToolkitWrappers, and
        all ToolkitWrappers in the ToolkitRegistry were tried.

        Other forms of exceptions are possible if raise_exception_types is specified.
        These are defined by the ToolkitWrapper method being called.

        Examples
        --------

        Create a molecule, and call the toolkit ``to_smiles()`` method directly

        >>> from openff.toolkit import Molecule
        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry([OpenEyeToolkitWrapper, RDKitToolkitWrapper])
        >>> smiles = toolkit_registry.call('to_smiles', molecule)

        """
        if raise_exception_types is None:
            raise_exception_types = [Exception]

        errors = list()
        for toolkit in self._toolkits:
            if hasattr(toolkit, method_name):
                method = getattr(toolkit, method_name)
                try:
                    return method(*args, **kwargs)
                except Exception as e:
                    for exception_type in raise_exception_types:
                        if isinstance(e, exception_type):
                            raise e
                    errors.append((toolkit, e))

        # No toolkit was found to provide the requested capability
        # TODO: Can we help developers by providing a check for typos in expected method names?
        msg = (
            f'No registered toolkits can provide the capability "{method_name}" '
            f'for args "{args}" and kwargs "{kwargs}"\n'
        )

        msg += "Available toolkits are: {}\n".format(self.registered_toolkits)
        # Append information about toolkits that implemented the method, but could not handle the provided parameters
        for toolkit, error in errors:
            msg += f" {toolkit} {type(error)} : {error}\n"
        raise ValueError(msg)

    def __repr__(self):
        return f"<ToolkitRegistry containing {', '.join([tk.toolkit_name for tk in self._toolkits])}>"


# Copied from https://github.com/openforcefield/openff-fragmenter/blob/4a290b866a8ed43eabcbd3231c62b01f0c6d7df6
# /openff/fragmenter/utils.py#L97-L123
@contextmanager
def toolkit_registry_manager(toolkit_registry: Union[ToolkitRegistry, ToolkitWrapper]):
    """
    A context manager that temporarily changes the ToolkitWrappers in the GLOBAL_TOOLKIT_REGISTRY.
    This can be useful in cases where one would otherwise need to otherwise manually specify the use of a specific
    ToolkitWrapper or ToolkitRegistry repeatedly in a block of code, or in cases where there isn't another way to
    switch the ToolkitWrappers used for a particular operation.

    Examples
    --------
    >>> from openff.toolkit import Molecule, RDKitToolkitWrapper, AmberToolsToolkitWrapper
    >>> from openff.toolkit.utils import toolkit_registry_manager, ToolkitRegistry
    >>> mol = Molecule.from_smiles("CCO")
    >>> print(mol.to_smiles()) # This will use the OpenEye backend (if installed and licensed)
    [H]C([H])([H])C([H])([H])O[H]
    >>> with toolkit_registry_manager(ToolkitRegistry([RDKitToolkitWrapper()])):
    ...    print(mol.to_smiles())
    [H][O][C]([H])([H])[C]([H])([H])[H]

    """
    from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY

    if isinstance(toolkit_registry, ToolkitRegistry):
        context_toolkits = toolkit_registry.registered_toolkits
    elif isinstance(toolkit_registry, ToolkitWrapper):
        context_toolkits = [toolkit_registry]
    else:
        raise NotImplementedError(
            "Only ``ToolkitRegistry`` and ``ToolkitWrapper`` are supported."
        )

    original_toolkits = GLOBAL_TOOLKIT_REGISTRY.registered_toolkits

    for toolkit in original_toolkits:
        GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkit)

    for toolkit in context_toolkits:
        GLOBAL_TOOLKIT_REGISTRY.register_toolkit(toolkit)

    try:
        yield
    finally:
        for toolkit in context_toolkits:
            GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkit)
        for toolkit in original_toolkits:
            print(len(original_toolkits), toolkit)
            GLOBAL_TOOLKIT_REGISTRY.register_toolkit(toolkit)


_toolkit_registry_manager = toolkit_registry_manager
