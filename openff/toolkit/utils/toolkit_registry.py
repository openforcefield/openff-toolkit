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
    A registry and precedence list for :py:class:`ToolkitWrapper` objects.

    ``ToolkitRegistry`` allows the OpenFF Toolkit to provide a concise,
    universal API that can call out to other well-established libraries rather
    than re-implement algorithms with well-established implementations, while
    still giving users control over dependencies. It contains a list of
    :py:class:`ToolkitWrapper` objects, each of which provides some collection
    of methods that may be requested from the registry. The :py:meth:`call`
    method takes a name and calls the method of that name on each wrapper until
    it finds a working implementation, whose result it returns. For details on
    how this search is conducted and what counts as a working implementation,
    see that method's API docs.

    Examples
    --------

    Register toolkits in a specified order, skipping if unavailable

    >>> from openff.toolkit.utils.toolkits import ToolkitRegistry
    >>> toolkit_precedence = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]
    >>> toolkit_registry = ToolkitRegistry(toolkit_precedence)
    >>> toolkit_registry
    <ToolkitRegistry containing OpenEye Toolkit, The RDKit, AmberTools>

    Register all available toolkits (in the order OpenEye, RDKit, AmberTools,
    built-in)

    >>> toolkits = [OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, BuiltInToolkitWrapper]
    >>> toolkit_registry = ToolkitRegistry(toolkit_precedence=toolkits)
    >>> toolkit_registry
    <ToolkitRegistry containing OpenEye Toolkit, The RDKit, AmberTools, Built-in Toolkit>

    Retrieve the global singleton toolkit registry, which is created when this
    module is imported from all available toolkits:

    >>> from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
    >>> GLOBAL_TOOLKIT_REGISTRY
    <ToolkitRegistry containing OpenEye Toolkit, The RDKit, AmberTools, Built-in Toolkit>

    Note that this will contain different ToolkitWrapper objects based on what
    toolkits are currently installed.

    Call a method on the registered toolkit wrapper with the highest precedence:

    >>> molecule = toolkit_registry.call('from_smiles', 'Cc1ccccc1')

    For more, see the :py:meth:`call` method.

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
        Execute a method with the first registered toolkits that supports it.

        This method searches the registry's precedence list for the first
        wrapper that has an attribute named ``method_name`` and attempts to
        call it as a method using the arguments in ``*args`` and ``**kwargs``.
        If that method raises no exception, its return value is returned.

        By default, if a wrapper with an appropriately-named method raises an
        exception of any type, then iteration over the registered toolkits
        stops early and that exception is raised. To limit this behavior to only
        certain exceptions and otherwise continue iteration, customize this
        behavior using the optional ``raise_exception_types`` keyword argument.
        If iteration finishes without finding a wrapper that can successfully
        call the requested method, a ``ValueError`` is raised, containing a
        message listing the registered toolkits and any exceptions that
        occurred.

        Parameters
        ----------
        method_name
            The name of the method to execute
        raise_exception_types
            A list of exception-derived types to catch and raise immediately. If
            ``None``, the first exception encountered will be raised. To ignore
            all exceptions, set this to the empty list ``[]``.

        Raises
        ------
        ValueError
            If no suitable toolkit wrapper was found in the registry.

        Exception
            Other forms of exceptions are possible if raise_exception_types is
            specified. These are defined by the ``ToolkitWrapper`` method being
            called.

        Examples
        --------

        Create a molecule, and call the toolkit ``to_smiles()`` method directly

        >>> from openff.toolkit import (
        ...     Molecule,
        ...     ToolkitRegistry,
        ...     AmberToolsToolkitWrapper,
        ...     RDKitToolkitWrapper,
        ... )
        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry([
        ...     AmberToolsToolkitWrapper,
        ...     RDKitToolkitWrapper,
        ... ])
        >>> smiles = toolkit_registry.call('to_smiles', molecule)

        Stop if a partial charge assignment method encounters an error during
        the partial charge calculation (:py:exc:`ChargeCalculationError`), but
        proceed to the next wrapper for any other exception such as the wrapper
        not supporting the charge method:

        >>> from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY
        >>> from openff.toolkit.utils.exceptions import ChargeCalculationError
        >>> GLOBAL_TOOLKIT_REGISTRY.call(
        ...     "assign_partial_charges",
        ...     molecule=Molecule.from_smiles('C'),
        ...     partial_charge_method="gasteiger",
        ...     raise_exception_types=[ChargeCalculationError],
        ... )

        Calling a method that exists on none of the wrappers raises a
        ``ValueError``:

        >>> from openff.toolkit import ToolkitRegistry, RDKitToolkitWrapper
        >>> toolkit_registry = ToolkitRegistry([RDKitToolkitWrapper])
        >>> toolkit_registry.call("there_is_no_spoon") # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        Traceback (most recent call last):
          ...
        ValueError: No registered toolkits can provide the capability
        "there_is_no_spoon" for args "()" and kwargs "{}"
        Available toolkits are: [ToolkitWrapper around The RDKit version ...]
        <BLANKLINE>

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

        msg += f"Available toolkits are: {self.registered_toolkits}\n"
        # Append information about toolkits that implemented the method, but could not handle the provided parameters
        for toolkit, error in errors:
            msg += f" {toolkit} {type(error)} : {error}\n"
        raise ValueError(msg)

    def __repr__(self):
        return f"<ToolkitRegistry containing {', '.join([tk.toolkit_name for tk in self._toolkits])}>"

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
                    raise ToolkitUnavailableException(f"Unable to load toolkit '{toolkit_wrapper}'. ")
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
            msg = f"Did not find {toolkit_wrapper} in registry. Currently registered toolkits are {self._toolkits}"
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
        Get a method with the given name from the first registered toolkit that provides it.

        Resolve the requested method name by checking all registered toolkits in
        order of precedence for one that provides the requested method. Note
        that this may not be the method used by :py:meth:`call` if the
        ``raise_exception_types`` argument is passed.


        Parameters
        ----------
        method_name
            The name of the method to resolve

        Returns
        -------
        method
            The method of the first registered toolkit that provides the
            requested method name

        Raises
        ------
        NotImplementedError
            if the requested method cannot be found among the registered toolkits

        Examples
        --------

        Create a molecule, and call the toolkit ``to_smiles()`` method directly

        >>> from openff.toolkit import Molecule
        >>> molecule = Molecule.from_smiles('Cc1ccccc1')
        >>> toolkit_registry = ToolkitRegistry([
        ...     OpenEyeToolkitWrapper,
        ...     RDKitToolkitWrapper,
        ...     AmberToolsToolkitWrapper
        ... ])
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


# Copied from https://github.com/openforcefield/openff-fragmenter/blob/4a290b866a8ed43eabcbd3231c62b01f0c6d7df6
# /openff/fragmenter/utils.py#L97-L123
@contextmanager
def toolkit_registry_manager(toolkit_registry: Union[ToolkitRegistry, ToolkitWrapper]):
    """
    A context manager that temporarily changes the :py:data:`GLOBAL_TOOLKIT_REGISTRY`.

    This can be useful in cases where one would otherwise need to otherwise
    manually specify the use of a specific :py:class:`ToolkitWrapper` or
    :py:class:`ToolkitRegistry` repeatedly in a block of code, or in cases where
    there isn't another way to switch the ``ToolkitWrapper`` used for a
    particular operation.

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
        raise NotImplementedError("Only ``ToolkitRegistry`` and ``ToolkitWrapper`` are supported.")

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
            GLOBAL_TOOLKIT_REGISTRY.register_toolkit(toolkit)


_toolkit_registry_manager = toolkit_registry_manager
