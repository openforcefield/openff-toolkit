"""
Base class for toolkit wrappers. Defines the public API and some shared methods
"""

__all__ = ("ToolkitWrapper",)

from functools import wraps
from typing import TYPE_CHECKING, Optional

from openff.toolkit.utils.constants import DEFAULT_AROMATICITY_MODEL
from openff.toolkit.utils.exceptions import (
    IncorrectNumConformersError,
    IncorrectNumConformersWarning,
    ToolkitUnavailableException,
)

if TYPE_CHECKING:
    from openff.toolkit.topology.molecule import Molecule


def _mol_to_ctab_and_aro_key(
    self, molecule: "Molecule", aromaticity_model=DEFAULT_AROMATICITY_MODEL
) -> str:
    return f"{molecule.ordered_connection_table_hash()}-{aromaticity_model}"


class ToolkitWrapper:
    """
    Toolkit wrapper base class.

    .. warning :: This API is experimental and subject to change.
    """

    _is_available: Optional[bool] = None  # True if toolkit is available
    _toolkit_version: Optional[str] = None
    _toolkit_name: Optional[str] = None  # Name of the toolkit
    _toolkit_installation_instructions: Optional[
        str
    ] = None  # Installation instructions for the toolkit

    # @staticmethod
    # TODO: Right now, to access the class definition, I have to make this a classmethod
    # and thereby call it with () on the outermost decorator. Is this wasting time? Are we caching
    # the is_available results?
    @classmethod
    def requires_toolkit(cls):  # remember cls is a ToolkitWrapper subclass here
        def decorator(func):
            @wraps(func)
            def wrapped_function(*args, **kwargs):
                if not cls.is_available():
                    msg = "This function requires the {} toolkit".format(
                        cls._toolkit_name
                    )
                    raise ToolkitUnavailableException(msg)
                value = func(*args, **kwargs)
                return value

            return wrapped_function

        return decorator

    @property
    # @classmethod
    def toolkit_name(self):
        """
        Return the name of the toolkit wrapped by this class as a str

        .. warning :: This API is experimental and subject to change.

        Returns
        -------
        toolkit_name : str
            The name of the wrapped toolkit

        """
        return self.__class__._toolkit_name

    @property
    # @classmethod
    def toolkit_installation_instructions(self):
        """
        Instructions on how to install the wrapped toolkit.
        """
        return self._toolkit_installation_instructions

    # @classmethod
    @property
    def toolkit_file_read_formats(self):
        """
        List of file formats that this toolkit can read.
        """
        return self._toolkit_file_read_formats

    # @classmethod
    @property
    def toolkit_file_write_formats(self):
        """
        List of file formats that this toolkit can write.
        """
        return self._toolkit_file_write_formats

    @classmethod
    def is_available(cls):
        """
        Check whether the corresponding toolkit can be imported

        Returns
        -------
        is_installed : bool
            True if corresponding toolkit is installed, False otherwise.

        """
        return NotImplementedError

    @property
    def toolkit_version(self):
        """
        Return the version of the wrapped toolkit as a str

        .. warning :: This API is experimental and subject to change.

        Returns
        -------
        toolkit_version : str
            The version of the wrapped toolkit

        """
        return self._toolkit_version

    def from_file(self, file_path, file_format, allow_undefined_stereo=False):
        """
        Return an openff.toolkit.topology.Molecule from a file using this toolkit.

        Parameters
        ----------
        file_path : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check
            ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if any molecules contain undefined stereochemistry.
        _cls : class
            Molecule constructor
        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.

        """
        return NotImplementedError

    def from_file_obj(
        self, file_obj, file_format, allow_undefined_stereo=False, _cls=None
    ):
        """
        Return an openff.toolkit.topology.Molecule from a file-like object (an object with
        a ".read()" method using this toolkit.


        Parameters
        ----------
        file_obj : file-like object
            The file-like object to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check
            ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if any molecules contain undefined stereochemistry.
            If false, the function skips loading the molecule.
        _cls : class
            Molecule constructor

        Returns
        -------
        molecules : Molecule or list of Molecules
            a list of Molecule objects is returned.
        """
        return NotImplementedError

    def _check_n_conformers(
        self,
        molecule: "Molecule",
        partial_charge_method: Optional[str] = None,
        min_confs: Optional[int] = None,
        max_confs: Optional[int] = None,
        strict_n_conformers: bool = False,
    ):
        """
        Private method for validating the number of conformers on a molecule prior to partial
        charge calculation

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method : str, optional, default=None
            The name of the charge method being used
        min_confs : int, optional, default=None
            The minimum number of conformers required to use this charge method
        max_confs : int, optional, default=None
            The maximum number of conformers required to use this charge method
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided.
            If this is False and an invalid number of conformers is found, a warning will be raised.

        Raises
        ------
        IncorrectNumConformersError
            If the wrong number of conformers is attached to the input molecule, and
            strict_n_conformers is True.
        """
        import warnings

        n_confs = molecule.n_conformers
        wrong_confs_msg = (
            f"Molecule '{molecule}' has {n_confs} conformers, "
            f"but charge method '{partial_charge_method}' expects"
        )
        exception_suffix = (
            "You can disable this error by setting `strict_n_conformers=False' "
            "when calling 'molecule.assign_partial_charges'."
        )
        # If there's no n_confs filter, then this molecule automatically passes
        if min_confs is None and max_confs is None:
            return
        # If there's constraints on both ends, check both limits
        elif min_confs is not None and max_confs is not None:
            if not (min_confs <= n_confs <= max_confs):
                if min_confs == max_confs:
                    wrong_confs_msg += f" exactly {min_confs}."
                else:
                    wrong_confs_msg += f" between {min_confs} and {max_confs}."

            else:
                return
        # If there's only a max constraint, check that
        elif min_confs is not None and max_confs is None:
            if not (min_confs <= n_confs):
                wrong_confs_msg += f" at least {min_confs}."
            else:
                return
        # If there's only a maximum constraint, check that
        elif min_confs is None and max_confs is not None:
            if not (n_confs <= max_confs):
                wrong_confs_msg += f" at most {max_confs}."
            else:
                return
        # If we've made it this far, the molecule has the wrong number of conformers
        if strict_n_conformers:
            wrong_confs_msg += exception_suffix
            raise IncorrectNumConformersError(wrong_confs_msg)
        else:
            warnings.warn(wrong_confs_msg, IncorrectNumConformersWarning)

    def __repr__(self):
        return (
            f"ToolkitWrapper around {self.toolkit_name} version {self.toolkit_version}"
        )
