import importlib
from typing import TYPE_CHECKING, List, Optional

from openff.toolkit.utils.base_wrapper import ToolkitWrapper
from openff.toolkit.utils.exceptions import ToolkitUnavailableException

if TYPE_CHECKING:
    from openff.toolkit.topology.molecule import Molecule


__all__ = ("_NAGLToolkitWrapper",)


class _NAGLToolkitWrapper(ToolkitWrapper):

    _toolkit_name = "OpenFF NAGL"
    _toolkit_installation_instructions = (
        "Installation instructions are in flux. See for updates: "
        "https://github.com/openforcefield/openff-nagl#installation"
    )

    def __init__(self):
        super().__init__()

        self._toolkit_file_read_formats = []
        self._toolkit_file_write_formats = []

        if not self.is_available():
            raise ToolkitUnavailableException(
                f"The required toolkit {self._toolkit_name} is not "
                f"available. {self._toolkit_installation_instructions}"
            )
        else:
            from openff.nagl import __version__ as nagl_version

            self._toolkit_version = nagl_version

    @classmethod
    def is_available(cls) -> bool:
        if cls._is_available is None:
            try:
                importlib.import_module("openff.nagl")
            except ImportError:
                cls._is_available = False
            else:
                cls._is_available = True
        return cls._is_available

    def assign_partial_charges(
        self,
        molecule: "Molecule",
        partial_charge_method: str = "_nagl_am1bccelf10",
        use_conformers: Optional[List] = None,
        strict_n_conformers: bool = False,
        normalize_partial_charges: bool = True,
        _cls=None,
    ):

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        if use_conformers:
            raise Exception

        if strict_n_conformers:
            raise Exception

        if partial_charge_method == "_nagl_am1bccelf10":
            # I can't yet tell from the source code where in NAGL a high-level
            # `def charge(Molecule, ...):` function is implemented, but I'm sure
            # it's there somewhere
            from openff.nagl import charge

            charge(molecule)

        else:
            raise Exception

        if normalize_partial_charges:
            molecule._normalize_partial_charges()
