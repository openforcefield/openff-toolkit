import importlib
import warnings
from typing import TYPE_CHECKING, List, Optional, Type

from openff.toolkit.utils.base_wrapper import ToolkitWrapper
from openff.toolkit.utils.exceptions import ToolkitUnavailableException

if TYPE_CHECKING:
    from openff.units import Quantity

    from openff.toolkit.topology.molecule import FrozenMolecule, Molecule


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
        use_conformers: Optional[List["Quantity"]] = None,
        strict_n_conformers: bool = False,
        normalize_partial_charges: bool = True,
        _cls: Optional[Type["FrozenMolecule"]] = None,
    ):
        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        if use_conformers:
            warnings.warn(
                "`_NAGLToolkitWrapper.assign_partial_charges` was passed optional orgument "
                "`use_conformers` which will not be used. OpenFF NAGL does not generate "
                "conformers as part of assigning partial charges.",
                UserWarning,
            )

        if strict_n_conformers:
            warnings.warn(
                "`_NAGLToolkitWrapper.assign_partial_charges` was passed optional orgument "
                "`strict_n_conformers` which will not be used. OpenFF NAGL does not generate "
                "conformers as part of assigning partial charges.",
                UserWarning,
            )

        if partial_charge_method == "_nagl_am1bccelf10":
            from openff.nagl.data.files import EXAMPLE_AM1BCC_MODEL
            from openff.nagl.nn._models import GNNModel

            model = GNNModel.load(EXAMPLE_AM1BCC_MODEL, eval_mode=True)

            model.compute_property(molecule, as_numpy=True)

        else:
            # This should be a more specific exception that inherits from ValueError?
            raise ValueError(
                f"Charge model {partial_charge_method} not supported by "
                f"{self.__class__.__name__}."
            )

        if normalize_partial_charges:
            molecule._normalize_partial_charges()
