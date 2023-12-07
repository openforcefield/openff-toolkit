import importlib
import warnings
from typing import TYPE_CHECKING, List, Optional, Type

from openff.toolkit import Quantity, unit
from openff.toolkit.utils.base_wrapper import ToolkitWrapper
from openff.toolkit.utils.exceptions import (
    ChargeMethodUnavailableError,
    ToolkitUnavailableException,
)

if TYPE_CHECKING:
    from openff.toolkit.topology.molecule import FrozenMolecule, Molecule


__all__ = ("NAGLToolkitWrapper",)


class NAGLToolkitWrapper(ToolkitWrapper):
    _toolkit_name = "OpenFF NAGL"
    _toolkit_installation_instructions = (
        "See https://docs.openforcefield.org/projects/nagl/en/latest/installation.html"
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
        partial_charge_method: str,
        use_conformers: Optional[List["Quantity"]] = None,
        strict_n_conformers: bool = False,
        normalize_partial_charges: bool = True,
        _cls: Optional[Type["FrozenMolecule"]] = None,
    ):
        from openff.nagl import GNNModel
        from openff.nagl_models import validate_nagl_model_path

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        if use_conformers:
            warnings.warn(
                "`NAGLToolkitWrapper.assign_partial_charges` was passed optional argument "
                "`use_conformers` which will not be used. OpenFF NAGL does not generate "
                "conformers as part of assigning partial charges.",
                UserWarning,
                stacklevel=2,
            )

        if strict_n_conformers:
            warnings.warn(
                "`NAGLToolkitWrapper.assign_partial_charges` was passed optional argument "
                "`strict_n_conformers` which will not be used. OpenFF NAGL does not generate "
                "conformers as part of assigning partial charges.",
                UserWarning,
                stacklevel=2,
            )

        try:
            model_path = validate_nagl_model_path(model=partial_charge_method)
        except FileNotFoundError as error:
            raise ChargeMethodUnavailableError(
                f"Charge model {partial_charge_method} not supported by "
                f"{self.__class__.__name__}."
            ) from error

        model = GNNModel.load(model_path, eval_mode=True)
        charges = model.compute_property(
            molecule,
            as_numpy=True,
            readout_name="am1bcc_charges",
            check_domains=True,
            error_if_unsupported=True,
        )

        molecule.partial_charges = Quantity(
            charges.astype(float),
            unit.elementary_charge,
        )

        if normalize_partial_charges:
            molecule._normalize_partial_charges()
