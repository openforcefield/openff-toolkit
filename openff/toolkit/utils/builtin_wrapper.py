"""
Built-in ToolkitWrapper for very basic functionality. Intended for testing and not much more.
"""

__all__ = ("BuiltInToolkitWrapper",)

from typing import TYPE_CHECKING, Optional

from openff.toolkit import Quantity, unit
from openff.toolkit.utils import base_wrapper
from openff.toolkit.utils.exceptions import ChargeMethodUnavailableError
from openff.toolkit.utils.utils import inherit_docstrings

if TYPE_CHECKING:
    from openff.toolkit.topology.molecule import FrozenMolecule


@inherit_docstrings
class BuiltInToolkitWrapper(base_wrapper.ToolkitWrapper):
    """
    Built-in ToolkitWrapper for very basic functionality. Intended for testing and not much more.

    .. warning :: This API is experimental and subject to change.
    """

    _toolkit_name = "Built-in Toolkit"
    _toolkit_installation_instructions = (
        "This toolkit is installed with the Open Force Field Toolkit and does "
        "not require additional dependencies."
    )
    _supported_charge_methods = {
        "zeros": {"rec_confs": 0, "min_confs": 0, "max_confs": 0},
        "formal_charge": {"rec_confs": 0, "min_confs": 0, "max_confs": 0},
    }

    PARTIAL_CHARGE_METHODS = _supported_charge_methods

    def __init__(self):
        super().__init__()

    def assign_partial_charges(
        self,
        molecule: "FrozenMolecule",
        partial_charge_method: Optional[str] = None,
        use_conformers: Optional[Quantity] = None,
        strict_n_conformers: bool = False,
        normalize_partial_charges: bool = True,
        _cls: Optional[type] = None,
    ):
        """

        Compute partial charges with the built-in toolkit using simple arithmetic operations,
        and assign the new values to the partial_charges attribute.


        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule
            Molecule for which partial charges are to be computed
        partial_charge_method
            The charge model to use. One of ['zeros', 'formal_charge']. If None, 'formal_charge'
            will be used.
        use_conformers
            (n_atoms, 3) and dimension of distance. Optional, default = None
            Coordinates to use for partial charge calculation. If None, an appropriate number
            of conformers will be generated.
        strict_n_conformers
            Whether to raise an exception if an invalid number of conformers is provided for the
            given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised
            instead of an Exception.
        normalize_partial_charges
            Whether to offset partial charges so that they sum to the total formal charge of the molecule.
            This is used to prevent accumulation of rounding errors when the partial charge generation method has
            low precision.
        _cls
            Molecule constructor

        Raises
        ------
        ChargeMethodUnavailableError if this toolkit cannot handle the requested charge method

        IncorrectNumConformersError if strict_n_conformers is True and use_conformers is provided
        and specifies an invalid number of conformers for the requested method

        ChargeCalculationError if the charge calculation is supported by this toolkit, but fails
        """

        if partial_charge_method is None:
            partial_charge_method = "formal_charge"

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        # Make a temporary copy of the molecule, since we'll be messing with its conformers
        mol_copy = _cls(molecule)

        partial_charge_method = partial_charge_method.lower()
        if partial_charge_method not in self._supported_charge_methods:
            raise ChargeMethodUnavailableError(
                f'Partial charge method "{partial_charge_method}"" is not supported by '
                f"the Built-in toolkit. Available charge methods are "
                f"{self._supported_charge_methods}"
            )

        if use_conformers is None:
            # Note that this refers back to the GLOBAL_TOOLKIT_REGISTRY by default, since
            # BuiltInToolkitWrapper can't generate conformers
            mol_copy.generate_conformers(
                n_conformers=self._supported_charge_methods[partial_charge_method][
                    "rec_confs"
                ]
            )
        else:
            mol_copy._conformers = None
            for conformer in use_conformers:  # type: ignore[attr-defined]
                mol_copy._add_conformer(conformer)
            self._check_n_conformers(
                mol_copy,
                partial_charge_method=partial_charge_method,
                min_confs=0,
                max_confs=0,
                strict_n_conformers=strict_n_conformers,
            )

        if partial_charge_method == "zeros":
            partial_charges = [0.0] * molecule.n_atoms

        elif partial_charge_method == "formal_charge":
            partial_charges = [float(atom.formal_charge.m) for atom in molecule.atoms]

        molecule.partial_charges = Quantity(partial_charges, unit.elementary_charge)

        if normalize_partial_charges:
            molecule._normalize_partial_charges()
