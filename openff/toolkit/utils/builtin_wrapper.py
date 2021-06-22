"""
Built-in ToolkitWrapper for very basic functionality. This is intended for use in testing and not much more.
"""
__all__ = (
    "BuiltInToolkitWrapper",
    )

# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import numpy as np
from simtk import unit

# =============================================================================================
# LOCAL IMPORTS
# =============================================================================================

from . import base_wrapper

from .exceptions import (
    ChargeMethodUnavailableError,
    )

from .utils import (
    inherit_docstrings,
    )

# =============================================================================================
# Implementation
# =============================================================================================


@inherit_docstrings
class BuiltInToolkitWrapper(base_wrapper.ToolkitWrapper):
    """
    Built-in ToolkitWrapper for very basic functionality. This is intended for use in testing and not much more.

    .. warning :: This API is experimental and subject to change.
    """

    _toolkit_name = "Built-in Toolkit"
    _toolkit_installation_instructions = (
        "This toolkit is installed with the Open Force Field Toolkit and does "
        "not require additional dependencies."
    )

    def __init__(self):
        super().__init__()

        self._toolkit_file_read_formats = []
        self._toolkit_file_write_formats = []

    def assign_partial_charges(
        self,
        molecule,
        partial_charge_method=None,
        use_conformers=None,
        strict_n_conformers=False,
        _cls=None,
    ):
        """
        Compute partial charges with the built-in toolkit using simple arithmetic operations, and assign
        the new values to the partial_charges attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            Molecule for which partial charges are to be computed
        partial_charge_method: str, optional, default=None
            The charge model to use. One of ['zeros', 'formal_charge']. If None, 'formal_charge' will be used.
        use_conformers : iterable of simtk.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3) and dimension of distance. Optional, default = None
            Coordinates to use for partial charge calculation. If None, an appropriate number of conformers
            will be generated.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised
            instead of an Exception.
        _cls : class
            Molecule constructor

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not be handled by this toolkit

        IncorrectNumConformersError if strict_n_conformers is True and use_conformers is provided and specifies an
        invalid number of conformers for the requested method

        ChargeCalculationError if the charge calculation is supported by this toolkit, but fails
        """

        PARTIAL_CHARGE_METHODS = {
            "zeros": {"rec_confs": 0, "min_confs": 0, "max_confs": 0},
            "formal_charge": {"rec_confs": 0, "min_confs": 0, "max_confs": 0},
        }

        if partial_charge_method is None:
            partial_charge_method = "formal_charge"

        if _cls is None:
            from openff.toolkit.topology.molecule import Molecule

            _cls = Molecule

        # Make a temporary copy of the molecule, since we'll be messing with its conformers
        mol_copy = _cls(molecule)

        partial_charge_method = partial_charge_method.lower()
        if partial_charge_method not in PARTIAL_CHARGE_METHODS:
            raise ChargeMethodUnavailableError(
                f'Partial charge method "{partial_charge_method}"" is not supported by '
                f"the Built-in toolkit. Available charge methods are "
                f"{list(PARTIAL_CHARGE_METHODS.keys())}"
            )

        if use_conformers is None:
            # Note that this refers back to the GLOBAL_TOOLKIT_REGISTRY by default, since
            # BuiltInToolkitWrapper can't generate conformers
            mol_copy.generate_conformers(
                n_conformers=PARTIAL_CHARGE_METHODS[partial_charge_method]["rec_confs"]
            )
        else:
            mol_copy._conformers = None
            for conformer in use_conformers:
                mol_copy._add_conformer(conformer)
            self._check_n_conformers(
                mol_copy,
                partial_charge_method=partial_charge_method,
                min_confs=0,
                max_confs=0,
                strict_n_conformers=strict_n_conformers,
            )

        partial_charges = unit.Quantity(
            np.zeros((molecule.n_particles)), unit.elementary_charge
        )
        if partial_charge_method == "zeroes":
            pass
        elif partial_charge_method == "formal_charge":
            for part_idx, particle in enumerate(molecule.particles):
                partial_charges[part_idx] = particle.formal_charge

        molecule.partial_charges = partial_charges

