"""
environment.py

Lightweight module for validating SMIRKS using Open Force Field ToolkitWrappers.

AUTHORS

Caitlin Bannan <bannanc@uci.edu> (Original author), Mobley Lab, University of California Irvine,
Jeff Wagner <jeffrey.wagner@openforcefield.org> (refactored to use ToolkitWrappers),
with contributions from John Chodera, Memorial Sloan Kettering Cancer Center
and David Mobley, UC Irvine.

"""

__all__ = [
    "SMIRKSMismatchError",
    "SMIRKSParsingError",
    "ChemicalEnvironment",
    "AtomChemicalEnvironment",
    "BondChemicalEnvironment",
    "AngleChemicalEnvironment",
    "TorsionChemicalEnvironment",
    "ImproperChemicalEnvironment",
]

import warnings
from typing import Optional

from openff.toolkit.utils.exceptions import SMIRKSMismatchError, SMIRKSParsingError
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, ToolkitWrapper


class ChemicalEnvironmentDeprecationWarning(UserWarning):
    """Warning for deprecated portions of the Molecule API."""


class ChemicalEnvironment:
    """Chemical environment abstract base class used for validating SMIRKS"""

    _expected_type: Optional[str] = None

    def __init__(
        self,
        smirks=None,
        label=None,
        validate_parsable=True,
        validate_valence_type=True,
        toolkit_registry=None,
    ):
        """Initialize a chemical environment abstract base class.

        smirks = string, optional
            if smirks is not None, a chemical environment is built
            from the provided SMIRKS string
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        validate_parsable: bool, optional, default=True
            If specified, ensure the provided smirks is parsable
        validate_valence_type : bool, optional, default=True
            If specified, ensure the tagged atoms are appropriate to the specified valence type
        toolkit_registry = string or ToolkitWrapper or ToolkitRegistry. Default = None
            Either a ToolkitRegistry, ToolkitWrapper, or the strings 'openeye' or 'rdkit',
            indicating the backend to use for validating the correct
            connectivity of the SMIRKS during initialization. If None,
            this function will use the GLOBAL_TOOLKIT_REGISTRY

        Raises
        ------
        SMIRKSParsingError
            if smirks was unparsable
        SMIRKSMismatchError
            if smirks did not have expected connectivity between tagged atoms
            and validate_valence_type=True
        """
        warnings.warn(
            "ChemicalEnvironment is deprecated and will be removed in a future "
            "release. If you rely on the functionality in this class, please "
            "open an issue on the openff-toolkit GitHub.",
            ChemicalEnvironmentDeprecationWarning,
        )

        # Support string input for toolkit names for legacy reasons
        if toolkit_registry == "openeye":
            from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper

            toolkit_registry = OpenEyeToolkitWrapper()
        elif toolkit_registry == "rdkit":
            from openff.toolkit.utils.toolkits import RDKitToolkitWrapper

            toolkit_registry = RDKitToolkitWrapper()

        self.smirks = smirks
        self.label = label
        if validate_parsable or validate_valence_type:
            self.validate(
                validate_valence_type=validate_valence_type,
                toolkit_registry=toolkit_registry,
            )

    def validate(self, validate_valence_type=True, toolkit_registry=None):
        """
        Returns True if the underlying smirks is the correct valence type, False otherwise. If the expected type
        is None, this method always returns True.

        validate_valence_type : bool, optional, default=True
            If specified, ensure the tagged atoms are appropriate to the specified valence type
        toolkit_registry = ToolkitWrapper or ToolkitRegistry. Default = None
            Either a ToolkitRegistry or ToolkitWrapper,
            indicating the backend to use for validating the correct
            connectivity of the SMIRKS during initialization. If None,
            this function will use the GLOBAL_TOOLKIT_REGISTRY

        Raises
        ------
        SMIRKSParsingError
            if smirks was unparsable
        SMIRKSMismatchError
            if smirks did not have expected connectivity between tagged atoms
            and validate_valence_type=True
        """
        perceived_type = self.get_type(toolkit_registry=toolkit_registry)
        if validate_valence_type and self._expected_type is not None:
            if perceived_type != self._expected_type:
                raise SMIRKSMismatchError(
                    f"{self.__class__} expected '{self._expected_type}' chemical environment, but "
                    f"smirks was set to '{self.smirks}', which is type '{perceived_type}'"
                )

    @classmethod
    def validate_smirks(
        cls,
        smirks,
        validate_parsable=True,
        validate_valence_type=True,
        toolkit_registry=None,
    ):
        """
        Check the provided SMIRKS string is valid, and if requested, tags atoms appropriate to the
        specified valence type.

        Parameters
        ----------
        smirks : str
            The SMIRKS expression to validate
        validate_parsable: bool, optional, default=True
            If specified, ensure the provided smirks is parsable
        validate_valence_type : bool, optional, default=True
            If specified, ensure the tagged atoms are appropriate to the specified valence type
        toolkit_registry = string or ToolkitWrapper or ToolkitRegistry. Default = None
            Either a ToolkitRegistry, ToolkitWrapper, or the strings 'openeye' or 'rdkit',
            indicating the backend to use for validating the correct
            connectivity of the SMIRKS during initialization. If None,
            this function will use the GLOBAL_TOOLKIT_REGISTRY

        Raises
        ------
        SMIRKSParsingError
            if smirks was unparsable
        SMIRKSMismatchError
            if smirks did not have expected connectivity between tagged atoms
            and validate_valence_type=True
        """
        cls(
            smirks,
            validate_parsable=validate_parsable,
            validate_valence_type=validate_valence_type,
            toolkit_registry=toolkit_registry,
        )

    def get_type(self, toolkit_registry=None):
        """
        Return the valence type implied by the connectivity of the bound atoms in this ChemicalEnvironment.

        Parameters
        -----------
        toolkit_registry : openff.toolkit.utils.ToolkitRegistry or openff.toolkit.utils.ToolkitWrapper
            The cheminformatics toolkit to use for parsing the smirks

        Returns
        -------
        valence_type : str
            One of "Atom", "Bond", "Angle", "ProperTorsion", "ImproperTorsion", or None.
            If tagged atoms are not connected in a known pattern this method will return None.

        Raises
        ------
        SMIRKSParsingError
            if smirks was unparsable
        """

        # Query a toolkit wrapper for substructure type
        if toolkit_registry is None:
            toolkit_registry = GLOBAL_TOOLKIT_REGISTRY

        if isinstance(toolkit_registry, ToolkitWrapper):
            unique_tags, conn = toolkit_registry.get_tagged_smarts_connectivity(
                self.smirks
            )
        else:
            unique_tags, conn = toolkit_registry.call(
                "get_tagged_smarts_connectivity", self.smirks
            )

        if unique_tags == (1,):
            if len(conn) == 0:
                return "Atom"
        if unique_tags == (1, 2):
            if (1, 2) in conn:
                return "Bond"
        elif unique_tags == (1, 2, 3):
            if (1, 2) in conn and (2, 3) in conn:
                return "Angle"
        elif unique_tags == (1, 2, 3, 4):
            if (1, 2) in conn and (2, 3) in conn and (2, 3) in conn and (3, 4) in conn:
                return "ProperTorsion"
            elif (1, 2) in conn and (2, 3) in conn and (2, 4) in conn:
                return "ImproperTorsion"
        else:
            return None


class AtomChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching one labeled atom."""

    _expected_type = "Atom"


class BondChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching two labeled atoms (or a bond)."""

    _expected_type = "Bond"


class AngleChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching three marked atoms (angle)."""

    _expected_type = "Angle"


class TorsionChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching four marked atoms (torsion)."""

    _expected_type = "ProperTorsion"


class ImproperChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching four marked atoms (improper)."""

    _expected_type = "ImproperTorsion"
