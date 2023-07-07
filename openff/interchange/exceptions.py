"""Custom exceptions used in Interchange."""


class SMIRNOFFParameterAttributeNotImplementedError(Exception):
    """
    Exception for when a parameter attribute is supported by the SMIRNOFF specification but not yet implemented.

    For example, this was raised when k_bondorder (used in bond order-based interpolation of force constants)
    before the behavior was supported.
    """


class SMIRNOFFHandlersNotImplementedError(Exception):
    """
    Exception for when some parameter handlers in the SMIRNOFF specification are not implemented here.
    """


class SMIRNOFFVersionNotSupportedError(Exception):
    """
    Exception for when a parameter handler's version is not supported.
    """


class InvalidParameterHandlerError(ValueError):
    """
    Generic exception for mismatch between expected and found ParameterHandler types.
    """


class InvalidBoxError(ValueError):
    """
    Generic exception for errors reading box data.
    """


class InvalidTopologyError(ValueError):
    """
    Generic exception for errors reading chemical topology data.
    """


class NonbondedEnergyError(AssertionError):
    """
    Exception for when non-bonded energies computed from different objects differ.
    """


class InvalidExpressionError(ValueError):
    """
    Exception for when an expression cannot safely be interpreted.
    """


class UnsupportedCutoffMethodError(BaseException):
    """
    Exception for a cutoff method that is invalid or not supported by an engine.
    """


class UnimplementedCutoffMethodError(BaseException):
    """
    Exception for a cutoff method that should be supported but it not yet implemented.
    """


class UnsupportedMixingRuleError(BaseException):
    """
    Raised when attempting to use a mixing rule is invalid, unsupported or otherwise not implemented.
    """


class UnsupportedParameterError(ValueError):
    """
    Exception for parameters having unsupported values, i.e. non-1.0 idivf.
    """


class UnsupportedBoxError(ValueError):
    """
    Exception for processing an unsupported box, probably non-orthogonal.
    """


class UnsupportedImportError(BaseException):
    """
    Generic exception for attempting to import from an unsupported data format.
    """


class UnsupportedExportError(BaseException):
    """
    Exception for attempting to write to an unsupported file format.
    """


class UnsupportedCombinationError(BaseException):
    """General exception for something going wrong in Interchange object combination."""


class PluginCompatibilityError(BaseException):
    """A plugin is incompatible with the current version of Interchange in the way it is called."""


class CannotSetSwitchingFunctionError(BaseException):
    """
    Unable to set a switching function.
    """


class CannotInferEnergyError(BaseException):
    """
    Failed to infer a physical interpretation of this energy term.
    """


class CannotInferNonbondedEnergyError(CannotInferEnergyError):
    """
    Failed to infer a physical interpretation of this non-bonded energy.
    """


class InvalidWriterError(BaseException):
    """
    An unknown or unimplemnted writer was specified.
    """


class ConversionError(BaseException):
    """
    Base exception for error handling during object conversion.
    """


class MissingBoxError(BaseException):
    """
    Exception for when box vectors are needed but missing.
    """


class MissingPositionsError(BaseException):
    """
    Exception for when positions are needed but missing.
    """


class MissingParameterHandlerError(BaseException):
    """
    Exception for when a parameter handler is requested but not found.
    """


class MissingParametersError(BaseException):
    """
    Exception for when parameters are needed but missing.
    """


class MissingBondOrdersError(BaseException):
    """
    Exception for when a parameter handler needs fractional bond orders but they are missing.
    """


class DuplicateMoleculeError(ValueError):
    """
    Exception for when molecules are not unique.
    """


class MissingUnitError(ValueError):
    """
    Exception for data missing a unit tag.
    """


class UnitValidationError(ValueError):
    """
    Exception for bad behavior when validating unit-tagged data.
    """


class NonbondedCompatibilityError(BaseException):
    """
    Exception for unsupported combination of nonbonded methods.
    """


class MissingNonbondedCompatibilityError(BaseException):
    """
    Exception for uncovered combination of nonbonded methods.
    """


class InternalInconsistencyError(BaseException):
    """
    Fallback exception for bad behavior releating to a self-inconsistent internal state.

    These should not be reached but are raised to safeguard against problematic edge cases silently passing.
    """


class AmberError(BaseException):
    """
    Base exception for handling Amber-related errors.
    """


class AmberExecutableNotFoundError(AmberError):
    """
    Exception for when an Amber-related excutable is not found.
    """


class SanderError(AmberError):
    """
    Exception for when a sander subprocess fails.
    """


class GMXError(BaseException):
    """
    Exception for when a GROMACS subprocess fails.
    """


class GMXNotFoundError(GMXError):
    """
    Exception for when no GROMACS executable is found.
    """


class GMXGromppError(GMXError):
    """
    Exception for when `gmx grompp` fails.
    """


class GMXMdrunError(GMXError):
    """
    Exception for when `gmx mdrun` fails.
    """


class LAMMPSError(BaseException):
    """
    Base exception for handling LAMMPS-related errors.
    """


class LAMMPSNotFoundError(LAMMPSError):
    """
    Exception for when no LAMMPS executable is found.
    """


class LAMMPSRunError(LAMMPSError):
    """
    Exception for when a LAMMPS subprocess fails.
    """


class EnergyError(BaseException):
    """
    Base class for energies in reports not matching.
    """


class InvalidEnergyError(BaseException):
    """
    Energy type not understood.
    """


class IncompatibleTolerancesError(BaseException):
    """
    Exception for when one report has a value for an energy group but the other does not.
    """


class VirtualSiteTypeNotImplementedError(BaseException):
    """
    Raised when this type of virtual site is not yet implemented.
    """


class NonIntegralMoleculeChargeError(BaseException):
    """
    Exception raised when the partial charges on a molecule do not sum up to its formal charge.
    """


class MissingPartialChargesError(BaseException):
    """
    A molecule is missing partial charges.
    """


class UnassignedValenceError(BaseException):
    """Exception raised when there are valence terms for which a ParameterHandler did not find matches."""


class UnassignedBondError(UnassignedValenceError):
    """Exception raised when there are bond terms for which a ParameterHandler did not find matches."""


class UnassignedAngleError(UnassignedValenceError):
    """Exception raised when there are angle terms for which a ParameterHandler did not find matches."""


class UnassignedTorsionError(UnassignedValenceError):
    """Exception raised when there are torsion terms for which a ParameterHandler did not find matches."""


class MissingValenceError(BaseException):
    """Exception raised when there are valence interactions for which no parameters are found."""


class MissingBondError(UnassignedValenceError):
    """Exception raised when there exists a bond for which no parameters are found."""


class MissingAngleError(UnassignedValenceError):
    """Exception raised when there exists an angle for which no parameters are found."""


class MissingTorsionError(UnassignedValenceError):
    """Exception raised when there exists a torsion for which no parameters are found."""


class PACKMOLRuntimeError(Exception):
    """Exception raised when PACKMOL fails to execute / converge."""


class PACKMOLValueError(Exception):
    """Exception raised when a bad input is passed to a PACKMOL wrapper."""


class ExperimentalFeatureException(Exception):
    """Exception raised when an experimental feature is used without opt-in."""
