import warnings

deprecated_names = ["MessageException", "ParseError"]

# TODO: Remove in January 2022
# TODO: When _DepcreatedMessageException is removed, this should be removed. Also
#       include in release notes that ParseError no longer exists
def __getattr__(name):
    if name in deprecated_names:
        warnings.filterwarnings("default", category=DeprecationWarning)
        warning_msg = f"{name} is DEPRECATED and will be removed in a future release of the OpenFF Toolkit."
        if name == "MessageException":
            warning_msg += (
                " All custom exceptions now inherit from OpenFFToolkitException, "
                "which should be used as a replacement for MessageException. Import it via "
                "`from openff.toolkit.utils.exceptions import OpenFFToolkitException`."
            )
        warnings.warn(warning_msg, DeprecationWarning)
        return globals()[f"_Deprecated{name}"]
    raise AttributeError(f"module {__name__} has no attribute {name}")


# TODO: Remove in January 2022
class _DeprecatedMessageException(Exception):
    """DEPRECATED: A base class for exceptions that print out a string given in their constructor"""

    def __init__(self, msg):
        super().__init__(self, msg)
        self.msg = msg

    def __str__(self):
        return self.msg


# TODO: Remove _DeprecatedMessageException (the import here and its definition above)
# and make this exception only import from Exception
class OpenFFToolkitException(_DeprecatedMessageException):
    """Base exception for custom exceptions raised by the OpenFF Toolkit"""

    def __init__(self, msg):
        super().__init__(msg)
        self.msg = msg

    def __str__(self):
        return self.msg


class IncompatibleUnitError(OpenFFToolkitException):
    """
    Exception for when a parameter is in the wrong units for a ParameterHandler's unit system
    """


class MissingDependencyError(OpenFFToolkitException):
    """
    Exception for when an optional dependency is needed but not installed

    """

    def __init__(self, package_name):
        self.msg = (
            f"Missing dependency {package_name}. Try installing it "
            f"with\n\n$ conda install {package_name} -c conda-forge"
        )

        super().__init__(self.msg)


class MissingPackageError(OpenFFToolkitException):
    """This function requires a package that is not installed."""


class ToolkitUnavailableException(OpenFFToolkitException):
    """The requested toolkit is unavailable."""

    # TODO: Allow toolkit to be specified and used in formatting/printing exception.


class LicenseError(ToolkitUnavailableException):
    """This function requires a license that cannot be found."""


class InvalidToolkitError(OpenFFToolkitException):
    """A non-toolkit object was received when a toolkit object was expected"""


class InvalidToolkitRegistryError(OpenFFToolkitException):
    """An object other than a ToolkitRegistry or toolkit wrapper was received"""


class UndefinedStereochemistryError(OpenFFToolkitException):
    """A molecule was attempted to be loaded with undefined stereochemistry"""


class GAFFAtomTypeWarning(RuntimeWarning):
    """A warning raised if a loaded mol2 file possibly uses GAFF atom types."""


class ChargeMethodUnavailableError(OpenFFToolkitException):
    """A toolkit does not support the requested partial_charge_method combination"""


class IncorrectNumConformersError(OpenFFToolkitException):
    """The requested partial_charge_method expects a different number of conformers than was provided"""


class IncorrectNumConformersWarning(Warning):
    """The requested partial_charge_method expects a different number of conformers than was provided"""


class ChargeCalculationError(OpenFFToolkitException):
    """An unhandled error occurred in an external toolkit during charge calculation"""


class ConformerGenerationError(OpenFFToolkitException):
    """Conformer generation via a wrapped toolkit failed."""


class InvalidIUPACNameError(OpenFFToolkitException):
    """Failed to parse IUPAC name"""


class AntechamberNotFoundError(OpenFFToolkitException):
    """The antechamber executable was not found"""


class SMILESParseError(OpenFFToolkitException, ValueError):
    """The record could not be parsed into the given format"""


class RadicalsNotSupportedError(OpenFFToolkitException):
    """The OpenFF Toolkit does not currently support parsing molecules with radicals."""


class InvalidConformerError(OpenFFToolkitException):
    """
    This error is raised when the conformer added to the molecule
    has a different connectivity to that already defined.
    or any other conformer related issues.
    """


# TODO: Remove in favor of SMILESParseError? They are used in different modules
class SmilesParsingError(OpenFFToolkitException):
    """
    This error is raised when parsing a SMILES string results in an error.
    """


class NotAttachedToMoleculeError(OpenFFToolkitException):
    """Exception for when a component does not belong to a Molecule object, but is queried"""


class DuplicateUniqueMoleculeError(OpenFFToolkitException):
    """
    Exception for when the user provides indistinguishable unique molecules when trying to identify atoms from a PDB
    """


class NotBondedError(OpenFFToolkitException):
    """
    Exception for when a function requires a bond between two atoms, but none is present
    """


class InvalidBoxVectorsError(OpenFFToolkitException):
    """
    Exception for setting invalid box vectors
    """


class InvalidPeriodicityError(OpenFFToolkitException):
    """
    Exception for setting invalid periodicity
    """


class MissingUniqueMoleculesError(OpenFFToolkitException):
    """
    Exception for when unique_molecules is required but not found
    """


class SMIRKSMismatchError(OpenFFToolkitException):
    """
    Exception for cases where smirks are inappropriate
    for the environment type they are being parsed into
    """


class SMIRKSParsingError(OpenFFToolkitException):
    """
    Exception for when SMIRKS are not parseable for any environment
    """


class ParameterHandlerRegistrationError(OpenFFToolkitException):
    """
    Exception for errors in ParameterHandler registration
    """


class SMIRNOFFVersionError(OpenFFToolkitException):
    """
    Exception thrown when an incompatible SMIRNOFF version data structure is attempted to be read.
    """


class SMIRNOFFAromaticityError(OpenFFToolkitException):
    """
    Exception thrown when an incompatible SMIRNOFF aromaticity model is checked for compatibility.
    """


# TODO: Remove in January 2022
class _DeprecatedParseError(_DeprecatedMessageException):
    """DEPRECATED: Error for when a SMIRNOFF data structure is not parseable by a ForceField"""


class SMIRNOFFParseError(OpenFFToolkitException, _DeprecatedParseError):
    """
    Error for when a SMIRNOFF data structure is not parseable by a ForceField
    """

    # TODO: Remove ParseError altogether by v0.11.0


class PartialChargeVirtualSitesError(OpenFFToolkitException):
    """
    Exception thrown when partial charges cannot be computed for a Molecule because the ForceField applies virtual sites.
    """


class SMIRNOFFSpecError(OpenFFToolkitException):
    """
    Exception for when data is noncompliant with the SMIRNOFF data specification.
    """


class SMIRNOFFSpecUnimplementedError(OpenFFToolkitException):
    """
    Exception for when a portion of the SMIRNOFF specification is not yet implemented.
    """


class FractionalBondOrderInterpolationMethodUnsupportedError(OpenFFToolkitException):
    """
    Exception for when an unsupported fractional bond order interpolation assignment method is called.
    """


class NotEnoughPointsForInterpolationError(OpenFFToolkitException):
    """Exception for when less than two points are provided for interpolation"""


class IncompatibleParameterError(OpenFFToolkitException):
    """
    Exception for when a set of parameters is scientifically/technically incompatible with another
    """


class UnassignedValenceParameterException(OpenFFToolkitException):
    """Exception raised when there are valence terms for which a ParameterHandler can't find parameters."""


class UnassignedBondParameterException(UnassignedValenceParameterException):
    """Exception raised when there are bond terms for which a ParameterHandler can't find parameters."""


class UnassignedAngleParameterException(UnassignedValenceParameterException):
    """Exception raised when there are angle terms for which a ParameterHandler can't find parameters."""


class UnassignedProperTorsionParameterException(UnassignedValenceParameterException):
    """Exception raised when there are proper torsion terms for which a ParameterHandler can't find parameters."""


class UnassignedMoleculeChargeException(OpenFFToolkitException):
    """Exception raised when no charge method is able to assign charges to a molecule."""


class NonintegralMoleculeChargeException(OpenFFToolkitException):
    """Exception raised when the partial charges on a molecule do not sum up to its formal charge."""


class DuplicateParameterError(OpenFFToolkitException):
    """Exception raised when trying to add a ParameterType that already exists"""


class ParameterLookupError(OpenFFToolkitException):
    """Exception raised when something goes wrong in a parameter lookup in
    ParameterHandler.__getitem__"""


class DuplicateVirtualSiteTypeException(OpenFFToolkitException):
    """Exception raised when trying to register two different virtual site classes with the same 'type'"""


class CallbackRegistrationError(OpenFFToolkitException, TypeError):
    """Error raised when callback registration fails."""


class MissingIndexedAttributeError(
    OpenFFToolkitException, IndexError, KeyError, AttributeError
):
    """Error raised when an indexed attribute does not exist"""
