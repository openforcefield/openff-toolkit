class MessageException(Exception):
    """DEPRECATED: A base class for exceptions that print out a string given in their constructor"""

    import warnings

    # Must "turn on" DeprecationWawrning, see https://stackoverflow.com/a/60486394/4248961
    warnings.filterwarnings("default", category=DeprecationWarning, module=__name__)

    warnings.warn(
        "MessageException is deprecated and will be removed in a future release of the OpenFF Toolkit."
        "All custom exceptions now inherit from OpenFFToolkitException, which should be used as a "
        "replacement for MessageException. (`from openff.toolkit.utils.exceptions import OpenFFToolkitException`",
        DeprecationWarning,
    )


class OpenFFToolkitException(Exception):
    """Base exception for custom exceptions raised by the OpenFF Toolkit"""


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
    """An unhandled error occured in an external toolkit during charge calculation"""


class InvalidIUPACNameError(OpenFFToolkitException):
    """Failed to parse IUPAC name"""


class AntechamberNotFoundError(OpenFFToolkitException):
    """The antechamber executable was not found"""


class ParseError(ValueError):
    """The record couple not be parsed into the given format"""


class InvalidConformerError(Exception):
    """
    This error is raised when the conformer added to the molecule
    has a different connectivity to that already defined.
    or anyother conformer related issues.
    """


class SmilesParsingError(Exception):
    """
    This error is rasied when parsing a smiles string results in an error.
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
    Exception for a when unique_molecules is required but not found
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
    Exception thrown when an incompatible SMIRNOFF version data structure in attempted to be read.
    """


class SMIRNOFFAromaticityError(OpenFFToolkitException):
    """
    Exception thrown when an incompatible SMIRNOFF aromaticity model is checked for compatibility.
    """


class ParseError(OpenFFToolkitException):
    """
    Error for when a SMIRNOFF data structure is not parseable by a ForceField
    """


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


class UnassignedValenceParameterException(Exception):
    """Exception raised when there are valence terms for which a ParameterHandler can't find parameters."""


class UnassignedBondParameterException(UnassignedValenceParameterException):
    """Exception raised when there are bond terms for which a ParameterHandler can't find parameters."""


class UnassignedAngleParameterException(UnassignedValenceParameterException):
    """Exception raised when there are angle terms for which a ParameterHandler can't find parameters."""


class UnassignedProperTorsionParameterException(UnassignedValenceParameterException):
    """Exception raised when there are proper torsion terms for which a ParameterHandler can't find parameters."""


class UnassignedMoleculeChargeException(Exception):
    """Exception raised when no charge method is able to assign charges to a molecule."""


class NonintegralMoleculeChargeException(Exception):
    """Exception raised when the partial charges on a molecule do not sum up to its formal charge."""


class DuplicateParameterError(OpenFFToolkitException):
    """Exception raised when trying to add a ParameterType that already exists"""


class ParameterLookupError(OpenFFToolkitException):
    """Exception raised when something goes wrong in a parameter lookup in
    ParameterHandler.__getitem__"""


class DuplicateVirtualSiteTypeException(Exception):
    """Exception raised when trying to register two different virtual site classes with the same 'type'"""


class CallbackRegistrationError(TypeError):
    """Error raised when callback registration fails."""
