__all__ = (
    "MessageException",
    "IncompatibleUnitError",
    "MissingDependencyError",
    "MissingPackageError",
    "ToolkitUnavailableException",
    "LicenseError",
    "InvalidToolkitError",
    "InvalidToolkitRegistryError",
    "UndefinedStereochemistryError",
    "GAFFAtomTypeWarning",
    "ChargeMethodUnavailableError",
    "IncorrectNumConformersError",
    "IncorrectNumConformersWarning",
    "ChargeCalculationError",
    "InvalidIUPACNameError",
    "AntechamberNotFoundError",
    "ParseError",
)

# =============================================================================================
# COMMON EXCEPTION TYPES
# =============================================================================================


class MessageException(Exception):
    """A base class for exceptions that print out a string given in their constructor"""

    def __init__(self, msg):
        super().__init__(self, msg)
        self.msg = msg

    def __str__(self):
        return self.msg


class IncompatibleUnitError(MessageException):
    """
    Exception for when a parameter is in the wrong units for a ParameterHandler's unit system
    """

    pass


class MissingDependencyError(MessageException):
    """
    Exception for when an optional dependency is needed but not installed

    """

    def __init__(self, package_name):
        self.msg = (
            f"Missing dependency {package_name}. Try installing it "
            f"with\n\n$ conda install {package_name} -c conda-forge"
        )

        super().__init__(self.msg)


class MissingPackageError(MessageException):
    """This function requires a package that is not installed."""


class ToolkitUnavailableException(MessageException):
    """The requested toolkit is unavailable."""

    # TODO: Allow toolkit to be specified and used in formatting/printing exception.


class LicenseError(ToolkitUnavailableException):
    """This function requires a license that cannot be found."""


class InvalidToolkitError(MessageException):
    """A non-toolkit object was received when a toolkit object was expected"""


class InvalidToolkitRegistryError(MessageException):
    """An object other than a ToolkitRegistry or toolkit wrapper was received"""


class UndefinedStereochemistryError(MessageException):
    """A molecule was attempted to be loaded with undefined stereochemistry"""


class GAFFAtomTypeWarning(RuntimeWarning):
    """A warning raised if a loaded mol2 file possibly uses GAFF atom types."""


class ChargeMethodUnavailableError(MessageException):
    """A toolkit does not support the requested partial_charge_method combination"""


class IncorrectNumConformersError(MessageException):
    """The requested partial_charge_method expects a different number of conformers than was provided"""


class IncorrectNumConformersWarning(Warning):
    """The requested partial_charge_method expects a different number of conformers than was provided"""


class ChargeCalculationError(MessageException):
    """An unhandled error occured in an external toolkit during charge calculation"""


class InvalidIUPACNameError(MessageException):
    """Failed to parse IUPAC name"""


class AntechamberNotFoundError(MessageException):
    """The antechamber executable was not found"""
    
class ParseError(MessageException, ValueError):
    """The record couple not be parsed into the given format"""
