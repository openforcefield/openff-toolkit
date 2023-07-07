class MissingUnitError(ValueError):
    """
    Exception for data missing a unit tag.
    """


class UnitValidationError(ValueError):
    """
    Exception for bad behavior when validating unit-tagged data.
    """


class UnsupportedExportError(BaseException):
    """
    Exception for attempting to write to an unsupported file format.
    """
