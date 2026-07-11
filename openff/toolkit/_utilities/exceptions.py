class OpenFFError(BaseException):
    """The base exception from which most custom exceptions in openff-utilities inherited."""


class MissingOptionalDependencyError(OpenFFError, ImportError):
    """An exception raised when an optional dependency is required
    but cannot be found.

    Attributes
    ----------
    library_name
        The name of the missing library.
    license_issue
        Whether the library was importable but was unusable due
        to a missing license.
    """

    def __init__(self, library_name: str, license_issue: bool = False):
        """

        Parameters
        ----------
        library_name
            The name of the missing library.
        license_issue
            Whether the library was importable but was unusable due
            to a missing license.
        """

        message = f"The required {library_name} module could not be imported."

        if license_issue:
            message = f"{message} This is due to a missing license."

        library_name_corrected = library_name.replace(".", "-")

        if "openeye" not in library_name:
            message = (
                f"{message} Try installing the package by running "
                f"`conda install -c conda-forge {library_name_corrected}`"
            )

        super().__init__(message)

        self.library_name = library_name
        self.license_issue = license_issue


class CondaExecutableNotFoundError(OpenFFError):
    """
    A conda (or mamba/micromamba) executable is not found.
    """
