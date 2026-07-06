import functools
import json
import os
import subprocess
import warnings


@functools.lru_cache
def _get_conda_list_package_versions() -> dict[str, str]:
    """
    Returns the versions of any packages found while executing `conda list`.
    If no conda executable is found, emits CondaExecutableNotFoundWarning
    """
    from openff.toolkit._utilities.warnings import CondaExecutableNotFoundWarning

    if os.environ.get("PIXI_IN_SHELL") == "1" and os.environ.get("PIXI_EXE"):
        conda_command = "{} list --json --manifest-path {}".format(
            os.environ["PIXI_EXE"], os.environ["PIXI_PROJECT_MANIFEST"]
        )
    elif os.environ.get("CONDA_SHLVL", "0") != "0" and os.environ.get("CONDA_EXE"):
        conda_command = "{} list --json".format(os.environ["CONDA_EXE"])
    elif os.environ.get("CONDA_SHLVL", "0") != "0" and os.environ.get("MAMBA_EXE"):
        conda_command = "{} list --json".format(os.environ["MAMBA_EXE"])
    else:
        warnings.warn(
            "No conda/mamba/micromamba executable found. Unable to determine package versions.",
            CondaExecutableNotFoundWarning,
        )
        return dict()

    output = json.loads(subprocess.check_output(conda_command.split()).decode())

    package_versions = {}
    for _package in output:
        package_versions[_package["name"]] = _package["version"]

    return package_versions


def get_ambertools_version() -> str | None:
    """
    Attempts to retrieve the version of the currently installed AmberTools.

    There are two soft failure modes, each of which cause this function to return `None`:
        1. If there is a failure calling `{conda|mamba|etc.} list`, the user is
            warned by `_get_conda_list_package_versions` and this function returns `None`.
        2. If there is a failure calling `{conda|mamba|etc.} list`, this function
            still returns `None`, but without a warning associated with the above failure.
    """

    try:
        return _get_conda_list_package_versions().get("ambertools", None)
    except (
        ValueError,  # Issue 98
        subprocess.CalledProcessError,  # Issue 101
    ):
        from openff.toolkit._utilities.warnings import CondaExecutableNotFoundWarning

        warnings.warn(
            "Something went wrong parsing the output of `conda list` or similar. Unable to "
            "determine AmberTools version, returning None.",
            CondaExecutableNotFoundWarning,
        )

        return None
