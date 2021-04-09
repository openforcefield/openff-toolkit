"""Python wrapper around the mamba and conda cli tools"""

import subprocess as sp
from pathlib import Path

try:
    if sp.run("mamba", check=True, capture_output=True).returncode == 0:
        CMD = "mamba"
    else:
        CMD = "conda"
except FileNotFoundError:
    CMD = "conda"


def cmd(*args, **kwargs):
    """Run conda or mamba with the given commands and arguments"""
    run_kwargs = dict(
        text=True,
        check=True
    )
    run_kwargs.update(kwargs)

    return sp.run([CMD, *args], **run_kwargs)  # pylint: disable=subprocess-run-check


def create_environment(prefix=None, name=None, dry_run=False, quiet=False):
    """Create a new conda environment with the given prefix or name"""
    conda_args = []
    if dry_run:
        conda_args.append("--dry-run")
    if quiet:
        conda_args.append("--quiet")

    if prefix and name:
        raise ValueError("prefix and name are mutually exclusive")
    if not (prefix or name):
        raise ValueError("prefix or name is required")
    if prefix:
        conda_args.append("--prefix")
        conda_args.append(str(prefix))
    if name:
        conda_args.append("--name")
        conda_args.append(str(name))

    cmd("create", *conda_args)


def update_envs(environments, prefix=None, name=None):
    """Update the conda environment given by prefix or name with the provided environment.yml files"""
    conda_args = []
    if prefix and name:
        raise ValueError("prefix and name are mutually exclusive")
    if prefix:
        conda_args.append("--prefix")
        conda_args.append(str(prefix))
    elif name:
        conda_args.append("--name")
        conda_args.append(str(name))

    cmd(
        "env",
        "update",
        *conda_args,
        *environments,
    )


def get_current_prefix():
    """Return the currently active conda environment

    Hacky as hell, but Conda doesn't seem to provide a robust solution:
    https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#determining-your-current-environment
    """
    proc = cmd("info", "--envs", capture_output=True)
    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            continue
        _name, _, prefix = line.partition(" * ")
        if prefix:
            return Path(prefix.strip())
    raise ValueError("Could not determine prefix")
