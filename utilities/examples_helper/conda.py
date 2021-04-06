import subprocess as sp
import click

if sp.run(["mamba", "--version"]).returncode == 0:
    CMD = "mamba"
else:
    CMD = "conda"


def cmd(*args, **kwargs):
    run_kwargs = dict(
        text=True,
        check=True
    )
    run_kwargs.update(kwargs)

    return sp.run([CMD, *args], **run_kwargs)


def create_environment(prefix=None, name=None, dry_run=False, quiet=False):
    conda_args = []
    if dry_run:
        conda_args.append("--dry-run")
    if quiet:
        conda_args.append("--quiet")

    if prefix and name:
        raise ValueError("prefix and name are mutually exclusive")
    elif not (prefix or name):
        raise ValueError("prefix or name is required")
    elif prefix:
        conda_args.append("--prefix")
        conda_args.append(str(prefix))
    elif name:
        conda_args.append("--name")
        conda_args.append(str(name))

    proc = cmd("create", *conda_args)


def update_envs(environments, prefix=None, name=None):
    conda_args = []
    if prefix and name:
        raise ValueError("prefix and name are mutually exclusive")
    elif prefix:
        conda_args.append("--prefix")
        conda_args.append(str(prefix))
    elif name:
        conda_args.append("--name")
        conda_args.append(str(name))

    proc = cmd(
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
    proc = conda("info", "--envs", capture_output=True)
    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            continue
        name, _, prefix = line.partition(" * ")
        if prefix:
            return Path(prefix.strip())
    raise ValueError("Could not determine prefix")
