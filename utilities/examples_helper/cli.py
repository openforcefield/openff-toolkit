import shutil
import subprocess as sp
import sys
from pathlib import Path

import click
import conda.cli.python_api
import conda.exceptions

OFFTK_ROOT = (Path(__file__) / "../../../").resolve()
EXAMPLES_DIR = OFFTK_ROOT / "examples"
EXAMPLES_ENV = OFFTK_ROOT / "devtools/conda-envs/examples.yaml"


def flatten(iterable):
    return (item for inner in iterable for item in inner)


def echo(*objects, sep=" ", end="\n", file=sys.stdout, color=None):
    """click.echo(), but like print()"""
    message = sep.join((str(object) for object in objects)) + end

    return click.echo(message, file=file, nl=False, err=False, color=color)


def conda_env(cmd, *args, **kwargs):
    """Run `conda env`; mimics signature of `conda.cli.python_api.run_command()`

    This function is necessary because `conda env` run through run_command()
    incorrectly consumes arguments from the parent process, rather than from
    the run_command() call.
    """
    kwargs = dict(stdin="STRING", stdout="STRING", use_exception_handler=False)
    kwargs.update(kwargs)

    if "search_path" in kwargs:
        raise ValueError("conda env does not support search_path")

    if kwargs["stdin"] == "STRING":
        del kwargs["stdin"]
        kwargs["text"] = True
    if kwargs["stdout"] == "STRING":
        del kwargs["stdout"]
        kwargs["text"] = True
    if not kwargs["use_exception_handler"]:
        del kwargs["use_exception_handler"]
        kwargs["check"] = True

    proc = sp.run(["conda", "env", cmd, *args], **kwargs)

    return (proc.stdout, proc.stderr, proc.returncode)


def conda_cmd(cmd, *args, **kwargs):
    """Wrapper around both run_command() and conda_env()"""
    args = [str(arg) for arg in args]
    if cmd == "env":
        if "--dry-run" in args:
            return (None, None, 0)
        return conda_env(*args, **kwargs)
    else:
        try:
            return conda.cli.python_api.run_command(cmd, *args, **kwargs)
        except conda.exceptions.DryRunExit:
            return (None, None, 0)


def get_current_conda_prefix():
    """Return the currently active conda environment

    Hacky as hell, but Conda doesn't seem to provide a robust solution:
    https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#determining-your-current-environment
    """
    stdout, stderr, code = conda_cmd("info", "--envs")
    for line in stdout.splitlines():
        if line.startswith("#"):
            continue
        name, _, prefix = line.partition(" * ")
        if prefix:
            return Path(prefix.strip())
    raise ValueError("Could not determine prefix")


class Example:
    """
    References to the files that make up an example
    """

    def __init__(self, path):
        """
        Initialize an `Example` from a `path`

        `path` should be the path to a readable directory including
        exactly one Jupyter notebook (*.ipynb). The path to this
        notebook, relative to the example directory, is available as
        `Example.notebook`. A human and machine readable name of the
        `Example` is taken from the directory's path and is available
        as `Example.name`. If a file named environment.yml or
        environment.yaml exists in the directory, it is available as
        `Example.environment`
        """
        self.path = Path(path).resolve()

        if not self.path.is_dir():
            raise ValueError("path should be a directory")

        # `list(path.glob("*.ipynb"))` doesn't work here
        notebooks = [p for p in self.path.glob("*.ipynb")]

        if len(notebooks) == 1:
            self.notebook = notebooks[0].relative_to(self.path)
        else:
            raise ValueError("path should be a folder with exactly one notebook")

        if (self.path / "environment.yaml").is_file():
            self.environment = self.path / "environment.yaml"
        elif (self.path / "environment.yml").is_file():
            self.environment = self.path / "environment.yml"
        else:
            self.environment = None

    def __str__(self):
        return f"Example('{self.path}')"

    @property
    def name(self):
        return self.path.name


class ExampleArg(click.ParamType):
    """Construct an `Example` from a click argument"""

    def convert(self, value, param, ctx):
        """
        Takes a path or name and looks in `EXAMPLES_DIR` or the working
        directory for the corresponding example.
        """

        if value is None:
            return None
        if isinstance(value, Example):
            return value

        value = Path(value)

        try:
            return Example(EXAMPLES_DIR / value)
        except ValueError:
            pass

        try:
            return Example(EXAMPLES_DIR / value)
        except ValueError:
            pass

        self.fail(
            f"Could not find example {value}. It should be a folder "
            f"containing exactly one Jupyter notebook.",
            param,
            ctx,
        )


def style(value, *args, **kwargs):
    return click.style(str(value), *args, **kwargs)


def strong(value):
    return style(value, bold=True)


def style_path(value):
    return strong(value)


def style_cmd(value):
    return (
        "\n"
        + style(" $ ", dim=True, bg="black")
        + style(value, bold=True, bg="black")
        + "\n"
    )


@click.group()
def main():
    """Easily run examples with the Open Force Field Toolkit."""
    pass


@main.command()
def list():
    """List the available examples."""
    for example in EXAMPLES_DIR.iterdir():
        try:
            example = Example(example)
        except ValueError:
            continue

        echo(example.name)


@main.command()
@click.argument("example", type=ExampleArg())
@click.option(
    "-t",
    "--target",
    type=click.Path(),
    default=Path("."),
    help="Target directory in which to install the example.",
)
@click.option(
    "-d",
    "--dry-run",
    default=False,
    is_flag=True,
    help="Describe what we're doing, but take no action.",
)
@click.option(
    "-q", "--quiet", default=False, is_flag=True, help="Produce less text output."
)
@click.option(
    "-u",
    "--update-current",
    default=False,
    is_flag=True,
    help="Update the current Conda environment, rather than creating a new one.",
)
def install(example, target, dry_run, quiet, update_current):
    """Copy an example to a convenient location."""
    target = target / example.name

    if dry_run:
        echo(style("Dry run; taking no file system actions.", fg="red"))

    if not quiet:
        echo(
            f"Taking example {style_path(example.name)} from {style_path(example.path)}"
        )
        echo(f"Installing to {style_path(target)}")

    if target.exists():
        raise click.BadParameter(
            f"Directory {style_path(target)} already exists", param_hint="--target"
        )
    if not dry_run:
        shutil.copytree(example.path, target)

    conda_args = []
    if dry_run:
        conda_args.append("--dry-run")
    if quiet:
        conda_args.append("--quiet")

    try:
        if not update_current:
            prefix = target / ".env"
            if not quiet:
                echo(f"Creating new conda environment at prefix {style_path(prefix)}")

            conda_cmd("create", "--prefix", prefix, *conda_args)
        else:
            prefix = get_current_conda_prefix()

        environments = [EXAMPLES_ENV, example.environment]
        environments = [("--file", e) for e in environments if e is not None]
        if not quiet:
            echo(
                f"Installing dependencies from the following files into environment at {style_path(prefix)}:"
            )
            for (_, path) in environments:
                echo(f"    {style_path(path)}")
        environments = flatten(environments)

        conda_cmd(
            "env",
            "update",
            "--prefix",
            prefix,
            *environments,
            *conda_args,
            # Pass through stdout and stderr
            stdout=None,
            stderr=None,
        )

        if dry_run:
            echo(style("Dry run complete.", fg="red"))
        else:
            echo(
                style("Installation complete.", fg="green"),
                f"Activate the example Conda environment with",
                style_cmd(f"conda activate {prefix.resolve()}"),
                f"With the environment active and working from the new "
                f"{style_path(target)} directory, run the example with",
                style_cmd(f"jupyter notebook {example.notebook}"),
                f"To get started right now:",
                style_cmd(
                    f"conda run -p {prefix} --cwd {target} jupyter notebook {example.notebook}"
                ),
                sep="\n",
            )
    except Exception as e:
        if not dry_run:
            echo(
                style(
                    "Exception raised; cleaning up installation in progress", fg="red"
                ),
                file=sys.stderr,
            )
            shutil.rmtree(target)
        raise e


if __name__ == "__main__":
    main()
