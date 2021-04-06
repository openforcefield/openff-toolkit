import shutil
import sys
from pathlib import Path

import click

# Make sure the examples_helper dir is in PYTHONPATH, even if run through setup.py
sys.path.insert(0, str((Path(__file__) / "..").resolve()))

from example import Example, ExampleArg, OFFTK_ROOT, EXAMPLES_DIR, EXAMPLES_ENV
from utils import *
from pathlib import Path
import conda



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
    "--update-environment",
    default=False,
    is_flag=True,
    help="Update a Conda environment with the example's dependencies.",
)
@click.option(
    "-e",
    "--create-environment",
    default=False,
    is_flag=True,
    help="Create a new Conda environment with the example's dependencies.",
)
@click.option(
    "-p",
    "--prefix",
    type=click.Path(),
    default=None,
    help="Prefix directory for the Conda environment to update.",
)
@click.option(
    "-n",
    "--name",
    help="Name of the Conda environment to update.",
)
def install(example, target, dry_run, quiet, update_environment, create_environment, prefix, name):
    """Copy EXAMPLE to a convenient location."""
    target = target / example.name

    if dry_run:
        echo(style("Dry run; taking no file system actions.", fg="red"))

    copy_example(example, target, dry_run=dry_run, quiet=quiet)

    if update_environment and create_environment:
        raise click.BadParameter(
            "Options --update-environment and --create-environment are mutually exclusive; specify one or the other but not both",
            param_hint="--update-environment",
        )
    if prefix and name:
        raise click.BadParameter(
            "Options --prefix and --name are mutually exclusive; specify one or the other but not both",
            param_hint="--name",
        )
    if update_environment and prefix is None and name is None:
        prefix = conda.get_current_prefix()
    if create_environment and prefix is None and name is None:
        prefix = target / ".env"

    try:
        if create_environment:
            create_env(prefix=prefix, name=name, dry_run=dry_run, quiet=quiet)
        if update_environment or create_environment:
            update_env(example, prefix=prefix, name=name, dry_run=dry_run, quiet=quiet)

        if dry_run:
            echo(style("Dry run complete.", fg="red"))
        else:
            echo(style("Installation complete.", fg="green"))
            if (update_environment or create_environment) and (prefix or name):
                echo(
                    f"Activate the example Conda environment:",
                )
                if prefix:
                    echo(
                        style_cmd(f"conda activate {prefix.resolve()}")
                    )
                else:
                    echo(
                        style_cmd(f"conda activate {name}")
                    )

            echo(
                f"From the new {style_path(target)} directory, run the example:",
                style_cmd(f"jupyter notebook {example.notebook}"),
            )
            if (update_environment or create_environment) and (prefix or name):
                echo("To get started right now:")
                if prefix:
                    echo(style_cmd(
                        f"conda run -p {prefix} --cwd {target} jupyter notebook {example.notebook}"
                    ))
                else:
                    echo(style_cmd(
                        f"conda run -n {name} --cwd {target} jupyter notebook {example.notebook}"
                    ))
            else:
                echo(f"If you need to install dependencies, try the --update-environment or --create-environment switches to this script")
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


def copy_example(example, target, dry_run=False, quiet=False):
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


def create_env(prefix=None, name=None, dry_run=False, quiet=False):
    if prefix and name:
        raise ValueError("prefix and name are mutually exclusive")
    if not quiet:
        if prefix:
            echo(f"Creating new conda environment at prefix {style_path(prefix)}")
        elif name:
            echo(f"Creating new conda environment with name {style_path(name)}")

    conda.create_environment(prefix=prefix, name=name, dry_run=dry_run, quiet=quiet)


def update_env(example, prefix=None, name=None, dry_run=False, quiet=False):
    if prefix and name:
        raise ValueError("prefix and name are mutually exclusive")

    environments = [EXAMPLES_ENV, example.environment]
    environments = [("--file", e) for e in environments if e is not None]
    if not quiet:
        echo(
            f"Installing dependencies from the following files into environment at {style_path(prefix)}:"
        )
        for (_, path) in environments:
            echo(f"    {style_path(path)}")
    environments = flatten(environments)

    if not dry_run:
        conda.update_envs(environments, prefix=prefix, name=name)


if __name__ == "__main__":
    main()
