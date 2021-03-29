import click
from click import echo

from pathlib import Path
import shutil

EXAMPLES_DIR = (Path(__file__) / "../../../examples").resolve()


class Example:
    """
    References to the files that make up an example
    """

    def __init__(self, path):
        """
        Initialize an `Example` from a `path`

        `path` should be the path to a readable directory including
        exactly one Jupyter notebook (*.ipynb). The path to this
        notebook is available as `Example.notebook`. A human and
        machine readable name of the `Example` is taken from the
        directory's path and is available as `Example.name`
        """
        self.path = Path(path).resolve()

        if self.path.is_dir():
            # `list(path.glob("*.ipynb"))` doesn't work here
            notebooks = [p for p in self.path.glob("*.ipynb")]

            if len(notebooks) == 1:
                self.notebook = notebooks[0]
            else:
                raise ValueError(
                    "path should be a folder with exactly one notebook"
                )
        else:
            raise ValueError(
                "path should be a folder with exactly one notebook"
            )

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
            ctx
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
@click.option("-t", "--target", type=click.Path(), default=Path('.'),
              help="Target directory in which to install the example.")
@click.option("-d", "--dry-run", default=False, is_flag=True,
              help="Describe what we're doing, but take no action.")
def install(example, target, dry_run):
    """Copy an example to a convenient location."""
    target = (target / example.name).resolve()

    if dry_run:
        echo("Dry run; taking no file system actions")

    echo(f"Sourcing '{example.name}' from {example.path}")
    echo(f"Installing to {target}")

    if not dry_run:
        shutil.copytree(example.path, target)


if __name__ == "__main__":
    main()
