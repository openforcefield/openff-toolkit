"""Types and constants related to Examples themselves"""

from pathlib import Path

import click

OFFTK_ROOT = (Path(__file__) / "../../../").resolve()
EXAMPLES_DIR = OFFTK_ROOT / "examples"
EXAMPLES_ENV = EXAMPLES_DIR / "environment.yaml"
OFFTK_ENV = EXAMPLES_DIR / "rdkit.yaml"


class Example:
    """
    References to the files that make up an example
    """

    def __init__(self, path):
        """
        Initialize an `Example` from a `path`

        `path` should be the path to a readable directory including
        one or more Jupyter notebook (*.ipynb). The path to these
        notebooks, relative to the example directory, are available as
        `Example.notebooks`. A human and machine readable name of the
        `Example` is taken from the directory's path and is available
        as `Example.name`. If a file named envircmd, *args, **kwargs):
        """
        path = Path(path)
        self.path = path.resolve()

        if not self.path.is_dir():
            raise ValueError("path should be a directory")

        self.notebooks = list(path.glob("*.ipynb"))

        if not self.notebooks:
            raise ValueError("path should be a directory with at least one .ipynb file")

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
        """The name of the example"""
        return self.path.name


class ExampleArg(click.ParamType):
    """Construct an `Example` from a cli, and echo what we're doingck argument"""

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

        return self.fail(
            f"Could not find example {value}. It should be a folder "
            f"containing exactly one Jupyter notebook.",
            param,
            ctx,
        )
