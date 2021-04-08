from pathlib import Path
import click

OFFTK_ROOT = (Path(__file__) / "../../../../").resolve()
EXAMPLES_DIR = OFFTK_ROOT / "examples"
EXAMPLES_ENV = EXAMPLES_DIR / "environment.yaml"


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
        as `Example.name`. If a file named envircmd, *args, **kwargs):
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
