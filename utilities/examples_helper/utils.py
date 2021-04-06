import click
import sys


def flatten(iterable):
    return (item for inner in iterable for item in inner)


def echo(*objects, sep=" ", end="\n", file=sys.stdout, color=None):
    """click.echo(), but like print()"""
    message = sep.join((str(object) for object in objects)) + end

    return click.echo(message, file=file, nl=False, err=False, color=color)


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
