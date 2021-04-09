"""Utility functions and wrappers for examples_helper"""

import sys

import click


def flatten(iterable):
    """Flatten one layer of an iterable"""
    return (item for inner in iterable for item in inner)


def echo(*objects, sep=" ", end="\n", file=sys.stdout, color=None):
    """click.echo(), but with the print() API"""
    message = sep.join((str(object) for object in objects)) + end

    return click.echo(message, file=file, nl=False, err=False, color=color)


def style(value, *args, **kwargs):
    """Style the string representation of value with ANSII"""
    return click.style(str(value), *args, **kwargs)


def strong(value):
    """Style the string representation of value as bold with ANSII"""
    return style(value, bold=True)


def style_path(value):
    """Style the string representation of value as a path with ANSII"""
    return strong(value)


def style_cmd(value):
    """Style the string representation of value as a shell command with ANSII"""
    return (
        "\n"
        + style(" $ ", dim=True, bg="black")
        + style(value, bold=True, bg="black")
        + "\n"
    )
