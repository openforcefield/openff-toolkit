"""
openff-test-parameter-plugins
A test package used to ensure that parameterhandler plugins are handled correctly
"""
from setuptools import setup

setup(
    name="openff-test-parameter-plugins",
    packages=["custom_plugins"],
    include_package_data=True,
    entry_points={
        "openff.toolkit.plugins.handlers": [
            "CustomHandler = custom_plugins.handler_plugins:CustomHandler",
            "WrongSubclass = custom_plugins.handler_plugins:WrongSubclass",
        ]
    },
)
