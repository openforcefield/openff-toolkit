"""
openff-toolkit
A modern, extensible library for molecular mechanics force field science from the Open Force Field Consortium.
"""
from setuptools import setup, find_namespace_packages
import versioneer

short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except IOError:
    long_description = "\n".join(short_description[2:]),


setup(
    name="openff-toolkit",

    author="Open Forcefield Consortium",
    author_email="info@openforcefield.org",

    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",

    keywords="molecular mechanics, forcefield, parameterization",
    url="http://github.com/openforcefield/openff-toolkit",

    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),

    license="MIT",

    packages=find_namespace_packages(include=['openff.*']),
    include_package_data=True,

    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT",
        'Programming Language :: Python :: 3',
    ],
    entry_points={'console_scripts': []},
)
