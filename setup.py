import os
from os.path import relpath, join
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files

setup(
    name = "openforcefield",
    version = "0.0.3",
    author = "John Chodera, David Mobley, and others",
    author_email = "john.chodera@choderalab.org",
    description = ("Open Force Field Group tools"),
    license = "MIT",
    keywords = "molecular mechanics, forcefield, Bayesian parameterization",
    url = "http://github.com/open-forcefield-group/openforcefield",
    packages=[
        'openforcefield',
        'openforcefield/tests',
        'openforcefield/data',
        'openforcefield/typing',
        'openforcefield/typing/chemistry',
        'openforcefield/typing/engines',
        'openforcefield/typing/engines/smirnoff',
        'openforcefield/utils',
        ],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT",
    ],
    entry_points={'console_scripts': []},
    package_data={'openforcefield': find_package_data('openforcefield/data', 'openforcefield')},
)
