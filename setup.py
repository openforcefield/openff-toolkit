import os
from os.path import relpath, join
from setuptools import setup
import versioneer

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
    author = "Open Forcefield Consortium",
    author_email = "john.chodera@choderalab.org",
    description = ("Open Forcefield Toolkit"),
    license = "MIT",
    keywords = "molecular mechanics, forcefield, Bayesian parameterization",
    url = "http://github.com/openforcefield/openforcefield",
    packages=[
        'openforcefield',
        'openforcefield/tests',
        'openforcefield/data',
        'openforcefield/topology',
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
        'Programming Language :: Python :: 3',
    ],
    entry_points={'console_scripts': []},
    package_data={'openforcefield': find_package_data('openforcefield/data', 'openforcefield')},
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
