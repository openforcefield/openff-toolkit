import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "smarty",
    version = "0.1.0",
    author = "John Chodera",
    author_email = "john.chodera@choderalab.org",
    description = ("Automated Bayesian atomtype sampling"),
    license = "GNU Lesser General Public License (LGPL), Version 3",
    keywords = "Bayesian atomtype sampling forcefield parameterization",
    url = "http://github.com/open-forcefield-group/smarty",
    packages=['smarty', 'smarty/tests'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU Lesser General Public License (LGPL), Version 3",
    ],
    entry_points={'console_scripts': ['smarty = smarty.cli:main']}),
)
