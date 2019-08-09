#!/usr/bin/env python

#======================================================================
# MODULE DOCSTRING
#======================================================================

"""
Test that the examples in the repo run without errors.
"""

#======================================================================
# GLOBAL IMPORTS
#======================================================================

import glob
import os
import re
import subprocess
import textwrap

import pytest

from openforcefield.utils import temporary_directory


#======================================================================
# TEST UTILITIES
#======================================================================

ROOT_DIR_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..')


def run_script_file(file_path):
    """Run through the shell a python script."""
    cmd = ['python', file_path]
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        raise Exception('Example {file_path} failed'.format(file_path=file_path))


def run_script_str(script_str):
    """Execute a Python string through the shell in a temporary directory.

    With respect to eval, this has the advantage of catching all import errors.

    """
    with temporary_directory() as tmp_dir:
        temp_file_path = os.path.join(tmp_dir, 'temp.py')
        # Create temporary python script.
        with open(temp_file_path, 'w') as f:
            f.write(script_str)
        # Run the Python script.
        try:
            run_script_file(temp_file_path)
        except:
            script_str = textwrap.indent(script_str, '    ')
            raise Exception(f'The following script failed:\n{script_str}')


def find_examples():
    """Find all examples in the examples folder.

    Returns
    -------
    example_file_paths : List[str]
        List of python script to execute.
    """
    slow_examples = {
        os.path.join('SMIRNOFF_comparison', 'compare_set_energies.py')
    }
    examples_dir_path = os.path.join(ROOT_DIR_PATH, 'examples')

    example_file_paths = []
    for example_file_path in glob.glob(os.path.join(examples_dir_path, '*', '*.py')):
        example_file_path = os.path.relpath(example_file_path)
        # Check if this is a slow test.
        for slow_example in slow_examples:
            if slow_example in example_file_path:
                # This is a slow example.
                example_file_path = pytest.param(example_file_path, marks=pytest.mark.slow)
        example_file_paths.append(example_file_path)
    return example_file_paths


def find_readme_examples():
    """Yield the Python scripts in the main README.md file.

    Returns
    -------
    readme_examples : List[str]
        The list of Python scripts included in the README.md files.
    """
    readme_file_path = os.path.join(ROOT_DIR_PATH, 'README.md')
    with open(readme_file_path, 'r') as f:
        readme_content = f.read()
    return re.findall('```python(.*?)```', readme_content, flags=re.DOTALL)


def find_readme_links():
    """Yield all the links in the main README.md file.

    Returns
    -------
    readme_examples : List[str]
        The list of links included in the README.md file.
    """
    readme_file_path = os.path.join(ROOT_DIR_PATH, 'README.md')
    with open(readme_file_path, 'r') as f:
        readme_content = f.read()
    return re.findall('http[s]?://(?:[0-9a-zA-Z]|[-/.%:_])+', readme_content)


#======================================================================
# TESTS
#======================================================================

@pytest.mark.parametrize('readme_example_str', find_readme_examples())
def test_readme_examples(readme_example_str):
    """Test the example"""
    run_script_str(readme_example_str)


@pytest.mark.parametrize('readme_link', find_readme_links())
def test_readme_links(readme_link):
    """Test the example"""
    from urllib.request import Request, urlopen
    # Some websites do not accept requests that don't specify the
    # client and the type of accepted documents so we add fake info
    # to avoid the response being an error.
    headers = {'User-Agent':'Mozilla/5.0',
               'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',}
    request = Request(readme_link, headers=headers)
    urlopen(request)


@pytest.mark.parametrize('example_file_path', find_examples())
def test_examples(example_file_path):
    """Test that the example run without errors."""
    run_script_file(example_file_path)

