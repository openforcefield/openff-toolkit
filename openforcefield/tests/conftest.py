#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Configuration file for pytest.

This adds the following command line options.
- runslow: Run tests marked as slow (default is False).

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import logging

import numpy as np
import pytest
from simtk import unit

logger = logging.getLogger(__name__)


def pytest_configure(config):
    """
    Initialization hook to register custom markers without a pytest.ini
    More info: https://docs.pytest.org/en/latest/reference.html#initialization-hooks
    """
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with `-m 'not slow'`)"
    )


# =============================================================================================
# UTILITY FUNCTIONS
# =============================================================================================


def untar_full_alkethoh_and_freesolv_set():
    """When running slow tests, we unpack the full AlkEthOH and FreeSolv test
    sets in advance to speed things up.

    See
        test_forcefield.py::test_alkethoh_parameters_assignment
        test_forcefield.py::test_freesolv_parameters_assignment
    """
    import os
    import tarfile

    from openforcefield.utils import get_data_file_path

    molecule_dir_path = get_data_file_path("molecules")
    for tarfile_name in ["AlkEthOH_tripos.tar.gz", "FreeSolv.tar.gz"]:
        tarfile_path = os.path.join(molecule_dir_path, tarfile_name)
        with tarfile.open(tarfile_path, "r:gz") as tar:
            tar.extractall(path=molecule_dir_path)


# =============================================================================================
# CONFIGURATION
# =============================================================================================


def pytest_addoption(parser):
    """Add the pytest command line option --runslow and --failwip.

    If --runslow is not given, tests marked with pytest.mark.slow are
    skipped.

    If --failwip is not given, tests marked with pytest.mark.wip are
    xfailed.

    Parameters
    ----------
    parser : argparsing.parser
        The parser used by pytest to process arguments

    Returns
    -------
    None

    """

    # Loaded pytest plugins define their own arguments, and in certain cases
    # our options use the same name. Although this can define two different
    # behaviors for the same argument, the two options below are fairly
    # unambiguous in their interpretation, so we skip any errors and allow
    # the testing to proceed

    try:
        parser.addoption(
            "--runslow", action="store_true", default=False, help="run slow tests"
        )
    except ValueError:
        logger.warning(
            "Option --runslow already added elsewhere (from a plugin possibly?). Skipping..."
        )

    try:
        parser.addoption(
            "--failwip",
            action="store_true",
            default=False,
            help="fail work in progress tests",
        )
    except ValueError:
        logger.warning(
            "Option --failwip already added elsewhere (from a plugin possibly?). Skipping..."
        )


def pytest_collection_modifyitems(config, items):

    if config.getoption("runslow"):
        # If --runslow is given, we don't have to mark items for skipping,
        # but we need to extract the whole AlkEthOH and FreeSolv sets (see
        # test_forcefield::test_alkethoh/freesolv_parameters_assignment).
        untar_full_alkethoh_and_freesolv_set()
    else:
        # Mark for skipping all items marked as slow.
        skip_slow = pytest.mark.skip(
            reason="specify --runslow pytest option to run this test."
        )
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)

    # Mark work-in-progress tests for xfail.
    if not config.getoption("failwip"):
        xfail_wip_reason = (
            "This is a work in progress test. Specify "
            "--failwip pytest option to make this test fail."
        )
        for item in items:
            if "wip" in item.keywords:
                # Augment original reason.
                reason = xfail_wip_reason + item.get_closest_marker("wip").kwargs.get(
                    "reason", ""
                )
                item.add_marker(pytest.mark.xfail(reason=reason))


@pytest.fixture(scope="function")
def ethane():
    """
    An ethane Molecule created from a SMILES string
    """
    from openforcefield.topology.molecule import Molecule

    return Molecule.from_smiles("CC")


@pytest.fixture(scope="function")
def ethene():
    """
    An ethene Molecule created from a SMILES string
    """
    from openforcefield.topology.molecule import Molecule

    return Molecule.from_smiles("C=C")


@pytest.fixture(scope="function")
def propane():
    """
    A propane Molecule created from a SMILES string
    """
    from openforcefield.topology.molecule import Molecule

    return Molecule.from_smiles("CCC")


@pytest.fixture(scope="function")
def ethanol():
    """
    Creates an openforcefield.topology.Molecule representation of
    ethanol without the use of a cheminformatics toolkit
    """
    from openforcefield.topology.molecule import Molecule

    # Create an ethanol molecule without using a toolkit

    ethanol = Molecule()
    ethanol.add_atom(6, 0, False)  # C0
    ethanol.add_atom(6, 0, False)  # C1
    ethanol.add_atom(8, 0, False)  # O2
    ethanol.add_atom(1, 0, False)  # H3
    ethanol.add_atom(1, 0, False)  # H4
    ethanol.add_atom(1, 0, False)  # H5
    ethanol.add_atom(1, 0, False)  # H6
    ethanol.add_atom(1, 0, False)  # H7
    ethanol.add_atom(1, 0, False)  # H8
    ethanol.add_bond(0, 1, 1, False, fractional_bond_order=1.33)  # C0 - C1
    ethanol.add_bond(1, 2, 1, False, fractional_bond_order=1.23)  # C1 - O2
    ethanol.add_bond(0, 3, 1, False, fractional_bond_order=1)  # C0 - H3
    ethanol.add_bond(0, 4, 1, False, fractional_bond_order=1)  # C0 - H4
    ethanol.add_bond(0, 5, 1, False, fractional_bond_order=1)  # C0 - H5
    ethanol.add_bond(1, 6, 1, False, fractional_bond_order=1)  # C1 - H6
    ethanol.add_bond(1, 7, 1, False, fractional_bond_order=1)  # C1 - H7
    ethanol.add_bond(2, 8, 1, False, fractional_bond_order=1)  # O2 - H8
    charges = unit.Quantity(
        np.array([-0.4, -0.3, -0.2, -0.1, 0.00001, 0.1, 0.2, 0.3, 0.4]),
        unit.elementary_charge,
    )
    ethanol.partial_charges = charges

    return ethanol


@pytest.fixture(scope="function")
def reversed_ethanol():
    """
    Creates an openforcefield.topology.Molecule representation of
    ethanol without the use of a cheminformatics toolkit. This function
    reverses the atom indexing of ethanol
    """
    from openforcefield.topology.molecule import Molecule

    # Create an ethanol molecule without using a toolkit
    ethanol = Molecule()
    ethanol.add_atom(1, 0, False)  # H0
    ethanol.add_atom(1, 0, False)  # H1
    ethanol.add_atom(1, 0, False)  # H2
    ethanol.add_atom(1, 0, False)  # H3
    ethanol.add_atom(1, 0, False)  # H4
    ethanol.add_atom(1, 0, False)  # H5
    ethanol.add_atom(8, 0, False)  # O6
    ethanol.add_atom(6, 0, False)  # C7
    ethanol.add_atom(6, 0, False)  # C8
    ethanol.add_bond(8, 7, 1, False, fractional_bond_order=1.33)  # C8 - C7
    ethanol.add_bond(7, 6, 1, False, fractional_bond_order=1.23)  # C7 - O6
    ethanol.add_bond(8, 5, 1, False, fractional_bond_order=1)  # C8 - H5
    ethanol.add_bond(8, 4, 1, False, fractional_bond_order=1)  # C8 - H4
    ethanol.add_bond(8, 3, 1, False, fractional_bond_order=1)  # C8 - H3
    ethanol.add_bond(7, 2, 1, False, fractional_bond_order=1)  # C7 - H2
    ethanol.add_bond(7, 1, 1, False, fractional_bond_order=1)  # C7 - H1
    ethanol.add_bond(6, 0, 1, False, fractional_bond_order=1)  # O6 - H0
    charges = unit.Quantity(
        np.array([0.4, 0.3, 0.2, 0.1, 0.00001, -0.1, -0.2, -0.3, -0.4]),
        unit.elementary_charge,
    )
    ethanol.partial_charges = charges
    return ethanol


@pytest.fixture(scope="function")
def benzene_no_aromatic():
    """
    Creates an openforcefield.topology.Molecule representation of benzene through the API with aromatic bonds
    not defied, used to test the levels of isomorphic matching.
    """
    from openforcefield.topology.molecule import Molecule

    benzene = Molecule()
    benzene.add_atom(6, 0, False)  # C0
    benzene.add_atom(6, 0, False)  # C1
    benzene.add_atom(6, 0, False)  # C2
    benzene.add_atom(6, 0, False)  # C3
    benzene.add_atom(6, 0, False)  # C4
    benzene.add_atom(6, 0, False)  # C5
    benzene.add_atom(1, 0, False)  # H6
    benzene.add_atom(1, 0, False)  # H7
    benzene.add_atom(1, 0, False)  # H8
    benzene.add_atom(1, 0, False)  # H9
    benzene.add_atom(1, 0, False)  # H10
    benzene.add_atom(1, 0, False)  # H11
    benzene.add_bond(0, 5, 1, False)  # C0 - C5
    benzene.add_bond(0, 1, 1, False)  # C0 - C1
    benzene.add_bond(1, 2, 1, False)  # C1 - C2
    benzene.add_bond(2, 3, 1, False)  # C2 - C3
    benzene.add_bond(3, 4, 1, False)  # C3 - C4
    benzene.add_bond(4, 5, 1, False)  # C4 - C5
    benzene.add_bond(0, 6, 1, False)  # C0 - H6
    benzene.add_bond(1, 7, 1, False)  # C1 - C7
    benzene.add_bond(2, 8, 1, False)  # C2 - C8
    benzene.add_bond(3, 9, 1, False)  # C3 - C9
    benzene.add_bond(4, 10, 1, False)  # C4 - C10
    benzene.add_bond(5, 11, 1, False)  # C5 - C11
    return benzene


@pytest.fixture(scope="function")
def acetaldehyde():
    """
    Creates an openforcefield.topology.Molecule representation of acetaldehyde through the API
    """
    from openforcefield.topology.molecule import Molecule

    acetaldehyde = Molecule()
    acetaldehyde.add_atom(6, 0, False)  # C0
    acetaldehyde.add_atom(6, 0, False)  # C1
    acetaldehyde.add_atom(8, 0, False)  # O2
    acetaldehyde.add_atom(1, 0, False)  # H3
    acetaldehyde.add_atom(1, 0, False)  # H4
    acetaldehyde.add_atom(1, 0, False)  # H5
    acetaldehyde.add_atom(1, 0, False)  # H6
    acetaldehyde.add_bond(0, 1, 1, False)  # C0 - C1
    acetaldehyde.add_bond(1, 2, 2, False)  # C1 = O2
    acetaldehyde.add_bond(0, 3, 1, False)  # C0 - H3
    acetaldehyde.add_bond(0, 4, 1, False)  # C0 - H4
    acetaldehyde.add_bond(0, 5, 1, False)  # C0 - H5
    acetaldehyde.add_bond(1, 6, 1, False)  # C1 - H6
    charges = unit.Quantity(
        np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), unit.elementary_charge
    )
    acetaldehyde.partial_charges = charges
    return acetaldehyde


@pytest.fixture(scope="function")
def water():
    """
    Creates an openforcefield.topology.Molecule representation of water through the API
    """
    from openforcefield.topology.molecule import Molecule

    mol = Molecule()
    mol.add_atom(1, 0, False)  # H1
    mol.add_atom(8, 0, False)  # O
    mol.add_atom(1, 0, False)  # H2
    mol.add_bond(0, 1, 1, False)  # H1 - O
    mol.add_bond(1, 2, 1, False)  # O - H2
    charges = unit.Quantity(np.array([0.0, 0.0, 0.0]), unit.elementary_charge)
    mol.partial_charges = charges
    return mol


@pytest.fixture(scope="function")
def ammonia():
    """
    Creates an openforcefield.topology.Molecule representation of ammonia through the API
    """
    from openforcefield.topology.molecule import Molecule

    mol = Molecule()
    mol.add_atom(1, 0, False)  # H1
    mol.add_atom(7, 0, False)  # N
    mol.add_atom(1, 0, False)  # H2
    mol.add_atom(1, 0, False)  # H3
    mol.add_bond(0, 1, 1, False)  # H1 - N
    mol.add_bond(1, 2, 1, False)  # N - H2
    mol.add_bond(1, 3, 1, False)  # N - H3
    charges = unit.Quantity(np.array([0.0, 0.0, 0.0, 0.0]), unit.elementary_charge)
    mol.partial_charges = charges
    return mol


@pytest.fixture(scope="function")
def acetate():
    """
    Creates an openforcefield.topology.Molecule representation of
    acetate without the use of a cheminformatics toolkit
    """
    from openforcefield.topology.molecule import Molecule

    # Create an acetate molecule without using a toolkit
    acetate = Molecule()
    acetate.add_atom(6, 0, False)  # C0
    acetate.add_atom(6, 0, False)  # C1
    acetate.add_atom(8, 0, False)  # O2
    acetate.add_atom(8, -1, False)  # O3
    acetate.add_atom(1, 0, False)  # H4
    acetate.add_atom(1, 0, False)  # H5
    acetate.add_atom(1, 0, False)  # H6
    acetate.add_bond(0, 1, 1, False)  # C0 - C1
    acetate.add_bond(1, 2, 2, False)  # C1 = O2
    acetate.add_bond(1, 3, 1, False)  # C1 - O3[-1]
    acetate.add_bond(0, 4, 1, False)  # C0 - H4
    acetate.add_bond(0, 5, 1, False)  # C0 - H5
    acetate.add_bond(0, 6, 1, False)  # C0 - H6
    return acetate


@pytest.fixture(scope="function")
def cyclohexane():
    """
    Creates an openforcefield.topology.Molecule representation of
    cyclohexane without the use of a cheminformatics toolkit
    """
    from openforcefield.topology.molecule import Molecule

    cyclohexane = Molecule()
    cyclohexane.add_atom(6, 0, False)  # C0
    cyclohexane.add_atom(6, 0, False)  # C1
    cyclohexane.add_atom(6, 0, False)  # C2
    cyclohexane.add_atom(6, 0, False)  # C3
    cyclohexane.add_atom(6, 0, False)  # C4
    cyclohexane.add_atom(6, 0, False)  # C5
    cyclohexane.add_atom(1, 0, False)  # H6
    cyclohexane.add_atom(1, 0, False)  # H7
    cyclohexane.add_atom(1, 0, False)  # H8
    cyclohexane.add_atom(1, 0, False)  # H9
    cyclohexane.add_atom(1, 0, False)  # H10
    cyclohexane.add_atom(1, 0, False)  # H11
    cyclohexane.add_atom(1, 0, False)  # H12
    cyclohexane.add_atom(1, 0, False)  # H13
    cyclohexane.add_atom(1, 0, False)  # H14
    cyclohexane.add_atom(1, 0, False)  # H15
    cyclohexane.add_atom(1, 0, False)  # H16
    cyclohexane.add_atom(1, 0, False)  # H17
    cyclohexane.add_bond(0, 1, 1, False)  # C0 - C1
    cyclohexane.add_bond(1, 2, 1, False)  # C1 - C2
    cyclohexane.add_bond(2, 3, 1, False)  # C2 - C3
    cyclohexane.add_bond(3, 4, 1, False)  # C3 - C4
    cyclohexane.add_bond(4, 5, 1, False)  # C4 - C5
    cyclohexane.add_bond(5, 0, 1, False)  # C5 - C0
    cyclohexane.add_bond(0, 6, 1, False)  # C0 - H6
    cyclohexane.add_bond(0, 7, 1, False)  # C0 - H7
    cyclohexane.add_bond(1, 8, 1, False)  # C1 - H8
    cyclohexane.add_bond(1, 9, 1, False)  # C1 - H9
    cyclohexane.add_bond(2, 10, 1, False)  # C2 - H10
    cyclohexane.add_bond(2, 11, 1, False)  # C2 - H11
    cyclohexane.add_bond(3, 12, 1, False)  # C3 - H12
    cyclohexane.add_bond(3, 13, 1, False)  # C3 - H13
    cyclohexane.add_bond(4, 14, 1, False)  # C4 - H14
    cyclohexane.add_bond(4, 15, 1, False)  # C4 - H15
    cyclohexane.add_bond(5, 16, 1, False)  # C5 - H16
    cyclohexane.add_bond(5, 17, 1, False)  # C5 - H17
    return cyclohexane


@pytest.fixture(scope="function")
def dioxygen():
    """
    Creates an openforcefield.topology.Molecule representation of
    dioxygen without the use of a cheminformatics toolkit
    """
    from openforcefield.topology.molecule import Molecule

    dioxygen = Molecule()
    dioxygen.add_atom(8, 0, False)  # O0
    dioxygen.add_atom(8, 0, False)  # O1
    dioxygen.add_bond(0, 1, 2, False)  # O0 # O1
    charges = unit.Quantity(np.array([0.0, 0.0]), unit.elementary_charge)
    dioxygen.partial_charges = charges

    return dioxygen


@pytest.fixture(scope="function")
def dinitrogen():
    """
    Creates an openforcefield.topology.Molecule representation of
    dinitrogen without the use of a cheminformatics toolkit
    """
    from openforcefield.topology.molecule import Molecule

    dinitrogen = Molecule()
    dinitrogen.add_atom(7, 0, False)  # N0
    dinitrogen.add_atom(7, 0, False)  # N1
    dinitrogen.add_bond(0, 1, 3, False)  # N0 - N1
    charges = unit.Quantity(np.array([0.0, 0.0]), unit.elementary_charge)
    dinitrogen.partial_charges = charges
    return dinitrogen


@pytest.fixture(scope="function")
def cis_1_2_dichloroethene():
    """
    Creates an openforcefield.topology.Molecule representation of cis-1,2-dichloroethene
    without the use of a cheminformatics toolkit.
    """
    from openforcefield.topology.molecule import Molecule

    cis_dichloroethene = Molecule()
    cis_dichloroethene.add_atom(17, 0, False)
    cis_dichloroethene.add_atom(6, 0, False)
    cis_dichloroethene.add_atom(6, 0, False)
    cis_dichloroethene.add_atom(17, 0, False)
    cis_dichloroethene.add_atom(1, 0, False)
    cis_dichloroethene.add_atom(1, 0, False)
    cis_dichloroethene.add_bond(0, 1, 1, False)
    cis_dichloroethene.add_bond(1, 2, 2, False, "Z")
    cis_dichloroethene.add_bond(2, 3, 1, False)
    cis_dichloroethene.add_bond(1, 4, 1, False)
    cis_dichloroethene.add_bond(2, 5, 1, False)
    return cis_dichloroethene
