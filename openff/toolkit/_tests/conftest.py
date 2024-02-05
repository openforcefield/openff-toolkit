"""
Configuration file for pytest.
"""

import logging

logger = logging.getLogger(__name__)

# Ensure QCPortal is imported before any OpenEye modules, see
# https://github.com/conda-forge/qcfractal-feedstock/issues/43
try:
    import qcportal  # noqa
except ImportError:
    pass


def untar_full_alkethoh_and_freesolv_set():
    """When running slow tests, we unpack the full AlkEthOH and FreeSolv test
    sets in advance to speed things up.

    See
        test_forcefield.py::test_alkethoh_parameters_assignment
        test_forcefield.py::test_freesolv_parameters_assignment
    """
    import os
    import tarfile

    from openff.toolkit.utils import get_data_file_path

    molecule_dir_path = get_data_file_path("molecules")
    for tarfile_name in ["AlkEthOH_tripos.tar.gz", "FreeSolv.tar.gz"]:
        tarfile_path = os.path.join(molecule_dir_path, tarfile_name)
        with tarfile.open(tarfile_path, "r:gz") as tar:
            tar.extractall(path=molecule_dir_path)
