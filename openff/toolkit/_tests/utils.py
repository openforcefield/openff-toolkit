"""
Utilities for testing.

"""
import os
from contextlib import contextmanager

import pytest
from openff.utilities import has_package, skip_if_missing

from openff.toolkit.utils import (
    AmberToolsToolkitWrapper,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    get_data_file_path,
)

requires_ambertools = pytest.mark.skipif(
    not AmberToolsToolkitWrapper.is_available(),
    reason="Test requires AmberTools",
)
requires_rdkit = pytest.mark.skipif(
    not RDKitToolkitWrapper.is_available(),
    reason="Test requires RDKit",
)
requires_openeye = pytest.mark.skipif(
    not OpenEyeToolkitWrapper.is_available(),
    reason="Test requires OE toolkit",
)
requires_openeye_mol2 = pytest.mark.skipif(
    not OpenEyeToolkitWrapper.is_available(),
    reason="Test requires OE toolkit to read mol2 files",
)

has_pkg = has_package
requires_pkg = skip_if_missing

if has_package("openmm"):
    from openff.toolkit._tests.openmm_utils import (
        _compare_parameters,
        _find_all_bonds,
        _get_angle_force_parameters,
        _get_bond_force_parameters,
        _get_force_parameters,
        _get_improper_torsion_canonical_order,
        _get_nonbonded_force_parameters,
        _get_proper_torsion_canonical_order,
        _get_torsion_force_parameters,
        _merge_impropers_folds,
        compare_amber_smirnoff,
        compare_context_energies,
        compare_partial_charges,
        compare_system_energies,
        compare_system_parameters,
        coords_from_off_mols,
        create_system_from_amber,
        evaluate_molecules_off,
        evaluate_molecules_omm,
        get_14_scaling_factors,
        get_context_potential_energy,
        insert_vsite_padding,
        openmm_evaluate_vsites_and_energy,
        quantities_allclose,
        reorder_openff_to_openmm,
        reorder_openmm_to_openff,
    )


@contextmanager
def does_not_raise():
    """A helpful context manager to use inplace of a pytest raise statement
    when no exception is expected."""
    yield


def get_packmol_pdb_file_path(prefix="cyclohexane_ethanol_0.4_0.6"):
    """Get PDB filename for a packmol-generated box

    Parameters
    ----------
    prefix : str, optional, default='cyclohexane_ethanol_0.4_0.6'
        The prefix of .pdb file to retrieve from testdata/systems/packmol_boxes

    Returns
    -------
    pdb_filename : str
        Absolute path to the PDB file
    """
    prefix = os.path.join("systems", "packmol_boxes", prefix)
    pdb_filename = get_data_file_path(prefix + ".pdb")
    return pdb_filename


def extract_compressed_molecules(tar_file_name, file_subpaths=None, filter_func=None):
    if (file_subpaths is None) == (filter_func is None):
        raise ValueError(
            "Only one between file_subpaths and filter_func must be specified."
        )

    # Find the path of the tarfile with respect to the data/molecules/ folder.
    molecules_dir_path = get_data_file_path("molecules")
    tar_file_path = os.path.join(molecules_dir_path, tar_file_name)
    tar_root_dir_name = tar_file_name.split(".")[0]

    # Return value: Paths to the extracted molecules.
    extracted_file_paths = None

    # Handle subpaths search.
    if file_subpaths is not None:
        # We can already check the paths of the extracted files
        # and skipping opening the tarball if not necessary.
        extracted_file_paths = [
            os.path.join(molecules_dir_path, tar_root_dir_name, file_subpath)
            for file_subpath in file_subpaths
        ]

        # Remove files that we have already extracted.
        # Also, we augument the subpath with its root directory.
        file_subpaths_set = {
            os.path.join(tar_root_dir_name, subpath)
            for subpath, fullpath in zip(file_subpaths, extracted_file_paths)
            if not os.path.isfile(fullpath)
        }

        # If everything was already extracted, we don't need to open the tarball.
        if len(file_subpaths_set) == 0:
            return extracted_file_paths

        # Otherwise, create a filter matching only the subpaths.
        def filter_func(x):
            return x in file_subpaths_set

    # If no filter was specified, just create one matching everything.
    if filter_func is None:

        def filter_func(x):
            return True

    # Determine opening mode.
    if ".gz" in tar_file_name:
        mode = "r:gz"
    else:
        mode = "r"

    # Extract required files.
    import tarfile

    with tarfile.open(tar_file_path, mode) as tar_file:
        # Gather the paths to extract. Remove the
        members = [m for m in tar_file.getmembers() if filter_func(m.name)]

        # Built the paths to the extracted molecules we didn't already.
        if extracted_file_paths is None:
            extracted_file_paths = [
                os.path.join(molecules_dir_path, m.name) for m in members
            ]

        # Extract only the members that we haven't already extracted.
        members = [
            member
            for member, fullpath in zip(members, extracted_file_paths)
            if not os.path.isfile(fullpath)
        ]
        tar_file.extractall(path=molecules_dir_path, members=members)

    return extracted_file_paths


def get_alkethoh_file_path(alkethoh_name, get_amber=False):
    """Retrieve the mol2, top and crd files of a molecule in the AlkEthOH set.

    Parameters
    ----------
    alkethoh_name : str
        The name of the AlkEthOH molecule (e.g. "AlkEthOH_r0", "AlkEthOH_c1266").
    get_amber : bool, optional
        If True, the paths to the top and crd files are returned.

    Returns
    -------
    molecule_file_paths : str or List[str]
        All the requested paths. If ``get_amber`` is False, only a single string
        pointing to the path of the mol2 file is returned, otherwise this is a
        list ``[mol2_path, top_path, crd_path]``.

    """
    # Determine if this is a ring or a chain molecule and the subfolder name.
    is_ring = alkethoh_name[9] == "r"
    alkethoh_subdir_name = "rings" if is_ring else "chain"
    alkethoh_subdir_name = "AlkEthOH_" + alkethoh_subdir_name + "_filt1"

    # Determine which paths have to be returned.
    file_base_subpath = os.path.join(alkethoh_subdir_name, alkethoh_name)
    # We always return the mol2 file.
    file_subpaths = [file_base_subpath + "_tripos.mol2"]
    # Check if we need to return also Amber files.
    if get_amber:
        file_subpaths.append(file_base_subpath + ".top")
        file_subpaths.append(file_base_subpath + ".crd")

    return extract_compressed_molecules(
        "AlkEthOH_tripos.tar.gz", file_subpaths=file_subpaths
    )


def get_freesolv_file_path(freesolv_id, ff_version):
    file_base_name = "mobley_" + freesolv_id
    mol2_file_subpath = os.path.join("mol2files_sybyl", file_base_name + ".mol2")
    xml_dir = "xml_" + ff_version.replace(".", "_")
    xml_file_subpath = os.path.join(xml_dir, file_base_name + "_vacuum.xml")

    # Extract the files if needed.
    file_subpaths = [mol2_file_subpath, xml_file_subpath]
    return extract_compressed_molecules("FreeSolv.tar.gz", file_subpaths=file_subpaths)
