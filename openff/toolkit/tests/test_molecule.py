#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Tests for molecular topology representations

At least one supported cheminformatics toolkit must be installed to run these tests.
Only the tests applicable to that toolkit will be run.

TODO:
- Will the ToolkitWrapper allow us to pare down importing each wrapper directly?
- Add tests comparing RDKit and OpenEye aromaticity perception
- Right now, the test database of TestMolecule is read from mol2, requiring the OE
  toolkit. Find a different test set that RDKit can read, or make a database of
  serialized OFFMols.

"""

# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import copy
import os
import pickle
from tempfile import NamedTemporaryFile

import numpy as np
import pytest
from simtk import unit

from openff.toolkit.tests.create_molecules import (
    create_acetaldehyde,
    create_benzene_no_aromatic,
    create_cis_1_2_dichloroethene,
    create_cyclohexane,
    create_ethanol,
    create_reversed_ethanol,
)
from openff.toolkit.tests.utils import (
    has_pkg,
    requires_openeye,
    requires_pkg,
    requires_rdkit,
)
from openff.toolkit.topology import NotBondedError
from openff.toolkit.topology.molecule import (
    Atom,
    FrozenMolecule,
    InvalidConformerError,
    Molecule,
    SmilesParsingError,
)
from openff.toolkit.utils import get_data_file_path
from openff.toolkit.utils.toolkits import (
    AmberToolsToolkitWrapper,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    ToolkitRegistry,
)

# =============================================================================================
# TEST UTILITIES
# =============================================================================================


def assert_molecule_is_equal(molecule1, molecule2, msg):
    """Compare whether two Molecule objects are equal

    Parameters
    ----------
    molecule1, molecule2 : openff.toolkit.topology.Molecule
        Molecules to be compared
    msg : str
        Message to include if molecules fail to match.

    """
    if not (molecule1.is_isomorphic_with(molecule2)):
        raise AssertionError(msg)


def is_four_membered_ring_torsion(torsion):
    """Check that three atoms in the given torsion form a four-membered ring."""
    # Push a copy of the first and second atom in the end to make the code simpler.
    torsion = list(torsion) + [torsion[0], torsion[1]]

    is_four_membered_ring = True
    for i in range(4):
        # The atom is bonded to the next one.
        is_four_membered_ring &= torsion[i].is_bonded_to(torsion[i + 1])
        # The atom is not bonded to the atom on its diagonal.
        is_four_membered_ring &= not torsion[i].is_bonded_to(torsion[i + 2])

    return is_four_membered_ring


def is_three_membered_ring_torsion(torsion):
    """Check that three atoms in the given torsion form a three-membered ring.

    In order to be 4 atoms with a three-membered ring, there must be
    1) A central atom connected to all other atoms.
    2) An atom outside the ring connected exclusively to the central atom.
    3) Two atoms in the ring connected to the central atom and to each other.

    """
    # A set of atom indices for the atoms in the torsion.
    torsion_atom_indices = set(a.molecule_atom_index for a in torsion)

    # Collect all the bonds involving exclusively atoms in the torsion.
    bonds_by_atom_idx = {i: set() for i in torsion_atom_indices}
    for atom in torsion:
        for bond in atom.bonds:
            # Consider the bond only if both atoms are in the torsion.
            if (
                bond.atom1_index in torsion_atom_indices
                and bond.atom2_index in torsion_atom_indices
            ):
                bonds_by_atom_idx[bond.atom1_index].add(bond.atom2_index)
                bonds_by_atom_idx[bond.atom2_index].add(bond.atom1_index)

    # Find the central atom, which is connected to all other atoms.
    atom_indices = [i for i in torsion_atom_indices if len(bonds_by_atom_idx[i]) == 3]
    if len(atom_indices) != 1:
        return False
    central_atom_idx = atom_indices[0]

    # Find the atom outside the ring.
    atom_indices = [i for i in torsion_atom_indices if len(bonds_by_atom_idx[i]) == 1]
    if (
        len(atom_indices) != 1
        or central_atom_idx not in bonds_by_atom_idx[atom_indices[0]]
    ):
        return False
    outside_atom_idx = atom_indices[0]

    # Check that the remaining two atoms are non-central atoms in the membered ring.
    atom1, atom2 = [
        i for i in torsion_atom_indices if i not in [central_atom_idx, outside_atom_idx]
    ]
    # The two atoms are bonded to each other.
    if atom2 not in bonds_by_atom_idx[atom1] or atom1 not in bonds_by_atom_idx[atom2]:
        return False
    # Check that they are both bonded to the central atom and none other.
    for atom_idx in [atom1, atom2]:
        if (
            central_atom_idx not in bonds_by_atom_idx[atom_idx]
            or len(bonds_by_atom_idx[atom_idx]) != 2
        ):
            return False

    # This is a torsion including a three-membered ring.
    return True


# =============================================================================================
# FIXTURES
# =============================================================================================


def mini_drug_bank(xfail_mols=None, wip_mols=None):
    """Load the full MiniDrugBank into Molecule objects.

    Parameters
    ----------
    xfail_mols : Dict[str, str or None]
        Dictionary mapping the molecule names that are allowed to
        failed to the failure reason.
    wip_mols : Dict[str, str or None]
        Dictionary mapping the molecule names that are work in progress
        to the failure reason.

    """
    # If we have already loaded the data set, return the cached one.
    if mini_drug_bank.molecules is not None:
        molecules = mini_drug_bank.molecules
    else:
        # Load the dataset.
        file_path = get_data_file_path("molecules/MiniDrugBank_tripos.mol2")
        try:
            # We need OpenEye to parse the molecules, but pytest execute this
            # whether or not the test class is skipped so if OE is not available
            # we just return an empty list of test cases as a workaround.
            molecules = Molecule.from_file(file_path, allow_undefined_stereo=True)
        except NotImplementedError as e:
            assert "No toolkits in registry can read file" in str(e)
            mini_drug_bank.molecules = []
            return []
        else:
            mini_drug_bank.molecules = molecules

    # Check if we need to mark anything.
    if xfail_mols is None and wip_mols is None:
        return molecules

    # Handle mutable default.
    if xfail_mols is None:
        xfail_mols = {}
    if wip_mols is None:
        wip_mols = {}
    # There should be no molecule in both dictionaries.
    assert len(set(xfail_mols).intersection(set(wip_mols))) == 0

    # Don't modify the cached molecules.
    molecules = copy.deepcopy(molecules)
    for i, mol in enumerate(molecules):
        if mol.name in xfail_mols:
            marker = pytest.mark.xfail(reason=xfail_mols[mol.name])
        elif mol.name in wip_mols:
            marker = pytest.mark.wip(reason=wip_mols[mol.name])
        else:
            marker = None

        if marker is not None:
            molecules[i] = pytest.param(mol, marks=marker)

    return molecules


# Use a "static" variable as a workaround as fixtures cannot be
# used inside pytest.mark.parametrize (see issue #349 in pytest).
mini_drug_bank.molecules = None

# All the molecules that raise UndefinedStereochemistryError when read by OETK()
openeye_drugbank_undefined_stereo_mols = {
    "DrugBank_1634",
    "DrugBank_1700",
    "DrugBank_1962",
    "DrugBank_2519",
    "DrugBank_2987",
    "DrugBank_3502",
    "DrugBank_3930",
    "DrugBank_4161",
    "DrugBank_4162",
    "DrugBank_5043",
    "DrugBank_5418",
    "DrugBank_6531",
}

# All the molecules that raise UndefinedStereochemistryError when read by OETK().
# Note that this list is different from that for OEMol,
# since the toolkits have different definitions of "stereogenic"
rdkit_drugbank_undefined_stereo_mols = {
    "DrugBank_1634",
    "DrugBank_1962",
    "DrugBank_2519",
    "DrugBank_3930",
    "DrugBank_5043",
    "DrugBank_5418",
    "DrugBank_7124",
    "DrugBank_6865",
}


# Missing stereo in OE but not RDK:  'DrugBank_2987', 'DrugBank_3502', 'DrugBank_4161',
# 'DrugBank_4162', 'DrugBank_6531', 'DrugBank_1700',
drugbank_stereogenic_in_rdkit_but_not_openeye = {
    "DrugBank_5329",
    "DrugBank_7124",
    "DrugBank_6865",
}

# Some molecules are _valid_ in both OETK and RDKit, but will fail if you try
# to convert from one to the other, since OE adds stereo that RDKit doesn't
drugbank_stereogenic_in_oe_but_not_rdkit = {
    "DrugBank_1598",
    "DrugBank_4346",
    "DrugBank_1849",
    "DrugBank_2141",
}

# =============================================================================================
# TESTS
# =============================================================================================


class TestAtom:
    """Test Atom class."""

    def test_atom_constructor(self):
        """Test Atom creation"""
        # Create a non-aromatic carbon atom
        atom1 = Atom(6, 0, False)
        assert atom1.atomic_number == 6
        assert atom1.formal_charge == 0 * unit.elementary_charge

        # Create a chiral carbon atom
        atom2 = Atom(6, 0, False, stereochemistry="R", name="CT")
        assert atom1.stereochemistry != atom2.stereochemistry

        # Ensure that formal charge can also be set as a Quantity
        atom1 = Atom(6, 1 * unit.elementary_charge, False)
        assert atom1.formal_charge == 1 * unit.elementary_charge

    def test_atom_properties(self):
        """Test that atom properties are correctly populated and gettable"""
        from simtk.openmm.app import element

        formal_charge = 0 * unit.elementary_charge
        is_aromatic = False
        # Attempt to create all elements supported by OpenMM
        elements = [
            getattr(element, name)
            for name in dir(element)
            if (type(getattr(element, name)) == element.Element)
        ]
        # The above runs into a problem with deuterium (fails name assertion)
        elements.remove(element.deuterium)
        for this_element in elements:
            atom = Atom(
                this_element.atomic_number,
                formal_charge,
                is_aromatic,
                name=this_element.name,
            )
            assert atom.atomic_number == this_element.atomic_number
            assert atom.element == this_element
            assert atom.mass == this_element.mass
            assert atom.formal_charge == formal_charge
            assert atom.is_aromatic == is_aromatic
            assert atom.name == this_element.name

    def test_set_molecule(self):
        """Test appropriately setting a molecule with no errors"""
        mol = Molecule.from_smiles("CCO")
        atom = Atom(6, 0, False)
        atom.molecule = mol

    def test_set_molecule_error(self):
        """Test setting molecule for atom with molecule raises error"""
        mol = Molecule.from_smiles("CCO")
        atom = Atom(6, 0, False)
        atom.molecule = mol
        with pytest.raises(AssertionError, match="already has an associated molecule"):
            atom.molecule = mol


class TestMolecule:
    """Test Molecule class."""

    # TODO: Test getstate/setstate

    # Test serialization {to|from}_{dict|yaml|toml|json|bson|xml|messagepack|pickle}

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_dict_serialization(self, molecule):
        """Test serialization of a molecule object to and from dict."""
        serialized = molecule.to_dict()
        molecule_copy = Molecule.from_dict(serialized)
        assert molecule == molecule_copy
        assert molecule_copy.n_conformers == molecule.n_conformers
        assert np.allclose(molecule_copy.conformers[0], molecule.conformers[0])

    @requires_pkg("yaml")
    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_yaml_serialization(self, molecule):
        """Test serialization of a molecule object to and from YAML."""
        serialized = molecule.to_yaml()
        molecule_copy = Molecule.from_yaml(serialized)
        assert molecule == molecule_copy
        assert molecule_copy.n_conformers == molecule.n_conformers
        assert np.allclose(molecule_copy.conformers[0], molecule.conformers[0])

    @requires_pkg("toml")
    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_toml_serialization(self, molecule):
        """Test serialization of a molecule object to and from TOML."""
        # TODO: Test round-trip, on mini_drug_bank, when implemented
        mol = Molecule.from_smiles("CCO")
        with pytest.raises(NotImplementedError):
            mol.to_toml()

    @requires_pkg("bson")
    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_bson_serialization(self, molecule):
        """Test serialization of a molecule object to and from BSON."""
        serialized = molecule.to_bson()
        molecule_copy = Molecule.from_bson(serialized)
        assert molecule == molecule_copy
        assert molecule_copy.n_conformers == molecule.n_conformers
        assert np.allclose(molecule_copy.conformers[0], molecule.conformers[0])

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_json_serialization(self, molecule):
        """Test serialization of a molecule object to and from JSON."""
        molecule_copy = Molecule.from_json(molecule.to_json())
        assert molecule_copy == molecule
        assert molecule_copy.n_conformers == molecule.n_conformers
        assert np.allclose(molecule_copy.conformers[0], molecule.conformers[0])

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_xml_serialization(self, molecule):
        """Test serialization of a molecule object to and from XML."""
        # TODO: Test round-trip, on mini_drug_bank, when from_xml is implemented
        mol = Molecule.from_smiles("CCO")
        serialized = mol.to_xml()
        with pytest.raises(NotImplementedError):
            Molecule.from_xml(serialized)

    @requires_pkg("msgpack")
    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_messagepack_serialization(self, molecule):
        """Test serialization of a molecule object to and from messagepack."""
        serialized = molecule.to_messagepack()
        molecule_copy = Molecule.from_messagepack(serialized)
        assert molecule == molecule_copy
        assert molecule_copy.n_conformers == molecule.n_conformers
        assert np.allclose(molecule_copy.conformers[0], molecule.conformers[0])

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_pickle_serialization(self, molecule):
        """Test round-trip pickling of a molecule object."""
        serialized = pickle.dumps(molecule)
        molecule_copy = pickle.loads(serialized)
        assert molecule == molecule_copy
        assert molecule_copy.n_conformers == molecule.n_conformers
        assert np.allclose(molecule_copy.conformers[0], molecule.conformers[0])

    @requires_pkg("yaml")
    @requires_pkg("toml")
    @requires_pkg("msgpack")
    def test_serialization_no_conformers(self):
        """Test round-trip serialization when molecules have no conformers or partial charges."""
        mol = Molecule.from_smiles("CCO")

        dict_copy = Molecule.from_dict(mol.to_dict())
        assert mol == dict_copy

        # TODO: yaml_copy = Molecule.from_yaml(mol.to_yaml())
        with pytest.raises(NotImplementedError):
            mol.to_toml()

        bson_copy = Molecule.from_bson(mol.to_bson())
        assert mol == bson_copy

        json_copy = Molecule.from_json(mol.to_json())
        assert mol == json_copy

        # TODO: round-trip when from_xml is implemented
        mol_as_xml = mol.to_xml()
        with pytest.raises(NotImplementedError):
            Molecule.from_xml(mol_as_xml)

        messagepack_copy = Molecule.from_messagepack(mol.to_messagepack())
        assert mol == messagepack_copy

        pickle_copy = pickle.loads(pickle.dumps(mol))
        assert mol == pickle_copy

    def test_json_numpy_roundtrips(self):
        """Ensure that array data survives several round-trips through JSON,
        which depends on list serialization instead of bytes."""
        mol = Molecule.from_smiles("CCO")
        mol.generate_conformers(n_conformers=1)
        initial_conformer = mol.conformers[0]

        for _ in range(10):
            mol = Molecule.from_json(mol.to_json())

        assert np.allclose(initial_conformer, mol.conformers[0])

    # ----------------------------------------------------
    # Test Molecule constructors and conversion utilities.
    # ----------------------------------------------------

    def test_create_empty(self):
        """Test empty constructor."""
        molecule = Molecule()
        assert len(molecule.atoms) == 0
        assert len(molecule.bonds) == 0

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_create_copy(self, molecule):
        """Test copy constructor."""
        molecule_copy = Molecule(molecule)
        assert molecule_copy == molecule

        # Test that the "properties" dict of both molecules is unique
        # (see https://github.com/openforcefield/openff-toolkit/pull/786)
        molecule_copy.properties["aaa"] = "bbb"
        assert "aaa" not in molecule.properties

    @pytest.mark.parametrize("toolkit", [OpenEyeToolkitWrapper, RDKitToolkitWrapper])
    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_to_from_smiles(self, molecule, toolkit):
        """Test round-trip creation from SMILES"""
        if not toolkit.is_available():
            pytest.skip("Required toolkit is unavailable")

        if toolkit == RDKitToolkitWrapper:
            # Skip the test if OpenEye assigns stereochemistry but RDKit doesn't (since then, the
            # OFF molecule will be loaded, but fail to convert in to_rdkit)
            if molecule.name in drugbank_stereogenic_in_oe_but_not_rdkit:
                pytest.skip(
                    "Molecule is stereogenic in OpenEye (which loaded this dataset), but not RDKit, so it "
                    "is impossible to make a valid RDMol in this test"
                )
            undefined_stereo_mols = rdkit_drugbank_undefined_stereo_mols
        elif toolkit == OpenEyeToolkitWrapper:
            undefined_stereo_mols = openeye_drugbank_undefined_stereo_mols

        toolkit_wrapper = toolkit()

        undefined_stereo = molecule.name in undefined_stereo_mols
        # Since OpenEye did the original reading of MiniDrugBank, if OPENEYE doesn't
        # think a feature is stereogenic, then "molecule" won't have stereochemistry defined
        stereogenic_in_rdk_but_not_openeye = (
            molecule.name in drugbank_stereogenic_in_rdkit_but_not_openeye
        )

        smiles1 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        if undefined_stereo or (
            (toolkit is RDKitToolkitWrapper) and stereogenic_in_rdk_but_not_openeye
        ):
            molecule2 = Molecule.from_smiles(
                smiles1, allow_undefined_stereo=True, toolkit_registry=toolkit_wrapper
            )
        else:
            molecule2 = Molecule.from_smiles(smiles1, toolkit_registry=toolkit_wrapper)

        smiles2 = molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles1 == smiles2

    @pytest.mark.parametrize(
        "smiles, expected", [("[Cl:1]Cl", {0: 1}), ("[Cl:1][Cl:2]", {0: 1, 1: 2})]
    )
    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_from_smiles_with_map(self, smiles, expected, toolkit_class):
        if not (toolkit_class.is_available()):
            pytest.skip(f"Required toolkit {toolkit_class} is unavailable")
        molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_class())
        assert molecule.properties["atom_map"] == expected

    smiles_types = [
        {"isomeric": True, "explicit_hydrogens": True, "mapped": True, "error": None},
        {"isomeric": False, "explicit_hydrogens": True, "mapped": True, "error": None},
        {
            "isomeric": True,
            "explicit_hydrogens": False,
            "mapped": True,
            "error": AssertionError,
        },
        {"isomeric": True, "explicit_hydrogens": True, "mapped": False, "error": None},
        {"isomeric": True, "explicit_hydrogens": False, "mapped": False, "error": None},
        {"isomeric": False, "explicit_hydrogens": True, "mapped": False, "error": None},
        {
            "isomeric": False,
            "explicit_hydrogens": False,
            "mapped": True,
            "error": AssertionError,
        },
        {
            "isomeric": False,
            "explicit_hydrogens": False,
            "mapped": False,
            "error": None,
        },
    ]

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    @pytest.mark.parametrize("data", smiles_types)
    def test_smiles_types(self, data, toolkit_class):
        """Test that the toolkit is passing the correct args to the toolkit backends across different combinations."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = create_cis_1_2_dichloroethene()
            isomeric, explicit_hs, mapped = (
                data["isomeric"],
                data["explicit_hydrogens"],
                data["mapped"],
            )
            if data["error"] is not None:
                with pytest.raises(data["error"]):
                    mol.to_smiles(
                        isomeric=isomeric,
                        explicit_hydrogens=explicit_hs,
                        mapped=mapped,
                        toolkit_registry=toolkit,
                    )

            else:

                # make the smiles then do some checks on it
                output_smiles = mol.to_smiles(
                    isomeric=isomeric,
                    explicit_hydrogens=explicit_hs,
                    mapped=mapped,
                    toolkit_registry=toolkit,
                )
                if isomeric:
                    assert "\\" in output_smiles
                if explicit_hs:
                    assert "H" in output_smiles
                if mapped:
                    for i in range(1, 7):
                        assert f":{i}" in output_smiles
                    # if the molecule is mapped make it using the mapping
                    mol2 = Molecule.from_mapped_smiles(
                        mapped_smiles=output_smiles,
                        toolkit_registry=toolkit,
                        allow_undefined_stereo=not isomeric,
                    )
                else:
                    # make a molecule from a standard smiles
                    mol2 = Molecule.from_smiles(
                        smiles=output_smiles,
                        allow_undefined_stereo=not isomeric,
                        toolkit_registry=toolkit,
                    )

                isomorphic, atom_map = Molecule.are_isomorphic(
                    mol,
                    mol2,
                    return_atom_map=True,
                    aromatic_matching=True,
                    formal_charge_matching=True,
                    bond_order_matching=True,
                    atom_stereochemistry_matching=isomeric,
                    bond_stereochemistry_matching=isomeric,
                )

                assert isomorphic is True
                if mapped:
                    assert {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5} == atom_map

        else:
            pytest.skip(
                f"The required toolkit ({toolkit_class.toolkit_name}) is not available."
            )

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_smiles_cache(self, toolkit_class):
        """Make sure that the smiles cache is being used correctly."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            # this uses no toolkit back end so no smiles should be saved
            mol = create_ethanol()

            # now lets populate the cache with a test result
            # first we need to make the cache key for the default input
            isomeric, explicit_hydrogens, mapped = True, True, False
            cache_key = (
                toolkit.to_smiles.__qualname__
                + str(isomeric)
                + str(explicit_hydrogens)
                + str(mapped)
            )
            cache_key += str(mol._properties.get("atom_map", None))
            mol._cached_smiles = {cache_key: None}
            assert (
                mol.to_smiles(
                    isomeric=isomeric,
                    toolkit_registry=toolkit,
                    explicit_hydrogens=explicit_hydrogens,
                    mapped=mapped,
                )
                is None
            )

            # now make sure the cache is not used if we change an input arg
            assert (
                mol.to_smiles(
                    isomeric=True,
                    explicit_hydrogens=True,
                    mapped=True,
                    toolkit_registry=toolkit,
                )
                is not None
            )

            # now make sure the cache was updated
            assert mol._cached_smiles != {cache_key: None}
            assert len(mol._cached_smiles) == 2

        else:
            pytest.skip(
                f"The required toolkit ({toolkit_class.toolkit_name}) is not available."
            )

    mapped_types = [
        {"atom_map": None},
        {"atom_map": {0: 0}},
        {"atom_map": {0: 0, 1: 0, 2: 0, 3: 0}},
        {"atom_map": {0: 0, 1: 1, 2: 2, 3: 3}},
        {"atom_map": {0: 1, 1: 2, 2: 3, 3: 4}},
    ]

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    @pytest.mark.parametrize("data", mapped_types)
    def test_partial_mapped_smiles(self, toolkit_class, data):

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = create_cis_1_2_dichloroethene()
            mol._properties["atom_map"] = data["atom_map"]

            smiles = mol.to_smiles(
                isomeric=True,
                explicit_hydrogens=True,
                mapped=True,
                toolkit_registry=toolkit,
            )

            # now we just need to check the smiles generated
            if data["atom_map"] is None:
                for i, atom in enumerate(mol.atoms, 1):
                    assert f"[{atom.element.symbol}:{i}]" in smiles
            else:
                if 0 in data["atom_map"].values():
                    increment = True
                else:
                    increment = False

                for atom, index in data["atom_map"].items():
                    assert f"[{mol.atoms[atom].element.symbol}:{index + 1 if increment else index}]"

        else:
            pytest.skip(
                f"The required toolkit ({toolkit_class.toolkit_name}) is not available."
            )

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_unique_atom_names(self, molecule):
        """Test molecules have unique atom names"""
        # The dataset we load in has atom names, so let's strip them first
        # to ensure that we can fail the uniqueness check
        for atom in molecule.atoms:
            atom.name = ""
        assert not (molecule.has_unique_atom_names)
        # Then genreate unique atom names using the built in algorithm
        molecule.generate_unique_atom_names()
        # Check that the molecule has unique atom names
        assert molecule.has_unique_atom_names
        # Check molecule.has_unique_atom_names is working correctly
        assert (
            len(set([atom.name for atom in molecule.atoms])) == molecule.n_atoms
        ) == molecule.has_unique_atom_names
        molecule.atoms[1].name = molecule.atoms[0].name  # no longer unique
        assert (
            len(set([atom.name for atom in molecule.atoms])) == molecule.n_atoms
        ) == molecule.has_unique_atom_names

    inchi_data = [
        {
            "molecule": create_ethanol(),
            "standard_inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            "fixed_hydrogen_inchi": "InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3",
        },
        {
            "molecule": create_reversed_ethanol(),
            "standard_inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            "fixed_hydrogen_inchi": "InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3",
        },
        {
            "molecule": create_acetaldehyde(),
            "standard_inchi": "InChI=1S/C2H4O/c1-2-3/h2H,1H3",
            "fixed_hydrogen_inchi": "InChI=1/C2H4O/c1-2-3/h2H,1H3",
        },
        {
            "molecule": create_cyclohexane(),
            "standard_inchi": "InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2",
            "fixed_hydrogen_inchi": "InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2",
        },
    ]

    @requires_openeye
    @pytest.mark.parametrize("data", inchi_data)
    def test_from_inchi(self, data):
        """Test building a molecule from standard and non-standard InChI strings."""

        toolkit = OpenEyeToolkitWrapper()
        ref_mol = data["molecule"]
        # make a molecule from inchi
        inchi_mol = Molecule.from_inchi(
            data["standard_inchi"], toolkit_registry=toolkit
        )
        assert inchi_mol.to_inchi(toolkit_registry=toolkit) == data["standard_inchi"]

        def compare_mols(ref_mol, inchi_mol):
            assert ref_mol.n_atoms == inchi_mol.n_atoms
            assert ref_mol.n_bonds == inchi_mol.n_bonds
            assert ref_mol.n_angles == inchi_mol.n_angles
            assert ref_mol.n_propers == inchi_mol.n_propers
            assert ref_mol.is_isomorphic_with(inchi_mol) is True

        compare_mols(ref_mol, inchi_mol)

        # now make the molecule from the non-standard inchi and compare
        nonstandard_inchi_mol = Molecule.from_inchi(
            data["fixed_hydrogen_inchi"], toolkit_registry=toolkit
        )
        assert (
            nonstandard_inchi_mol.to_inchi(
                fixed_hydrogens=True, toolkit_registry=toolkit
            )
            == data["fixed_hydrogen_inchi"]
        )

        compare_mols(ref_mol, nonstandard_inchi_mol)

    # TODO: Should there be an equivalent toolkit test and leave this as an integration test?
    @requires_openeye
    @pytest.mark.slow
    def test_create_from_file(self):
        """Test standard constructor taking a filename or file-like object."""
        # TODO: Expand test to both openeye and rdkit toolkits
        filename = get_data_file_path("molecules/toluene.mol2")

        molecule1 = Molecule(filename, allow_undefined_stereo=True)
        with open(filename, "r") as infile:
            molecule2 = Molecule(
                infile, file_format="MOL2", allow_undefined_stereo=True
            )
        assert molecule1 == molecule2

        import gzip

        with gzip.GzipFile(filename + ".gz", "r") as infile:
            molecule3 = Molecule(
                infile, file_format="MOL2", allow_undefined_stereo=True
            )
        assert molecule3 == molecule1

        # Ensure that attempting to initialize a single Molecule from a file
        # containing multiple molecules raises a ValueError
        with pytest.raises(ValueError) as exc_info:
            filename = get_data_file_path("molecules/zinc-subset-tripos.mol2.gz")
            molecule = Molecule(filename, allow_undefined_stereo=True)

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_create_from_serialized(self, molecule):
        """Test standard constructor taking the output of __getstate__()."""
        serialized_molecule = molecule.__getstate__()
        molecule_copy = Molecule(serialized_molecule)
        assert molecule == molecule_copy

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_to_from_dict(self, molecule):
        """Test that conversion/creation of a molecule to and from a dict is consistent."""
        serialized = molecule.to_dict()
        molecule_copy = Molecule.from_dict(serialized)
        assert molecule == molecule_copy

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_to_networkx(self, molecule):
        """Test conversion to NetworkX graph."""
        graph = molecule.to_networkx()

    @requires_rdkit
    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_to_from_rdkit(self, molecule):
        """Test that conversion/creation of a molecule to and from an RDKit rdmol is consistent."""
        # import pickle
        from openff.toolkit.utils.toolkits import UndefinedStereochemistryError

        undefined_stereo = molecule.name in rdkit_drugbank_undefined_stereo_mols

        toolkit_wrapper = RDKitToolkitWrapper()

        rdmol = molecule.to_rdkit()
        molecule_smiles = molecule.to_smiles(toolkit_registry=toolkit_wrapper)

        # First test making a molecule using the Molecule(oemol) method

        # If this is a known failure, check that it raises UndefinedStereochemistryError
        # and proceed with the test ignoring it.
        test_mol = None
        if undefined_stereo:
            with pytest.raises(UndefinedStereochemistryError):
                Molecule(rdmol)
            test_mol = Molecule(rdmol, allow_undefined_stereo=True)
        else:
            test_mol = Molecule(rdmol)

        test_mol_smiles = test_mol.to_smiles(toolkit_registry=toolkit_wrapper)
        assert molecule_smiles == test_mol_smiles

        # Check that the two topologies are isomorphic.
        assert_molecule_is_equal(
            molecule, test_mol, "Molecule.to_rdkit()/Molecule(rdmol) round trip failed"
        )

        # Second, test making a molecule using the Molecule.from_openeye(oemol) method

        # If this is a known failure, check that it raises UndefinedStereochemistryError
        # and proceed with the test.
        if undefined_stereo:
            with pytest.raises(UndefinedStereochemistryError):
                Molecule.from_rdkit(rdmol)
            test_mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
        else:
            test_mol = Molecule.from_rdkit(rdmol)

        test_mol_smiles = test_mol.to_smiles(toolkit_registry=toolkit_wrapper)
        assert molecule_smiles == test_mol_smiles

        # Check that the two topologies are isomorphic.
        assert_molecule_is_equal(
            molecule, test_mol, "Molecule.to_rdkit()/from_rdkit() round trip failed"
        )

    @requires_openeye
    def test_to_from_iupac(self):
        """
        Test basic behavior of the IUPAC conversion functions. More rigorous
        testing of the toolkity wrapper behavior is in test_toolkits.py
        """
        from openff.toolkit.utils.toolkits import (
            InvalidIUPACNameError,
            UndefinedStereochemistryError,
        )

        with pytest.raises(InvalidIUPACNameError):
            Molecule.from_iupac(".BETA.-PINENE")

        # DrugBank_977, tagged as a problem molecule in earlier tests
        bad_stereo_iupac = (
            "(~{E},3~{R},5~{S})-7-[4-(4-fluorophenyl)-6-isopropyl-2-"
            "[methyl(methylsulfonyl)amino]pyrimidin-5-yl]-3,5-"
            "dihydroxy-hept-6-enoic acid"
        )
        with pytest.raises(UndefinedStereochemistryError):
            Molecule.from_iupac(bad_stereo_iupac)

        cholesterol = Molecule.from_smiles(
            "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
        )

        cholesterol_iupac = cholesterol.to_iupac()

        assert Molecule.are_isomorphic(
            cholesterol, Molecule.from_iupac(cholesterol_iupac)
        )

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_to_from_topology(self, molecule):
        """Test that conversion/creation of a molecule to and from a Topology is consistent."""
        topology = molecule.to_topology()
        molecule_copy = Molecule.from_topology(topology)
        assert molecule == molecule_copy

    @requires_openeye
    def test_to_multiframe_xyz_openeye(self):
        """
        Test writing out a molecule with multiple conformations to an xyz file

        This test is backend-specific because of precision/rounding differences between RDKit and OpenEye
        """
        from openff.toolkit.utils import OpenEyeToolkitWrapper

        tkw = OpenEyeToolkitWrapper()
        # load in an SDF of butane with multiple conformers in it
        molecules = Molecule.from_file(
            get_data_file_path("molecules/butane_multi.sdf"),
            "sdf",
            toolkit_registry=tkw,
        )
        # now we want to combine the conformers to one molecule
        butane = molecules[0]
        for mol in molecules[1:]:
            butane.add_conformer(mol._conformers[0])

        # make sure we have the 7 conformers
        assert butane.n_conformers == 7
        with NamedTemporaryFile(suffix=".xyz") as iofile:
            # try and write out the xyz file
            butane.to_file(iofile.name, "xyz", toolkit_registry=tkw)

            # now lets check whats in the file
            with open(iofile.name) as xyz_data:
                data = xyz_data.readlines()
                # make sure we have the correct amount of lines writen
                assert len(data) == 112
                # make sure all headers and frame data was writen
                assert data.count("14\n") == 7
                for i in range(1, 8):
                    assert f"C4H10 Frame {i}\n" in data

                # now make sure the first line of the coordinates are correct in every frame
                coords = [
                    "C        1.8902000189    0.0425999984    0.2431000024\n",
                    "C        1.8976000547   -0.0232999995    0.2845999897\n",
                    "C       -1.8794000149   -0.1792999953   -0.2565000057\n",
                    "C       -1.5205999613   -0.0164999999    0.2786999941\n",
                    "C       -1.4889999628   -0.2619000077    0.4871000051\n",
                    "C       -1.4940999746   -0.2249000072   -0.0957999974\n",
                    "C       -1.8826999664   -0.0372000001    0.1937000006\n",
                ]
                for coord in coords:
                    assert coord in data

    @requires_openeye
    def test_to_single_xyz_openeye(self):
        """
        Test writing to a single frame xyz file

        This test is backend-specific because of precision/rounding differences between RDKit and OpenEye
        """
        from openff.toolkit.utils import OpenEyeToolkitWrapper

        tkw = OpenEyeToolkitWrapper()

        # load a molecule with a single conformation
        toluene = Molecule.from_file(
            get_data_file_path("molecules/toluene.sdf"), "sdf", toolkit_registry=tkw
        )
        # make sure it has one conformer
        assert toluene.n_conformers == 1

        with NamedTemporaryFile(suffix=".xyz") as iofile:
            # try and write out the xyz file
            toluene.to_file(iofile.name, "xyz", toolkit_registry=tkw)

            # now lets check the file contents
            with open(iofile.name) as xyz_data:
                data = xyz_data.readlines()
                # make sure we have the correct amount of lines writen
                assert len(data) == 17
                # make sure all headers and frame data was writen
                assert data.count("15\n") == 1
                assert data.count("C7H8\n") == 1
                # now check that we can find the first and last coords
                coords = [
                    "C        0.0000000000    0.0000000000    0.0000000000\n",
                    "H       -0.0000000000    3.7604000568    0.0000000000\n",
                ]
                for coord in coords:
                    assert coord in data

    @requires_rdkit
    def test_to_multiframe_xyz_rdkit(self):
        """
        Test writing out a molecule with multiple conformations to an xyz file

        This test is backend-specific because of precision/rounding differences between RDKit and OpenEye
        """
        from openff.toolkit.utils import RDKitToolkitWrapper

        tkw = RDKitToolkitWrapper()
        # load in an SDF of butane with multiple conformers in it
        molecules = Molecule.from_file(
            get_data_file_path("molecules/butane_multi.sdf"),
            "sdf",
            toolkit_registry=tkw,
        )
        # now we want to combine the conformers to one molecule
        butane = molecules[0]
        for mol in molecules[1:]:
            butane.add_conformer(mol._conformers[0])

        # make sure we have the 7 conformers
        assert butane.n_conformers == 7
        with NamedTemporaryFile(suffix=".xyz") as iofile:
            # try and write out the xyz file
            butane.to_file(iofile.name, "xyz", toolkit_registry=tkw)

            # now lets check whats in the file
            with open(iofile.name) as xyz_data:
                data = xyz_data.readlines()
                # make sure we have the correct amount of lines writen
                assert len(data) == 112
                # make sure all headers and frame data was writen
                assert data.count("14\n") == 7
                for i in range(1, 8):
                    assert f"C4H10 Frame {i}\n" in data

                # now make sure the first line of the coordinates are correct in every frame
                coords = [
                    "C        1.8902000000    0.0426000000    0.2431000000\n",
                    "C        1.8976000000   -0.0233000000    0.2846000000\n",
                    "C       -1.8794000000   -0.1793000000   -0.2565000000\n",
                    "C       -1.5206000000   -0.0165000000    0.2787000000\n",
                    "C       -1.4890000000   -0.2619000000    0.4871000000\n",
                    "C       -1.4941000000   -0.2249000000   -0.0958000000\n",
                    "C       -1.8827000000   -0.0372000000    0.1937000000\n",
                ]
                for coord in coords:
                    assert coord in data

    @requires_rdkit
    def test_to_single_xyz_rdkit(self):
        """
        Test writing to a single frame xyz file

        This test is backend-specific because of precision/rounding differences between RDKit and OpenEye
        """
        from openff.toolkit.utils import RDKitToolkitWrapper

        tkw = RDKitToolkitWrapper()

        # load a molecule with a single conformation
        toluene = Molecule.from_file(
            get_data_file_path("molecules/toluene.sdf"), "sdf", toolkit_registry=tkw
        )
        # make sure it has one conformer
        assert toluene.n_conformers == 1

        with NamedTemporaryFile(suffix=".xyz") as iofile:
            # try and write out the xyz file
            toluene.to_file(iofile.name, "xyz", toolkit_registry=tkw)

            # now lets check the file contents
            with open(iofile.name) as xyz_data:
                data = xyz_data.readlines()
                # make sure we have the correct amount of lines writen
                assert len(data) == 17
                # make sure all headers and frame data was writen
                assert data.count("15\n") == 1
                assert data.count("C7H8\n") == 1
                # now check that we can find the first and last coords
                coords = [
                    "C        0.0000000000    0.0000000000    0.0000000000\n",
                    "H       -0.0000000000    3.7604000000    0.0000000000\n",
                ]
                for coord in coords:
                    assert coord in data

    def test_to_xyz_no_conformers(self):
        """Test writing a molecule out when it has no conformers here all coords should be 0."""

        # here we want to make a molecule with no coordinates
        ethanol = create_ethanol()
        assert ethanol.n_conformers == 0

        with NamedTemporaryFile(suffix=".xyz") as iofile:
            # try and write out the xyz file
            ethanol.to_file(iofile.name, "xyz")

            # now lets check the file contents
            with open(iofile.name) as xyz_data:
                data = xyz_data.readlines()
                # make sure we have the correct amount of lines writen
                assert len(data) == 11
                # make sure all headers and frame data was writen
                assert data.count("9\n") == 1
                assert data.count("C2H6O\n") == 1
                # now check that all coords are 0
                coords = ["0.0000000000", "0.0000000000", "0.0000000000"]
                for atom_coords in data[2:]:
                    assert atom_coords.split()[1:] == coords

    # TODO: Should there be an equivalent toolkit test and leave this as an integration test?
    @pytest.mark.parametrize("molecule", mini_drug_bank())
    @pytest.mark.parametrize(
        "format",
        [
            "mol2",
            "sdf",
            pytest.param(
                "pdb",
                marks=pytest.mark.wip(
                    reason="Read from pdb has not been implemented properly yet"
                ),
            ),
        ],
    )
    def test_to_from_file(self, molecule, format):
        """Test that conversion/creation of a molecule to and from a file is consistent."""
        from openff.toolkit.utils.toolkits import UndefinedStereochemistryError

        # TODO: Test all file capabilities; the current test is minimal
        # TODO: This is only for OE. Expand to both OE and RDKit toolkits.
        # Molecules that are known to raise UndefinedStereochemistryError.
        undefined_stereo_mols = {
            "DrugBank_1700",
            "DrugBank_2987",
            "DrugBank_3502",
            "DrugBank_4161",
            "DrugBank_4162",
            "DrugBank_6531",
        }
        undefined_stereo = molecule.name in undefined_stereo_mols

        # The file is automatically deleted outside the with-clause.
        with NamedTemporaryFile(suffix="." + format) as iofile:
            # If this has undefined stereo, check that the exception is raised.
            extension = os.path.splitext(iofile.name)[1][1:]
            molecule.to_file(iofile.name, extension)
            if undefined_stereo:
                with pytest.raises(UndefinedStereochemistryError):
                    Molecule.from_file(iofile.name)
            molecule2 = Molecule.from_file(
                iofile.name, allow_undefined_stereo=undefined_stereo
            )
            assert molecule == molecule2
            # TODO: Test to make sure properties are preserved?
            # NOTE: We can't read pdb files and expect chemical information to be preserved

    @requires_openeye
    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_to_from_oemol(self, molecule):
        """Test that conversion/creation of a molecule to and from a OEMol is consistent."""
        from openff.toolkit.utils.toolkits import UndefinedStereochemistryError

        # Known failures raise an UndefinedStereochemistryError, but
        # the round-trip SMILES representation with the OpenEyeToolkit
        # doesn't seem to be affected.
        # ZINC test set known failures.
        # known_failures = {'ZINC05964684', 'ZINC05885163', 'ZINC05543156', 'ZINC17211981',
        #                   'ZINC17312986', 'ZINC06424847', 'ZINC04963126'}

        undefined_stereo = molecule.name in openeye_drugbank_undefined_stereo_mols

        toolkit_wrapper = OpenEyeToolkitWrapper()

        oemol = molecule.to_openeye()
        molecule_smiles = molecule.to_smiles(toolkit_registry=toolkit_wrapper)

        # First test making a molecule using the Molecule(oemol) method

        # If this is a known failure, check that it raises UndefinedStereochemistryError
        # and proceed with the test ignoring it.
        test_mol = None
        if undefined_stereo:
            with pytest.raises(UndefinedStereochemistryError):
                Molecule(oemol)
            test_mol = Molecule(oemol, allow_undefined_stereo=True)
        else:
            test_mol = Molecule(oemol)

        test_mol_smiles = test_mol.to_smiles(toolkit_registry=toolkit_wrapper)
        assert molecule_smiles == test_mol_smiles

        # Check that the two topologies are isomorphic.
        assert_molecule_is_equal(
            molecule,
            test_mol,
            "Molecule.to_openeye()/Molecule(oemol) round trip failed",
        )

        # Second, test making a molecule using the Molecule.from_openeye(oemol) method

        # If this is a known failure, check that it raises UndefinedStereochemistryError
        # and proceed with the test.
        if undefined_stereo:
            with pytest.raises(UndefinedStereochemistryError):
                Molecule.from_openeye(oemol)
            test_mol = Molecule.from_openeye(oemol, allow_undefined_stereo=True)
        else:
            test_mol = Molecule.from_openeye(oemol)

        test_mol_smiles = test_mol.to_smiles(toolkit_registry=toolkit_wrapper)
        assert molecule_smiles == test_mol_smiles

        # Check that the two topologies are isomorphic.
        assert_molecule_is_equal(
            molecule, test_mol, "Molecule.to_openeye()/from_openeye() round trip failed"
        )

    # ----------------------------------------------------
    # Test properties.
    # ----------------------------------------------------

    def test_name(self):
        """Test Molecule name property"""
        molecule1 = Molecule()
        molecule1.name = None

        molecule2 = Molecule()
        molecule2.name = ""
        assert molecule1.name == molecule2.name

        name = "benzene"
        molecule = Molecule()
        molecule.name = name
        assert molecule.name == name

    def test_hill_formula(self):
        """Test that making the hill formula is consistent between input methods and ordering"""
        # make sure smiles match reference
        molecule_smiles = create_ethanol()
        assert molecule_smiles.hill_formula == "C2H6O"
        # make sure is not order dependent
        molecule_smiles_reverse = create_reversed_ethanol()
        assert molecule_smiles.hill_formula == molecule_smiles_reverse.hill_formula
        # make sure single element names are put first
        order_mol = Molecule.from_smiles("C(Br)CB")
        assert order_mol.hill_formula == "C2H6BBr"
        # test molecule with no carbon
        no_carb_mol = Molecule.from_smiles("OS(=O)(=O)O")
        assert no_carb_mol.hill_formula == "H2O4S"
        # test no carbon and hydrogen
        br_i = Molecule.from_smiles("BrI")
        assert br_i.hill_formula == "BrI"
        # make sure files and smiles match
        molecule_file = Molecule.from_file(get_data_file_path("molecules/ethanol.sdf"))
        assert molecule_smiles.hill_formula == molecule_file.hill_formula
        # make sure the topology molecule gives the same formula
        from openff.toolkit.topology.topology import Topology, TopologyMolecule

        topology = Topology.from_molecules(molecule_smiles)
        topmol = TopologyMolecule(molecule_smiles, topology)
        assert molecule_smiles.hill_formula == Molecule.to_hill_formula(topmol)
        # make sure the networkx matches
        assert molecule_smiles.hill_formula == Molecule.to_hill_formula(
            molecule_smiles.to_networkx()
        )

    def test_isomorphic_general(self):
        """Test the matching using different input types"""
        # check that hill formula fails are caught
        ethanol = create_ethanol()
        acetaldehyde = create_acetaldehyde()
        assert ethanol.is_isomorphic_with(acetaldehyde) is False
        assert acetaldehyde.is_isomorphic_with(ethanol) is False
        # check that different orderings work with full matching
        ethanol_reverse = create_reversed_ethanol()
        assert ethanol.is_isomorphic_with(ethanol_reverse) is True
        # check a reference mapping between ethanol and ethanol_reverse matches that calculated
        ref_mapping = {0: 8, 1: 7, 2: 6, 3: 3, 4: 4, 5: 5, 6: 1, 7: 2, 8: 0}
        assert (
            Molecule.are_isomorphic(ethanol, ethanol_reverse, return_atom_map=True)[1]
            == ref_mapping
        )
        # check matching with nx.Graph atomic numbers and connectivity only
        assert (
            Molecule.are_isomorphic(
                ethanol,
                ethanol_reverse.to_networkx(),
                aromatic_matching=False,
                formal_charge_matching=False,
                bond_order_matching=False,
                atom_stereochemistry_matching=False,
                bond_stereochemistry_matching=False,
            )[0]
            is True
        )
        # check matching with nx.Graph with full matching
        assert ethanol.is_isomorphic_with(ethanol_reverse.to_networkx()) is True
        # check matching with a TopologyMolecule class
        from openff.toolkit.topology.topology import Topology, TopologyMolecule

        topology = Topology.from_molecules(ethanol)
        topmol = TopologyMolecule(ethanol, topology)
        assert (
            Molecule.are_isomorphic(
                ethanol,
                topmol,
                aromatic_matching=False,
                formal_charge_matching=False,
                bond_order_matching=False,
                atom_stereochemistry_matching=False,
                bond_stereochemistry_matching=False,
            )[0]
            is True
        )
        # test hill formula passes but isomorphic fails
        mol1 = Molecule.from_smiles("Fc1ccc(F)cc1")
        mol2 = Molecule.from_smiles("Fc1ccccc1F")
        assert mol1.is_isomorphic_with(mol2) is False
        assert mol2.is_isomorphic_with(mol1) is False

    isomorphic_permutations = [
        {
            "aromatic_matching": True,
            "formal_charge_matching": True,
            "bond_order_matching": True,
            "atom_stereochemistry_matching": True,
            "bond_stereochemistry_matching": True,
            "result": False,
        },
        {
            "aromatic_matching": False,
            "formal_charge_matching": True,
            "bond_order_matching": True,
            "atom_stereochemistry_matching": True,
            "bond_stereochemistry_matching": True,
            "result": False,
        },
        {
            "aromatic_matching": True,
            "formal_charge_matching": False,
            "bond_order_matching": True,
            "atom_stereochemistry_matching": True,
            "bond_stereochemistry_matching": True,
            "result": False,
        },
        {
            "aromatic_matching": True,
            "formal_charge_matching": True,
            "bond_order_matching": False,
            "atom_stereochemistry_matching": True,
            "bond_stereochemistry_matching": True,
            "result": False,
        },
        {
            "aromatic_matching": True,
            "formal_charge_matching": True,
            "bond_order_matching": True,
            "atom_stereochemistry_matching": False,
            "bond_stereochemistry_matching": True,
            "result": False,
        },
        {
            "aromatic_matching": True,
            "formal_charge_matching": True,
            "bond_order_matching": True,
            "atom_stereochemistry_matching": True,
            "bond_stereochemistry_matching": False,
            "result": False,
        },
        {
            "aromatic_matching": False,
            "formal_charge_matching": False,
            "bond_order_matching": False,
            "atom_stereochemistry_matching": False,
            "bond_stereochemistry_matching": False,
            "result": True,
        },
        {
            "aromatic_matching": False,
            "formal_charge_matching": True,
            "bond_order_matching": False,
            "atom_stereochemistry_matching": True,
            "bond_stereochemistry_matching": True,
            "result": True,
        },
        {
            "aromatic_matching": False,
            "formal_charge_matching": False,
            "bond_order_matching": False,
            "atom_stereochemistry_matching": True,
            "bond_stereochemistry_matching": True,
            "result": True,
        },
    ]

    @pytest.mark.parametrize("inputs", isomorphic_permutations)
    def test_isomorphic_perumtations(self, inputs):
        """Test all of the different combinations of matching levels between benzene with and without the aromatic bonds
        defined"""
        # get benzene with all aromatic atoms/bonds labeled
        benzene = Molecule.from_smiles("c1ccccc1")
        # get benzene with no aromatic labels
        benzene_no_aromatic = create_benzene_no_aromatic()
        # now test all of the variations
        assert (
            Molecule.are_isomorphic(
                benzene,
                benzene_no_aromatic,
                aromatic_matching=inputs["aromatic_matching"],
                formal_charge_matching=inputs["formal_charge_matching"],
                bond_order_matching=inputs["bond_order_matching"],
                atom_stereochemistry_matching=inputs["atom_stereochemistry_matching"],
                bond_stereochemistry_matching=inputs["bond_stereochemistry_matching"],
            )[0]
            is inputs["result"]
        )

    @requires_openeye
    def test_strip_atom_stereochemistry(self):
        """Test the basic behavior of strip_atom_stereochemistry"""
        mol = Molecule.from_smiles("CCC[N@@](C)CC")

        nitrogen_idx = [
            atom.molecule_atom_index for atom in mol.atoms if atom.element.symbol == "N"
        ][0]

        # TODO: This fails with RDKitToolkitWrapper because it perceives
        # the stereochemistry of this nitrogen as None
        assert mol.atoms[nitrogen_idx].stereochemistry == "S"
        mol.strip_atom_stereochemistry(smarts="[N+0X3:1](-[*])(-[*])(-[*])")
        assert mol.atoms[nitrogen_idx].stereochemistry is None

        mol = Molecule.from_smiles("CCC[N@@](C)CC")

        assert mol.atoms[nitrogen_idx].stereochemistry == "S"
        mol.strip_atom_stereochemistry(smarts="[N+0X3:1](-[*])(-[*])(-[*])")
        assert mol.atoms[nitrogen_idx].stereochemistry is None

    @pytest.mark.parametrize("mol_smiles", [("C[N+](CC)(CC)CC"), ("c1cnc[nH]1")])
    def test_do_not_strip_atom_stereochemsitry(self, mol_smiles):
        """
        Test special cases in which we do not want nitrogen stereochemistry stripped:
        NH in planar rings and tetra-coordinatined N+
        """
        mol = Molecule.from_smiles(mol_smiles)
        mol_mod = copy.deepcopy(mol)
        mol_mod.strip_atom_stereochemistry("[N+0X3:1](-[*])(-[*])(-[*])")

        # Inspect each atom's stereochemistry instead of relying on __eq__
        for before, after in zip(mol.atoms, mol_mod.atoms):
            assert before.stereochemistry == after.stereochemistry

    def test_isomorphic_stripped_stereochemistry(self):
        """Test that the equality operator disregards an edge case of nitrogen stereocenters"""
        mol1 = Molecule.from_smiles("CCC[N@](C)CC")
        mol2 = Molecule.from_smiles("CCC[N@@](C)CC")

        # Ensure default value is respected and order does not matter
        assert Molecule.are_isomorphic(mol1, mol2, strip_pyrimidal_n_atom_stereo=True)
        assert Molecule.are_isomorphic(mol1, mol2)
        assert Molecule.are_isomorphic(mol2, mol1)

        assert mol1 == mol2
        assert Molecule.from_smiles("CCC[N@](C)CC") == Molecule.from_smiles(
            "CCC[N@@](C)CC"
        )

    def test_remap(self):
        """Test the remap function which should return a new molecule in the requested ordering"""
        # the order here is CCO
        ethanol = create_ethanol()
        # get ethanol in reverse order OCC
        ethanol_reverse = create_reversed_ethanol()
        # get the mapping between the molecules
        mapping = Molecule.are_isomorphic(ethanol, ethanol_reverse, True)[1]

        atom0, atom1 = list([atom for atom in ethanol.atoms])[:2]
        ethanol.add_bond_charge_virtual_site([atom0, atom1], 0.3 * unit.angstrom)

        # make sure that molecules with virtual sites raises an error
        with pytest.raises(NotImplementedError):
            remapped = ethanol.remap(mapping, current_to_new=True)

        # remake with no virtual site and remap to match the reversed ordering
        ethanol = create_ethanol()

        new_ethanol = ethanol.remap(mapping, current_to_new=True)

        def assert_molecules_match_after_remap(mol1, mol2):
            """Check all of the attributes in a molecule match after being remapped"""
            for atoms in zip(mol1.atoms, mol2.atoms):
                assert atoms[0].to_dict() == atoms[1].to_dict()
            # bonds will not be in the same order in the molecule and the atom1 and atom2 indecies could be out of order
            # make a dict to compare them both
            remapped_bonds = dict(
                ((bond.atom1_index, bond.atom2_index), bond) for bond in mol2.bonds
            )
            for bond in mol1.bonds:
                key = (bond.atom1_index, bond.atom2_index)
                if key not in remapped_bonds:
                    key = tuple(reversed(key))
                assert key in remapped_bonds
                # now compare each attribute of the bond except the atom indexes
                bond_dict = bond.to_dict()
                del bond_dict["atom1"]
                del bond_dict["atom2"]
                remapped_bond_dict = remapped_bonds[key].to_dict()
                del remapped_bond_dict["atom1"]
                del remapped_bond_dict["atom2"]
            assert mol1.n_bonds == mol2.n_bonds
            assert mol1.n_angles == mol2.n_angles
            assert mol1.n_propers == mol2.n_propers
            assert mol1.n_impropers == mol2.n_impropers
            assert mol1.total_charge == mol2.total_charge
            assert mol1.partial_charges.all() == mol2.partial_charges.all()

        # check all of the properties match as well, torsions and impropers will be in a different order
        # due to the bonds being out of order
        assert_molecules_match_after_remap(new_ethanol, ethanol_reverse)

        # test round trip (double remapping a molecule)
        new_ethanol = ethanol.remap(mapping, current_to_new=True)
        isomorphic, round_trip_mapping = Molecule.are_isomorphic(
            new_ethanol, ethanol, return_atom_map=True
        )
        assert isomorphic is True
        round_trip_ethanol = new_ethanol.remap(round_trip_mapping, current_to_new=True)
        assert_molecules_match_after_remap(round_trip_ethanol, ethanol)

    @requires_openeye
    def test_canonical_ordering_openeye(self):
        """Make sure molecules are returned in canonical ordering of openeye"""
        from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper

        openeye = OpenEyeToolkitWrapper()
        # get ethanol in canonical order
        ethanol = create_ethanol()
        # get reversed non canonical ethanol
        reversed_ethanol = create_reversed_ethanol()
        # get the canonical ordering
        canonical_ethanol = reversed_ethanol.canonical_order_atoms(openeye)
        # make sure the mapping between the ethanol and the openeye ref canonical form is the same
        assert (
            True,
            {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8},
        ) == Molecule.are_isomorphic(canonical_ethanol, ethanol, True)

    @requires_rdkit
    def test_canonical_ordering_rdkit(self):
        """Make sure molecules are returned in canonical ordering of the RDKit"""
        from openff.toolkit.utils.toolkits import RDKitToolkitWrapper

        rdkit = RDKitToolkitWrapper()
        # get ethanol in canonical order
        ethanol = create_ethanol()
        # get reversed non canonical ethanol
        reversed_ethanol = create_reversed_ethanol()
        # get the canonical ordering
        canonical_ethanol = reversed_ethanol.canonical_order_atoms(rdkit)
        # make sure the mapping between the ethanol and the rdkit ref canonical form is the same
        assert (
            True,
            {0: 2, 1: 0, 2: 1, 3: 8, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7},
        ) == Molecule.are_isomorphic(canonical_ethanol, ethanol, True)

    def test_too_small_remap(self):
        """Make sure remap fails if we do not supply enough indexes"""
        ethanol = Molecule.from_file(get_data_file_path("molecules/ethanol.sdf"))
        # catch mappings that are the wrong size
        too_small_mapping = {0: 1}
        with pytest.raises(ValueError):
            new_ethanol = ethanol.remap(too_small_mapping, current_to_new=True)

    def test_wrong_index_mapping(self):
        """Make sure the remap fails when the indexing starts from the wrong value"""
        ethanol = Molecule.from_file(get_data_file_path("molecules/ethanol.sdf"))
        mapping = {0: 2, 1: 1, 2: 0, 3: 6, 4: 7, 5: 8, 6: 4, 7: 5, 8: 3}
        wrong_index_mapping = dict(
            (i + 10, new_id) for i, new_id in enumerate(mapping.values())
        )
        with pytest.raises(IndexError):
            new_ethanol = ethanol.remap(wrong_index_mapping, current_to_new=True)

    tautomer_data = [
        {"molecule": "Oc1c(cccc3)c3nc2ccncc12", "tautomers": 2},
        {"molecule": "CN=c1nc[nH]cc1", "tautomers": 2},
        {"molecule": "c1[nH]c2c(=O)[nH]c(nc2n1)N", "tautomers": 14},
    ]

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    @pytest.mark.parametrize("molecule_data", tautomer_data)
    def test_enumerating_tautomers(self, molecule_data, toolkit_class):
        """Test the ability of each toolkit to produce tautomers of an input molecule."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = Molecule.from_smiles(
                molecule_data["molecule"],
                allow_undefined_stereo=True,
                toolkit_registry=toolkit,
            )

            tautomers = mol.enumerate_tautomers(toolkit_registry=toolkit)

            assert len(tautomers) == molecule_data["tautomers"]
            assert mol not in tautomers
            # check that the molecules are not isomorphic of the input
            for taut in tautomers:
                assert taut.n_conformers == 0
                assert mol.is_isomorphic_with(taut) is False

        else:
            pytest.skip("Required toolkit is unavailable")

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_enumerating_tautomers_options(self, toolkit_class):
        """Test the enumeration options"""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            # test the max molecules option
            mol = Molecule.from_smiles(
                "c1[nH]c2c(=O)[nH]c(nc2n1)N",
                toolkit_registry=toolkit,
                allow_undefined_stereo=True,
            )

            tauts_no = 5
            tautomers = mol.enumerate_tautomers(
                max_states=tauts_no, toolkit_registry=toolkit
            )
            assert len(tautomers) <= tauts_no
            assert mol not in tautomers

    @pytest.mark.parametrize(
        "toolkit_class", [RDKitToolkitWrapper, OpenEyeToolkitWrapper]
    )
    def test_enumerating_no_tautomers(self, toolkit_class):
        """Test that the toolkits return an empty list if there are no tautomers to enumerate."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = Molecule.from_smiles("CC", toolkit_registry=toolkit)

            tautomers = mol.enumerate_tautomers(toolkit_registry=toolkit)
            assert tautomers == []

        else:
            pytest.skip("Required toolkit is unavailable")

    @requires_openeye
    def test_enumerating_no_protomers(self):
        """Make sure no protomers are returned."""

        mol = Molecule.from_smiles("CC")

        assert mol.enumerate_protomers() == []

    @requires_openeye
    def test_enumerating_protomers(self):
        """Test enumerating the formal charges."""

        mol = Molecule.from_smiles("Oc2ccc(c1ccncc1)cc2")

        # there should be three protomers for this molecule so restrict the output
        protomers = mol.enumerate_protomers(max_states=2)

        assert mol not in protomers
        assert len(protomers) == 2

        # now make sure we can generate them all
        protomers = mol.enumerate_protomers(max_states=10)

        assert mol not in protomers
        assert len(protomers) == 3

        # make sure each protomer is unique
        unique_protomers = set(protomers)
        assert len(protomers) == len(unique_protomers)

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_enumerating_stereobonds(self, toolkit_class):
        """Test the backend toolkits in enumerating the stereo bonds in a molecule."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = Molecule.from_smiles(
                "ClC=CCl", allow_undefined_stereo=True, toolkit_registry=toolkit
            )

            # use the default options
            isomers = mol.enumerate_stereoisomers()
            assert len(isomers) == 2

            assert mol not in isomers
            # make sure the input molecule is only different by bond stereo
            for ismol in isomers:
                assert (
                    Molecule.are_isomorphic(
                        mol,
                        ismol,
                        return_atom_map=False,
                        bond_stereochemistry_matching=False,
                    )[0]
                    is True
                )
                assert mol.is_isomorphic_with(ismol) is False

            # make sure the isomers are different
            assert isomers[0].is_isomorphic_with(isomers[1]) is False

        else:
            pytest.skip("Required toolkit is unavailable")

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_enumerating_stereocenters(self, toolkit_class):
        """Test the backend toolkits in enumerating the stereo centers in a molecule."""

        if toolkit_class.is_available():
            toolkit = toolkit_class()
            mol = Molecule.from_smiles(
                "NC(Cl)(F)O", toolkit_registry=toolkit, allow_undefined_stereo=True
            )

            isomers = mol.enumerate_stereoisomers(toolkit_registry=toolkit)

            assert len(isomers) == 2
            # make sure the mol is not in the isomers and that they only differ by stereo chem
            assert mol not in isomers
            for ismol in isomers:
                assert ismol.n_conformers != 0
                assert (
                    Molecule.are_isomorphic(
                        mol,
                        ismol,
                        return_atom_map=False,
                        atom_stereochemistry_matching=False,
                    )[0]
                    is True
                )
                assert mol.is_isomorphic_with(ismol) is False

            # make sure the two isomers are different
            assert isomers[0].is_isomorphic_with(isomers[1]) is False

        else:
            pytest.skip("Required toolkit is unavailable")

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_enumerating_stereo_options(self, toolkit_class):
        """Test the enumerating stereo chem options"""

        if toolkit_class.is_available():
            toolkit = toolkit_class()

            # test undefined only
            mol = Molecule.from_smiles(
                "ClC=CCl", toolkit_registry=toolkit, allow_undefined_stereo=True
            )
            isomers = mol.enumerate_stereoisomers(
                undefined_only=True, rationalise=False
            )

            assert len(isomers) == 2
            for isomer in isomers:
                assert isomer.n_conformers == 0

            mol = Molecule.from_smiles(
                "Cl/C=C\Cl", toolkit_registry=toolkit, allow_undefined_stereo=True
            )
            isomers = mol.enumerate_stereoisomers(
                undefined_only=True, rationalise=False
            )

            assert isomers == []

            mol = Molecule.from_smiles(
                "Cl/C=C\Cl", toolkit_registry=toolkit, allow_undefined_stereo=True
            )
            isomers = mol.enumerate_stereoisomers(
                undefined_only=False, rationalise=False
            )

            assert len(isomers) == 1

            # test max isomers
            mol = Molecule.from_smiles(
                "BrC=C[C@H]1OC(C2)(F)C2(Cl)C1",
                toolkit_registry=toolkit,
                allow_undefined_stereo=True,
            )
            isomers = mol.enumerate_stereoisomers(
                max_isomers=5,
                undefined_only=True,
                toolkit_registry=toolkit,
                rationalise=True,
            )

            assert len(isomers) <= 5
            for isomer in isomers:
                assert isomer.n_conformers == 1

        else:
            pytest.skip("Required toolkit is unavailable")

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    @pytest.mark.parametrize(
        "smiles, undefined_only, expected",
        [
            (
                "FC(Br)(Cl)[C@@H](Br)(Cl)",
                False,
                [
                    "F[C@](Br)(Cl)[C@@H](Br)(Cl)",
                    "F[C@](Br)(Cl)[C@H](Br)(Cl)",
                    "F[C@@](Br)(Cl)[C@@H](Br)(Cl)",
                    "F[C@@](Br)(Cl)[C@H](Br)(Cl)",
                ],
            ),
            (
                "FC(Br)(Cl)[C@@H](Br)(Cl)",
                True,
                ["F[C@](Br)(Cl)[C@@H](Br)(Cl)", "F[C@@](Br)(Cl)[C@@H](Br)(Cl)"],
            ),
            ("F[C@H](Cl)Br", False, ["F[C@@H](Cl)Br"]),
            ("F[C@H](Cl)Br", True, []),
        ],
    )
    def test_enumerating_stereo_partially_defined(
        self, toolkit_class, smiles, undefined_only, expected
    ):
        """Test the enumerating stereo of molecules with partially defined chirality"""

        if not toolkit_class.is_available():
            pytest.skip("Required toolkit is unavailable")

        toolkit = toolkit_class()

        # test undefined only
        mol = Molecule.from_smiles(
            smiles, toolkit_registry=toolkit, allow_undefined_stereo=True
        )
        stereoisomers = mol.enumerate_stereoisomers(
            undefined_only=undefined_only, rationalise=False
        )

        # Ensure that the results of the enumeration are what the test expects.
        # This roundtrips the expected output from SMILES --> OFFMol --> SMILES,
        # since the SMILES for stereoisomers generated in this test may change depending
        # on which cheminformatics toolkit is used.
        expected = {
            Molecule.from_smiles(stereoisomer, allow_undefined_stereo=True).to_smiles(
                explicit_hydrogens=True, isomeric=True, mapped=False
            )
            for stereoisomer in expected
        }
        actual = {
            stereoisomer.to_smiles(explicit_hydrogens=True, isomeric=True, mapped=False)
            for stereoisomer in stereoisomers
        }

        assert expected == actual

    @requires_rdkit
    def test_from_pdb_and_smiles(self):
        """Test the ability to make a valid molecule using RDKit and SMILES together"""
        # try and make a molecule from a pdb and smiles that don't match
        with pytest.raises(InvalidConformerError):
            mol = Molecule.from_pdb_and_smiles(
                get_data_file_path("molecules/toluene.pdb"), "CC"
            )

        # make a molecule from the toluene pdb file and the correct smiles
        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("molecules/toluene.pdb"), "Cc1ccccc1"
        )

        # make toluene from the sdf file
        mol_sdf = Molecule.from_file(get_data_file_path("molecules/toluene.sdf"))
        # get the mapping between them and compare the properties
        isomorphic, atom_map = Molecule.are_isomorphic(
            mol, mol_sdf, return_atom_map=True
        )
        assert isomorphic is True
        for pdb_atom, sdf_atom in atom_map.items():
            assert mol.atoms[pdb_atom].to_dict() == mol_sdf.atoms[sdf_atom].to_dict()
        # check bonds match, however there order might not
        sdf_bonds = dict(
            ((bond.atom1_index, bond.atom2_index), bond) for bond in mol_sdf.bonds
        )
        for bond in mol.bonds:
            key = (atom_map[bond.atom1_index], atom_map[bond.atom2_index])
            if key not in sdf_bonds:
                key = tuple(reversed(key))
            assert key in sdf_bonds
            # now compare the attributes
            assert bond.is_aromatic == sdf_bonds[key].is_aromatic
            assert bond.stereochemistry == sdf_bonds[key].stereochemistry

    @requires_pkg("qcportal")
    def test_to_qcschema(self):
        """Test the ability to make and validate qcschema with extras"""
        # the molecule has no coordinates so this should fail
        ethanol = Molecule.from_smiles("CCO")
        with pytest.raises(InvalidConformerError):
            qcschema = ethanol.to_qcschema()

        # now remake the molecule from the sdf
        ethanol = Molecule.from_file(get_data_file_path("molecules/ethanol.sdf"))
        # make sure that requests to missing conformers are caught
        with pytest.raises(InvalidConformerError):
            qcschema = ethanol.to_qcschema(conformer=1)
        # now make a valid qcschema and check its properties
        qcschema = ethanol.to_qcschema(extras={"test_tag": "test"})
        # make sure the properties match
        charge = 0
        connectivity = [
            (0, 1, 1.0),
            (0, 4, 1.0),
            (0, 5, 1.0),
            (0, 6, 1.0),
            (1, 2, 1.0),
            (1, 7, 1.0),
            (1, 8, 1.0),
            (2, 3, 1.0),
        ]
        symbols = ["C", "C", "O", "H", "H", "H", "H", "H", "H"]

        def assert_check():
            assert charge == qcschema.molecular_charge
            assert connectivity == qcschema.connectivity
            assert symbols == qcschema.symbols.tolist()
            assert (
                qcschema.geometry.all()
                == ethanol.conformers[0].in_units_of(unit.bohr).all()
            )

        assert_check()
        assert qcschema.extras["test_tag"] == "test"
        assert qcschema.extras[
            "canonical_isomeric_explicit_hydrogen_mapped_smiles"
        ] == ethanol.to_smiles(mapped=True)
        # # now run again with no extras passed, only cmiles entry will be present with fix-720
        qcschema = ethanol.to_qcschema()
        assert_check()
        assert qcschema.extras[
            "canonical_isomeric_explicit_hydrogen_mapped_smiles"
        ] == ethanol.to_smiles(mapped=True)

    @requires_pkg("qcportal")
    def test_from_qcschema_no_client(self):
        """Test the ability to make molecules from QCArchive record instances and dicts"""

        import json

        # As the method can take a record instance or a dict with JSON encoding test both
        # test incomplete dict
        example_dict = {"name": "CH4"}
        with pytest.raises(KeyError):
            mol = Molecule.from_qcschema(example_dict)

        # test an object that is not a record
        wrong_object = "CH4"
        with pytest.raises(AttributeError):
            mol = Molecule.from_qcschema(wrong_object)

        with open(get_data_file_path("molecules/qcportal_molecules.json")) as json_file:
            # test loading the dict representation from a json file
            json_mol = json.load(json_file)
            mol = Molecule.from_qcschema(json_mol)
            # now make a molecule from the canonical smiles and make sure they are isomorphic
            can_mol = Molecule.from_smiles(
                json_mol["attributes"]["canonical_isomeric_smiles"]
            )
            assert mol.is_isomorphic_with(can_mol) is True

    client_examples = [
        {
            "dataset": "TorsionDriveDataset",
            "name": "Fragment Stability Benchmark",
            "index": "CC(=O)Nc1cc2c(cc1OC)nc[n:4][c:3]2[NH:2][c:1]3ccc(c(c3)Cl)F",
        },
        {
            "dataset": "TorsionDriveDataset",
            "name": "OpenFF Fragmenter Phenyl Benchmark",
            "index": "c1c[ch:1][c:2](cc1)[c:3](=[o:4])o",
        },
        {
            "dataset": "TorsionDriveDataset",
            "name": "OpenFF Full TorsionDrive Benchmark 1",
            "index": "0",
        },
        {
            "dataset": "TorsionDriveDataset",
            "name": "OpenFF Group1 Torsions",
            "index": "c1c[ch:1][c:2](cc1)[ch2:3][c:4]2ccccc2",
        },
        {
            "dataset": "OptimizationDataset",
            "name": "Kinase Inhibitors: WBO Distributions",
            "index": "cc1ccc(cc1nc2nccc(n2)c3cccnc3)nc(=o)c4ccc(cc4)cn5ccn(cc5)c-0",
        },
        {
            "dataset": "OptimizationDataset",
            "name": "SMIRNOFF Coverage Set 1",
            "index": "coc(o)oc-0",
        },
        {
            "dataset": "GridOptimizationDataset",
            "name": "OpenFF Trivalent Nitrogen Set 1",
            "index": "b1(c2c(ccs2)-c3ccsc3n1)c4c(c(c(c(c4f)f)f)f)f",
        },
        {
            "dataset": "GridOptimizationDataset",
            "name": "OpenFF Trivalent Nitrogen Set 1",
            "index": "C(#N)N",
        },
    ]

    @requires_pkg("qcportal")
    @pytest.mark.parametrize("input_data", client_examples)
    def test_from_qcschema_with_client(self, input_data):
        """For each of the examples try and make a offmol using the instance and dict and check they match"""

        import qcportal as ptl

        client = ptl.FractalClient()
        ds = client.get_collection(input_data["dataset"], input_data["name"])
        entry = ds.get_entry(input_data["index"])
        # now make the molecule from the record instance with and without the geometry
        mol_from_dict = Molecule.from_qcschema(entry.dict(encoding="json"))
        # make the molecule again with the geometries attached
        mol_from_instance = Molecule.from_qcschema(entry, client)
        if hasattr(entry, "initial_molecules"):
            assert mol_from_instance.n_conformers == len(entry.initial_molecules)
        else:
            # opt records have one initial molecule
            assert mol_from_instance.n_conformers == 1

        # now make a molecule from the smiles and make sure they are isomorphic
        mol_from_smiles = Molecule.from_smiles(
            entry.attributes["canonical_explicit_hydrogen_smiles"], True
        )

        assert mol_from_dict.is_isomorphic_with(mol_from_smiles) is True

    @requires_pkg("qcportal")
    def test_qcschema_round_trip(self):
        """Test making a molecule from qcschema then converting back
        Checking whether qca_mol and mol created from/to qcschema are the same or not"""

        # get a molecule qcschema
        import qcportal as ptl

        client = ptl.FractalClient()
        ds = client.get_collection("OptimizationDataset", "SMIRNOFF Coverage Set 1")
        # grab an entry from the optimization data set
        entry = ds.get_entry("coc(o)oc-0")
        # now make the molecule from the record instance with the geometry
        mol = Molecule.from_qcschema(entry, client)
        # now grab the initial molecule record
        qca_mol = client.query_molecules(id=entry.initial_molecule)[0]
        # mow make sure the majority of the qcschema attributes are the same
        # note we can not compare the full dict due to qcelemental differences
        qcschema = mol.to_qcschema()
        assert qcschema.atom_labels.tolist() == qca_mol.atom_labels.tolist()
        assert qcschema.symbols.tolist() == qca_mol.symbols.tolist()
        # due to conversion using different programs there is a slight difference here
        assert qcschema.geometry.flatten().tolist() == pytest.approx(
            qca_mol.geometry.flatten().tolist(), rel=1.0e-5
        )
        assert qcschema.connectivity == qca_mol.connectivity
        assert qcschema.atomic_numbers.tolist() == qca_mol.atomic_numbers.tolist()
        assert qcschema.fragment_charges == qca_mol.fragment_charges
        assert qcschema.fragment_multiplicities == qca_mol.fragment_multiplicities
        assert qcschema.fragments[0].tolist() == qca_mol.fragments[0].tolist()
        assert qcschema.mass_numbers.tolist() == qca_mol.mass_numbers.tolist()
        assert qcschema.name == qca_mol.name
        assert qcschema.masses.all() == qca_mol.masses.all()
        assert qcschema.molecular_charge == qca_mol.molecular_charge
        assert qcschema.molecular_multiplicity == qca_mol.molecular_multiplicity
        assert qcschema.real.all() == qca_mol.real.all()

    @requires_pkg("qcportal")
    def test_qcschema_round_trip_from_to_from(self):
        """Test making a molecule from qca record using from_qcschema,
        then converting back to qcschema using to_qcschema,
         and then reading that again using from_qcschema"""

        # get a molecule qcschema
        import qcportal as ptl

        client = ptl.FractalClient()
        ds = client.get_collection(
            "TorsionDriveDataset", "OpenFF-benchmark-ligand-fragments-v1.0"
        )
        # grab an entry from the torsiondrive data set
        entry = ds.get_entry(
            "[H]c1[c:1]([c:2](c(c(c1[H])N([H])C(=O)[H])[H])[C:3]2=C(C(=C([S:4]2)[H])OC([H])([H])[H])Br)[H]"
        )
        # now make the molecule from the record instance with the geometry
        mol_qca_record = Molecule.from_qcschema(entry, client)
        off_qcschema = mol_qca_record.to_qcschema()
        mol_using_from_off_qcschema = Molecule.from_qcschema(off_qcschema)
        assert_molecule_is_equal(
            mol_qca_record,
            mol_using_from_off_qcschema,
            "Molecule roundtrip to/from_qcschema failed",
        )

    @requires_pkg("qcportal")
    def test_qcschema_round_trip_raise_error(self):
        """Test making a molecule from qcschema,
        reaching inner except block where everything fails"""

        # get a molecule qcschema
        import qcportal as ptl

        client = ptl.FractalClient()
        ds = client.get_collection(
            "TorsionDriveDataset", "OpenFF-benchmark-ligand-fragments-v1.0"
        )
        # grab an entry from the torsiondrive data set
        entry = ds.get_entry(
            "[H]c1[c:1]([c:2](c(c(c1[H])N([H])C(=O)[H])[H])[C:3]2=C(C(=C([S:4]2)[H])OC([H])([H])[H])Br)[H]"
        )
        del entry.attributes["canonical_isomeric_explicit_hydrogen_mapped_smiles"]
        # now make the molecule from the record instance with the geometry
        with pytest.raises(KeyError):
            mol_qca_record = Molecule.from_qcschema(entry, client)

    @requires_pkg("qcportal")
    def test_qcschema_molecule_record_round_trip_from_to_from(self):
        """Test making a molecule from qca record using from_qcschema,
        then converting back to qcschema using to_qcschema,
         and then reading that again using from_qcschema"""

        # get a molecule qcschema
        import qcportal as ptl

        client = ptl.FractalClient()

        record = client.query_molecules(molecular_formula="C16H20N3O5")[0]

        # now make the molecule from the record instance with the geometry
        mol_qca_record = Molecule.from_qcschema(record, client)
        off_qcschema = mol_qca_record.to_qcschema()
        mol_using_from_off_qcschema = Molecule.from_qcschema(off_qcschema)

        assert_molecule_is_equal(
            mol_qca_record,
            mol_using_from_off_qcschema,
            "Molecule roundtrip to/from_qcschema failed",
        )

    def test_from_mapped_smiles(self):
        """Test making the molecule from issue #412 using both toolkits to ensure the issue
        is fixed."""

        # there should be no undefined sterochmeistry error when making the molecule
        mol = Molecule.from_mapped_smiles(
            "[H:14][c:1]1[c:3]([c:7]([c:11]([c:8]([c:4]1[H:17])[H:21])[C:13]([H:24])([H:25])[c:12]2[c:9]([c:5]([c:2]([c:6]([c:10]2[H:23])[H:19])[H:15])[H:18])[H:22])[H:20])[H:16]"
        )
        assert mol.n_atoms == 25
        # make sure the atom map is not exposed
        with pytest.raises(KeyError):
            mapping = mol._properties["atom_map"]

    @pytest.mark.parametrize(
        "toolkit_class", [OpenEyeToolkitWrapper, RDKitToolkitWrapper]
    )
    def test_from_mapped_smiles_partial(self, toolkit_class):
        """Test that creating a molecule from a partially mapped SMILES raises an
        exception."""
        if not (toolkit_class.is_available()):
            pytest.skip(f"Required toolkit {toolkit_class} is unavailable")
        with pytest.raises(
            SmilesParsingError,
            match="The mapped smiles does not contain enough indexes",
        ):
            Molecule.from_mapped_smiles("[Cl:1][Cl]", toolkit_registry=toolkit_class())

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_n_particles(self, molecule):
        """Test n_particles property"""
        n_particles = sum([1 for particle in molecule.particles])
        assert n_particles == molecule.n_particles

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_n_atoms(self, molecule):
        """Test n_atoms property"""
        n_atoms = sum([1 for atom in molecule.atoms])
        assert n_atoms == molecule.n_atoms

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_n_virtual_sites(self, molecule):
        """Test n_virtual_sites property"""
        n_virtual_sites = sum([1 for virtual_site in molecule.virtual_sites])
        assert n_virtual_sites == molecule.n_virtual_sites

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_n_bonds(self, molecule):
        """Test n_bonds property"""
        n_bonds = sum([1 for bond in molecule.bonds])
        assert n_bonds == molecule.n_bonds

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_angles(self, molecule):
        """Test angles property"""
        for angle in molecule.angles:
            assert angle[0].is_bonded_to(angle[1])
            assert angle[1].is_bonded_to(angle[2])

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_propers(self, molecule):
        """Test propers property"""
        for proper in molecule.propers:
            # The bonds should be in order 0-1-2-3 unless the
            # atoms form a three- or four-membered ring.
            is_chain = proper[0].is_bonded_to(proper[1])
            is_chain &= proper[1].is_bonded_to(proper[2])
            is_chain &= proper[2].is_bonded_to(proper[3])
            is_chain &= not proper[0].is_bonded_to(proper[2])
            is_chain &= not proper[0].is_bonded_to(proper[3])
            is_chain &= not proper[1].is_bonded_to(proper[3])

            assert (
                is_chain
                or is_three_membered_ring_torsion(proper)
                or is_four_membered_ring_torsion(proper)
            )

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_impropers(self, molecule):
        """Test impropers property"""
        for improper in molecule.impropers:
            assert improper[0].is_bonded_to(improper[1])
            assert improper[1].is_bonded_to(improper[2])
            assert improper[1].is_bonded_to(improper[3])

            # The non-central atoms can be connected only if
            # the improper atoms form a three-membered ring.
            is_not_cyclic = not (
                (improper[0].is_bonded_to(improper[2]))
                or (improper[0].is_bonded_to(improper[3]))
                or (improper[2].is_bonded_to(improper[3]))
            )
            assert is_not_cyclic or is_three_membered_ring_torsion(improper)

    @pytest.mark.parametrize(
        ("molecule", "n_impropers", "n_pruned"),
        [
            ("C", 24, 0),
            ("CC", 48, 0),
            ("N", 6, 6),
        ],
    )
    def test_pruned_impropers(self, molecule, n_impropers, n_pruned):
        """Test the amber_impropers and smirnoff_impropers properties"""
        mol = Molecule.from_smiles(molecule)
        assert mol.n_impropers == n_impropers
        assert len(mol.smirnoff_impropers) == n_pruned
        assert len(mol.amber_impropers) == n_pruned

        # Order not guaranteed, so cannot zip and compare directly
        for smirnoff_imp in mol.smirnoff_impropers:
            # Convert SMIRNOFF-style improper into AMBER-style
            mod_imp = (
                smirnoff_imp[1],
                smirnoff_imp[0],
                smirnoff_imp[2],
                smirnoff_imp[3],
            )
            assert mod_imp in mol.amber_impropers

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_torsions(self, molecule):
        """Test torsions property"""
        # molecule.torsions should be exactly equal to the union of propers and impropers.
        assert set(molecule.torsions) == set(molecule.propers) | set(molecule.impropers)

        # The intersection of molecule.propers and molecule.impropers should be largely null.
        # The only exception is for molecules containing 3-membered rings (e.g., DrugBank_5514).
        common_torsions = molecule.propers & molecule.impropers
        if len(common_torsions) > 0:
            for torsion in common_torsions:
                assert is_three_membered_ring_torsion(torsion)

    def test_nth_degree_neighbors_basic(self):
        benzene = Molecule.from_smiles("c1ccccc1")

        with pytest.raises(ValueError, match="0 or fewer atoms"):
            benzene.nth_degree_neighbors(n_degrees=0)

        with pytest.raises(ValueError, match="0 or fewer atoms"):
            benzene.nth_degree_neighbors(n_degrees=-1)

        # 3 proper torsions have 1-4 pairs that are duplicates of other torsions, i.e.
        # torsions 0-1-2-3 0-5-4-3 where #s are the indices of carbons in the ring
        n_14_pairs = len([*benzene.nth_degree_neighbors(n_degrees=3)])
        assert n_14_pairs == benzene.n_propers - 3

        for pair in benzene.nth_degree_neighbors(n_degrees=3):
            assert pair[0].molecule_particle_index < pair[1].molecule_particle_index

    @pytest.mark.parametrize(
        ("smiles", "n_degrees", "num_pairs"),
        [
            ("c1ccccc1", 6, 0),
            ("c1ccccc1", 5, 3),
            ("c1ccccc1", 4, 12),
            ("c1ccccc1", 3, 21),
            ("N1ONON1", 4, 2),
            ("N1ONON1", 3, 7),
            ("C1#CO1", 2, 0),
            ("C1#CO1", 1, 3),
            ("C1CO1", 4, 0),
            ("C1CO1", 2, 10),
            ("C1=C=C=C=1", 4, 0),
            ("C1=C=C=C=1", 3, 0),
            ("C1=C=C=C=1", 2, 2),
        ],
    )
    def test_nth_degree_neighbors_rings(self, smiles, n_degrees, num_pairs):
        molecule = Molecule.from_smiles(smiles)

        num_pairs_found = len([*molecule.nth_degree_neighbors(n_degrees=n_degrees)])
        assert num_pairs_found == num_pairs

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_total_charge(self, molecule):
        """Test total charge"""
        charge_sum = 0 * unit.elementary_charge
        for atom in molecule.atoms:
            charge_sum += atom.formal_charge
        assert charge_sum == molecule.total_charge

    # ----------------------------------------------------
    # Test magic methods.
    # ----------------------------------------------------

    def test_equality(self):
        """Test equality operator"""
        molecules = mini_drug_bank()
        nmolecules = len(molecules)
        # TODO: Performance improvements should let us un-restrict this test
        for i in range(nmolecules):
            for j in range(i, min(i + 3, nmolecules)):
                assert (molecules[i] == molecules[j]) == (i == j)

    # ----------------------
    # Test Molecule methods.
    # ----------------------

    def test_add_conformers(self):
        """Test addition of conformers to a molecule"""
        import numpy as np
        from simtk import unit

        # Define a methane molecule
        molecule = Molecule()
        molecule.name = "methane"
        C = molecule.add_atom(6, 0, False)
        H1 = molecule.add_atom(1, 0, False)
        H2 = molecule.add_atom(1, 0, False)
        H3 = molecule.add_atom(1, 0, False)
        H4 = molecule.add_atom(1, 0, False)
        molecule.add_bond(C, H1, 1, False)
        molecule.add_bond(C, H2, 1, False)
        molecule.add_bond(C, H3, 1, False)
        molecule.add_bond(C, H4, 1, False)

        assert molecule.n_conformers == 0
        # Add a conformer that should work
        conf1 = unit.Quantity(
            np.array(
                [
                    [1.0, 2.0, 3.0],
                    [4.0, 5.0, 6.0],
                    [7.0, 8.0, 9.0],
                    [10.0, 11.0, 12.0],
                    [13.0, 14.0, 15],
                ]
            ),
            unit.angstrom,
        )
        molecule.add_conformer(conf1)
        assert molecule.n_conformers == 1

        conf2 = unit.Quantity(
            np.array(
                [
                    [101.0, 102.0, 103.0],
                    [104.0, 105.0, 106.0],
                    [107.0, 108.0, 109.0],
                    [110.0, 111.0, 112.0],
                    [113.0, 114.0, 115],
                ]
            ),
            unit.angstrom,
        )
        molecule.add_conformer(conf2)
        assert molecule.n_conformers == 2

        # Add conformers with too few coordinates
        conf_missing_z = unit.Quantity(
            np.array(
                [
                    [101.0, 102.0, 103.0],
                    [104.0, 105.0, 106.0],
                    [107.0, 108.0, 109.0],
                    [110.0, 111.0, 112.0],
                    [113.0, 114.0],
                ]
            ),
            unit.angstrom,
        )
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_missing_z)

        conf_too_few_atoms = unit.Quantity(
            np.array(
                [
                    [101.0, 102.0, 103.0],
                    [104.0, 105.0, 106.0],
                    [107.0, 108.0, 109.0],
                    [110.0, 111.0, 112.0],
                ]
            ),
            unit.angstrom,
        )
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_too_few_atoms)

        # Add a conformer with too many coordinates
        conf_too_many_atoms = unit.Quantity(
            np.array(
                [
                    [101.0, 102.0, 103.0],
                    [104.0, 105.0, 106.0],
                    [107.0, 108.0, 109.0],
                    [110.0, 111.0, 112.0],
                    [113.0, 114.0, 115.0],
                    [116.0, 117.0, 118.0],
                ]
            ),
            unit.angstrom,
        )
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_too_many_atoms)

        # Add a conformer with no coordinates
        conf_no_coordinates = unit.Quantity(np.array([]), unit.angstrom)
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_no_coordinates)

        # Add a conforer with units of nanometers
        conf3 = unit.Quantity(
            np.array(
                [
                    [1.0, 2.0, 3.0],
                    [4.0, 5.0, 6.0],
                    [7.0, 8.0, 9.0],
                    [10.0, 11.0, 12.0],
                    [13.0, 14.0, 15],
                ]
            ),
            unit.nanometer,
        )
        molecule.add_conformer(conf3)
        assert molecule.n_conformers == 3
        assert molecule.conformers[2][0][0] == 10.0 * unit.angstrom

        # Add a conformer with units of nanometers
        conf_nonsense_units = unit.Quantity(
            np.array(
                [
                    [1.0, 2.0, 3.0],
                    [4.0, 5.0, 6.0],
                    [7.0, 8.0, 9.0],
                    [10.0, 11.0, 12.0],
                    [13.0, 14.0, 15],
                ]
            ),
            unit.joule,
        )
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_nonsense_units)

        # Add a conformer with no units
        conf_unitless = np.array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15],
            ]
        )
        with pytest.raises(Exception) as excinfo:
            molecule.add_conformer(conf_unitless)

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_add_atoms_and_bonds(self, molecule):
        """Test the creation of a molecule from the addition of atoms and bonds"""
        molecule_copy = Molecule()
        for atom in molecule.atoms:
            molecule_copy.add_atom(
                atom.atomic_number,
                atom.formal_charge,
                atom.is_aromatic,
                stereochemistry=atom.stereochemistry,
            )
        for bond in molecule.bonds:
            molecule_copy.add_bond(
                bond.atom1_index,
                bond.atom2_index,
                bond.bond_order,
                bond.is_aromatic,
                stereochemistry=bond.stereochemistry,
                fractional_bond_order=bond.fractional_bond_order,
            )
        # Try to add the final bond twice, which should raise an Exception
        with pytest.raises(Exception) as excinfo:
            molecule_copy.add_bond(
                bond.atom1_index,
                bond.atom2_index,
                bond.bond_order,
                bond.is_aromatic,
                stereochemistry=bond.stereochemistry,
                fractional_bond_order=bond.fractional_bond_order,
            )

        assert molecule == molecule_copy

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_add_virtual_site_units(self, molecule):
        """
        Tests the unit type checking of the VirtualSite base class
        """

        # TODO: Should these be using BondChargeVirtualSite, or should we just call the base class (which does the unit checks) directly?

        # Prepare values for unit checks
        distance_unitless = 0.4
        sigma_unitless = 0.1
        rmin_half_unitless = 0.2
        epsilon_unitless = 0.3
        charge_increments_unitless = [0.1, 0.2]
        distance = distance_unitless * unit.angstrom
        sigma = sigma_unitless * unit.angstrom
        rmin_half = rmin_half_unitless * unit.angstrom
        epsilon = epsilon_unitless * (unit.kilojoule / unit.mole)
        charge_increments = charge_increments_unitless * unit.elementary_charge

        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]

        atoms = (atom1, atom2)

        # Try to feed in unitless sigma
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site(
                atoms, distance, epsilon=epsilon, sigma=sigma_unitless, replace=True
            )

        # Try to feed in unitless rmin_half
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site(
                atoms,
                distance,
                epsilon=epsilon,
                rmin_half=rmin_half_unitless,
                replace=True,
            )

        # Try to feed in unitless epsilon
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site(
                atoms,
                distance,
                epsilon=epsilon_unitless,
                sigma=sigma,
                rmin_half=rmin_half,
                replace=True,
            )

        # Try to feed in unitless charges
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site(
                atoms,
                distance,
                charge_incrtements=charge_increments_unitless,
                replace=True,
            )

        # We shouldn't be able to give both rmin_half and sigma VdW parameters.
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site(
                [atom1, atom2],
                distance,
                epsilon=epsilon,
                sigma=sigma,
                rmin_half=rmin_half,
                replace=True,
            )

        # Try creating virtual site from sigma+epsilon; replace=False since this
        # should be the first vsite to succeed
        vsite1_index = molecule.add_bond_charge_virtual_site(
            atoms, distance, epsilon=epsilon, sigma=sigma, replace=False
        )
        # Try creating virutal site from rmin_half+epsilon
        vsite2_index = molecule.add_bond_charge_virtual_site(
            atoms, distance, epsilon=epsilon, rmin_half=rmin_half, replace=True
        )

        # TODO: Test the @property getters for sigma, epsilon, and rmin_half

        # We should have to give as many charge increments as atoms (len(charge_increments)) = 4
        with pytest.raises(Exception) as excinfo:
            molecule.add_bond_charge_virtual_site(
                atoms, distance, charge_increments=[0.0], replace=True
            )

        vsite3_index = molecule.add_bond_charge_virtual_site(
            atoms, distance, charge_increments=charge_increments, replace=True
        )

    def test_virtual_particle(self):
        """Test setter, getters, and methods of VirtualParticle"""
        mol = create_ethanol()
        # Add a SYMMETRIC (default) VirtualSite, which should have two particles
        mol.add_bond_charge_virtual_site(
            (mol.atoms[1], mol.atoms[2]),
            0.5 * unit.angstrom,
            charge_increments=[-0.1, 0.1] * unit.elementary_charge,
        )
        # Add a NON SYMMETRIC (specified in kwarg) VirtualSite, which should have one particle
        mol.add_bond_charge_virtual_site(
            (mol.atoms[0], mol.atoms[1]),
            0.5 * unit.angstrom,
            charge_increments=[-0.1, 0.1] * unit.elementary_charge,
            symmetric=False,
        )

        assert len(mol.virtual_sites) == 2

        # Ensure the first vsite has two particles
        vps0 = [vp for vp in mol.virtual_sites[0].particles]
        assert len(vps0) == 2
        assert vps0[0].virtual_site_particle_index == 0
        assert vps0[1].virtual_site_particle_index == 1
        orientations0 = [vp.orientation for vp in vps0]
        assert len(set(orientations0) & {(0, 1), (1, 0)}) == 2

        # Ensure the second vsite has one particle
        vps1 = [vp for vp in mol.virtual_sites[1].particles]
        assert len(vps1) == 1
        assert vps1[0].virtual_site_particle_index == 0
        orientations1 = [vp.orientation for vp in vps1]
        assert len(set(orientations1) & {(0, 1)}) == 1

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_add_bond_charge_virtual_site(self, molecule):
        """Test the addition of a BondChargeVirtualSite to a molecule.
        Also tests many of the inputs of the parent VirtualSite class
        """
        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]

        # Prepare values for unit checks
        distance_unitless = 0.4
        distance = distance_unitless * unit.angstrom

        # Try to feed in a unitless distance
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_bond_charge_virtual_site(
                [atom1, atom2], distance_unitless
            )

        vsite1_index = molecule.add_bond_charge_virtual_site([atom1, atom2], distance)
        vsite1 = molecule.virtual_sites[vsite1_index]
        assert atom1 in vsite1.atoms
        assert atom2 in vsite1.atoms
        assert vsite1 in atom1.virtual_sites
        assert vsite1 in atom2.virtual_sites
        assert vsite1.distance == distance

        # Make a virtual site using all arguments
        vsite2_index = molecule.add_bond_charge_virtual_site(
            [atom2, atom3],
            distance,
            sigma=0.1 * unit.angstrom,
            epsilon=1.0 * unit.kilojoule_per_mole,
            charge_increments=unit.Quantity(
                np.array([0.1, 0.2]), unit.elementary_charge
            ),
        )
        vsite2 = molecule.virtual_sites[vsite2_index]
        assert atom2 in vsite2.atoms
        assert atom3 in vsite2.atoms
        assert vsite2 in atom2.virtual_sites
        assert vsite2 in atom3.virtual_sites
        assert vsite2.distance == distance

        # test serialization
        molecule_dict = molecule.to_dict()
        molecule2 = Molecule.from_dict(molecule_dict)

        # test that the molecules are the same
        assert molecule.is_isomorphic_with(molecule2)
        # lets also make sure that the vsites were constructed correctly
        for site1, site2 in zip(molecule.virtual_sites, molecule2.virtual_sites):
            assert site1.to_dict() == site2.to_dict()

    # TODO: Make a test for to_dict and from_dict for VirtualSites (even though they're currently just unloaded using
    #      (for example) Molecule._add_bond_virtual_site functions
    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_add_monovalent_lone_pair_virtual_site(self, molecule):
        """Test addition of a MonovalentLonePairVirtualSite to the Molecule"""
        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]

        # Prepare values for unit checks
        distance_unitless = 0.3
        out_of_plane_angle_unitless = 30
        in_plane_angle_unitless = 0.2
        distance = distance_unitless * unit.angstrom
        out_of_plane_angle = out_of_plane_angle_unitless * unit.degree
        in_plane_angle = in_plane_angle_unitless * unit.radian

        atoms = (atom1, atom2, atom3)

        # Try passing in a unitless distance
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site(
                atoms,
                distance_unitless,
                out_of_plane_angle,
                in_plane_angle,
                replace=True,
            )

        # Try passing in a unitless out_of_plane_angle
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site(
                atoms,
                distance,
                out_of_plane_angle_unitless,
                in_plane_angle,
                replace=True,
            )

        # Try passing in a unitless in_plane_angle
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site(
                atoms,
                distance,
                out_of_plane_angle,
                in_plane_angle_unitless,
                replace=True,
            )

        # Try giving two atoms
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_monovalent_lone_pair_virtual_site(
                [atom1, atom2],
                distance,
                out_of_plane_angle,
                in_plane_angle,
                replace=True,
            )

        # Successfully make a virtual site; throw exception is already present (replace=False)
        vsite1_index = molecule.add_monovalent_lone_pair_virtual_site(
            atoms, distance, out_of_plane_angle, in_plane_angle, replace=False
        )
        # TODO: Check if we get the same values back out from the @properties
        molecule_dict = molecule.to_dict()
        molecule2 = Molecule.from_dict(molecule_dict)
        assert molecule.to_dict() == molecule2.to_dict()

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_add_divalent_lone_pair_virtual_site(self, molecule):
        """Test addition of a DivalentLonePairVirtualSite to the Molecule"""
        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]
        distance = 0.3 * unit.angstrom
        out_of_plane_angle = 30 * unit.degree
        vsite1_index = molecule.add_divalent_lone_pair_virtual_site(
            [atom1, atom2, atom3],
            distance,
            out_of_plane_angle,
            replace=False,
        )

        # test giving too few atoms
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_divalent_lone_pair_virtual_site(
                [atom1, atom2],
                distance,
                out_of_plane_angle,
                replace=False,
            )
        molecule_dict = molecule.to_dict()
        molecule2 = Molecule.from_dict(molecule_dict)
        assert molecule_dict == molecule2.to_dict()

    @pytest.mark.parametrize("molecule", mini_drug_bank())
    def test_add_trivalent_lone_pair_virtual_site(self, molecule):
        """Test addition of a TrivalentLonePairVirtualSite to the Molecule"""
        # Do not modify the original molecule.
        molecule = copy.deepcopy(molecule)

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]
        distance = 0.3 * unit.angstrom

        vsite1_index = molecule.add_trivalent_lone_pair_virtual_site(
            [atom1, atom2, atom3, atom4],
            distance,
            replace=False,
        )
        # Test for assertion when giving too few atoms
        with pytest.raises(AssertionError) as excinfo:
            vsite1_index = molecule.add_trivalent_lone_pair_virtual_site(
                [atom1, atom2, atom3],
                distance,
                replace=True,
            )
        molecule_dict = molecule.to_dict()
        molecule2 = Molecule.from_dict(molecule_dict)
        assert molecule.to_dict() == molecule2.to_dict()

    def test_bond_charge_virtual_site_position_symmetric(self):
        """
        Test for expected positions of virtual sites when created through the molecule
        API. This uses symmetric, which will add both orientations of the atoms to
        create two particles for the vsite.

        The order of the positions is the same as the order of the particles as they
        appear in the molecule using Molecule.particles
        """

        import openff.toolkit.tests.test_forcefield

        molecule = openff.toolkit.tests.test_forcefield.create_water()
        # molecule = create_water()

        # a 90 degree water molecule on the xy plane
        molecule.add_conformer(
            unit.Quantity(
                np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                unit.angstrom,
            )
        )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        distance = 0.1 * unit.angstrom
        vsite1_index = molecule.add_bond_charge_virtual_site(
            [atom1, atom2],
            distance,
            replace=False,
            symmetric=True,
        )
        expected = np.array([[1.1, 0.0, 0.0], [-0.1, 0.0, 0.0]])

        vsite_pos = molecule.compute_virtual_site_positions_from_conformer(0)
        vsite_pos = vsite_pos / vsite_pos.unit

        assert vsite_pos.shape == (2, 3)
        assert np.allclose(vsite_pos, expected)

        conformer = molecule.conformers[0]
        vsite_pos = molecule.compute_virtual_site_positions_from_atom_positions(
            conformer
        )
        vsite_pos = vsite_pos / vsite_pos.unit

        assert vsite_pos.shape == (2, 3)
        assert np.allclose(vsite_pos, expected)

    def test_bond_charge_virtual_site_position(self):
        """
        Test for expected positions of virtual sites when created through the molecule
        API.

        The order of the positions is the same as the order of the particles as they
        appear in the molecule using Molecule.particles
        """

        import openff.toolkit.tests.test_forcefield

        molecule = openff.toolkit.tests.test_forcefield.create_water()
        # molecule = create_water()

        # a 90 degree water molecule on the xy plane
        molecule.add_conformer(
            unit.Quantity(
                np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                unit.angstrom,
            )
        )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        distance = 0.1 * unit.angstrom
        vsite1_index = molecule.add_bond_charge_virtual_site(
            [atom1, atom2],
            distance,
            replace=False,
            symmetric=False,
        )
        expected = np.array([[1.1, 0.0, 0.0]])

        vsite_pos = molecule.compute_virtual_site_positions_from_conformer(0)
        vsite_pos = vsite_pos / vsite_pos.unit

        assert vsite_pos.shape == (1, 3)
        assert np.allclose(vsite_pos, expected)

        conformer = molecule.conformers[0]
        vsite_pos = molecule.compute_virtual_site_positions_from_atom_positions(
            conformer
        )
        vsite_pos = vsite_pos / vsite_pos.unit

        assert vsite_pos.shape == (1, 3)
        assert np.allclose(vsite_pos, expected)

    def test_monovalent_lone_pair_virtual_site_position_symmetric(self):
        """
        Test for expected positions of virtual sites when created through the molecule
        API. This uses symmetric, which will add both orientations of the atoms to
        create two particles for the vsite.

        Note: symmetric with Monovalent is a bit counterintuitive. Since the particles
        origin is placed on the first atom, using a "symmetric" definition will flip
        the ordering, causing another particle to be placed on the third atom, reflected
        on the xy plane between the first particle.

        The order of the positions is the same as the order of the particles as they
        appear in the molecule using Molecule.particles
        """

        import openff.toolkit.tests.test_forcefield

        molecule = openff.toolkit.tests.test_forcefield.create_water()

        # a 90 degree water molecule on the xy plane
        molecule.add_conformer(
            unit.Quantity(
                np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                unit.angstrom,
            )
        )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        distance = 0.2 * unit.angstrom
        out_of_plane_angle = 90 * unit.degree
        in_plane_angle = 180 * unit.degree
        vsite1_index = molecule.add_monovalent_lone_pair_virtual_site(
            [atom1, atom2, atom3],
            distance,
            out_of_plane_angle=out_of_plane_angle,
            in_plane_angle=in_plane_angle,
            symmetric=True,
        )
        expected = np.array([[1.0, 0.0, -0.2], [0.0, 1.0, 0.2]])

        vsite_pos = molecule.compute_virtual_site_positions_from_conformer(0)
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (2, 3)
        assert np.allclose(vsite_pos, expected)

        conformer = molecule.conformers[0]
        vsite_pos = molecule.compute_virtual_site_positions_from_atom_positions(
            conformer
        )
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (2, 3)
        assert np.allclose(vsite_pos, expected)

    def test_monovalent_lone_pair_virtual_site_position(self):
        """
        Test for expected positions of virtual sites when created through the molecule
        API.

        The order of the positions is the same as the order of the particles as they
        appear in the molecule using Molecule.particles
        """

        import openff.toolkit.tests.test_forcefield

        molecule = openff.toolkit.tests.test_forcefield.create_water()

        # a 90 degree water molecule on the xy plane
        molecule.add_conformer(
            unit.Quantity(
                np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                unit.angstrom,
            )
        )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        distance = 0.2 * unit.angstrom
        out_of_plane_angle = 90 * unit.degree
        in_plane_angle = 180 * unit.degree
        vsite1_index = molecule.add_monovalent_lone_pair_virtual_site(
            [atom1, atom2, atom3],
            distance,
            out_of_plane_angle=out_of_plane_angle,
            in_plane_angle=in_plane_angle,
            symmetric=False,
        )
        expected = np.array([[1.0, 0.0, -0.2]])

        vsite_pos = molecule.compute_virtual_site_positions_from_conformer(0)
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (1, 3)
        assert np.allclose(vsite_pos, expected)

        conformer = molecule.conformers[0]
        vsite_pos = molecule.compute_virtual_site_positions_from_atom_positions(
            conformer
        )
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (1, 3)
        assert np.allclose(vsite_pos, expected)

    def test_divalent_lone_pair_virtual_site_position(self):
        """
        Test for expected positions of virtual sites when created through the molecule
        API.

        The order of the positions is the same as the order of the particles as they
        appear in the molecule using Molecule.particles
        """

        import openff.toolkit.tests.test_forcefield

        molecule = openff.toolkit.tests.test_forcefield.create_water()

        # a 90 degree water molecule on the xy plane
        molecule.add_conformer(
            unit.Quantity(
                np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                unit.angstrom,
            )
        )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        distance = 0.2 * unit.angstrom
        out_of_plane_angle = 90 * unit.degree
        vsite1_index = molecule.add_divalent_lone_pair_virtual_site(
            [atom1, atom2, atom3],
            distance,
            out_of_plane_angle=out_of_plane_angle,
            symmetric=False,
        )
        expected = np.array([[0.0, 0.0, -0.2]])

        vsite_pos = molecule.compute_virtual_site_positions_from_conformer(0)
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (1, 3)
        assert np.allclose(vsite_pos, expected)

        conformer = molecule.conformers[0]
        vsite_pos = molecule.compute_virtual_site_positions_from_atom_positions(
            conformer
        )
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (1, 3)
        assert np.allclose(vsite_pos, expected)

    def test_divalent_lone_pair_virtual_site_position_symmetric(self):
        """
        Test for expected positions of virtual sites when created through the molecule
        API. This uses symmetric, which will add both orientations of the atoms to
        create two particles for the vsite.

        Note: symmetric with Monovalent is a bit counterintuitive. Since the particles
        origin is placed on the first atom, using a "symmetric" definition will flip
        the ordering, causing another particle to be placed on the third atom, reflected
        on the xy plane between the first particle.

        The order of the positions is the same as the order of the particles as they
        appear in the molecule using Molecule.particles
        """

        import openff.toolkit.tests.test_forcefield

        molecule = openff.toolkit.tests.test_forcefield.create_water()

        # a 90 degree water molecule on the xy plane
        molecule.add_conformer(
            unit.Quantity(
                np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                unit.angstrom,
            )
        )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        distance = 0.2 * unit.angstrom
        out_of_plane_angle = 90 * unit.degree
        vsite1_index = molecule.add_divalent_lone_pair_virtual_site(
            [atom1, atom2, atom3],
            distance,
            out_of_plane_angle=out_of_plane_angle,
            symmetric=True,
        )
        expected = np.array([[0.0, 0.0, -0.2], [0.0, 0.0, 0.2]])

        vsite_pos = molecule.compute_virtual_site_positions_from_conformer(0)
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (2, 3)
        assert np.allclose(vsite_pos, expected)

        conformer = molecule.conformers[0]
        vsite_pos = molecule.compute_virtual_site_positions_from_atom_positions(
            conformer
        )
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (2, 3)
        assert np.allclose(vsite_pos, expected)

    def test_trivalent_lone_pair_virtual_site_position(self):
        """
        Test for expected positions of virtual sites when created through the molecule
        API.

        The order of the positions is the same as the order of the particles as they
        appear in the molecule using Molecule.particles
        """

        import openff.toolkit.tests.test_forcefield

        molecule = openff.toolkit.tests.test_forcefield.create_ammonia()

        # an ammonia; the central nitrogen .5 above the hydrogen, which all lie one
        # an xy plane
        molecule.add_conformer(
            unit.Quantity(
                np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.5],
                        [-0.5, 0.5, 0.0],
                        [-0.5, -0.5, 0.0],
                    ]
                ),
                unit.angstrom,
            )
        )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]
        atom3 = molecule.atoms[2]
        atom4 = molecule.atoms[3]
        distance = 1.0 * unit.angstrom
        vsite1_index = molecule.add_trivalent_lone_pair_virtual_site(
            [atom1, atom2, atom3, atom4], distance, symmetric=False
        )
        expected = np.array([[0.0, 0.0, 1.5]])

        vsite_pos = molecule.compute_virtual_site_positions_from_conformer(0)
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (1, 3)
        assert np.allclose(vsite_pos, expected)

        conformer = molecule.conformers[0]
        vsite_pos = molecule.compute_virtual_site_positions_from_atom_positions(
            conformer
        )
        vsite_pos = vsite_pos / vsite_pos.unit
        assert vsite_pos.shape == (1, 3)
        assert np.allclose(vsite_pos, expected)

    @requires_openeye
    def test_chemical_environment_matches_OE(self):
        """Test chemical environment matches"""
        # TODO: Move this to test_toolkits, test all available toolkits
        # Create chiral molecule
        from simtk.openmm.app import element

        toolkit_wrapper = OpenEyeToolkitWrapper()
        molecule = Molecule()
        atom_C = molecule.add_atom(
            element.carbon.atomic_number, 0, False, stereochemistry="R", name="C"
        )
        atom_H = molecule.add_atom(element.hydrogen.atomic_number, 0, False, name="H")
        atom_Cl = molecule.add_atom(element.chlorine.atomic_number, 0, False, name="Cl")
        atom_Br = molecule.add_atom(element.bromine.atomic_number, 0, False, name="Br")
        atom_F = molecule.add_atom(element.fluorine.atomic_number, 0, False, name="F")
        molecule.add_bond(atom_C, atom_H, 1, False)
        molecule.add_bond(atom_C, atom_Cl, 1, False)
        molecule.add_bond(atom_C, atom_Br, 1, False)
        molecule.add_bond(atom_C, atom_F, 1, False)
        # Test known cases
        matches = molecule.chemical_environment_matches(
            "[#6:1]", toolkit_registry=toolkit_wrapper
        )
        assert (
            len(matches) == 1
        )  # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 1  # it should have one tagged atom
        assert set(matches[0]) == set([atom_C])
        matches = molecule.chemical_environment_matches(
            "[#6:1]~[#1:2]", toolkit_registry=toolkit_wrapper
        )
        assert (
            len(matches) == 1
        )  # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 2  # it should have two tagged atoms
        assert set(matches[0]) == set([atom_C, atom_H])
        matches = molecule.chemical_environment_matches(
            "[Cl:1]-[C:2]-[H:3]", toolkit_registry=toolkit_wrapper
        )
        assert (
            len(matches) == 1
        )  # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 3  # it should have three tagged atoms
        assert set(matches[0]) == set([atom_Cl, atom_C, atom_H])
        matches = molecule.chemical_environment_matches(
            "[#6:1]~[*:2]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 4  # there should be four matches
        for match in matches:
            assert len(match) == 2  # each match should have two tagged atoms
        # Test searching for stereo-specific SMARTS
        matches = molecule.chemical_environment_matches(
            "[#6@:1](-[F:2])(-[Cl])(-[Br])(-[H])", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 1  # there should be one match
        for match in matches:
            assert len(match) == 2  # each match should have two tagged atoms
        matches = molecule.chemical_environment_matches(
            "[#6@@:1](-[F:2])(-[Cl])(-[Br])(-[H])", toolkit_registry=toolkit_wrapper
        )
        assert (
            len(matches) == 0
        )  # this is the wrong stereochemistry, so there shouldn't be any matches

    # TODO: Test forgive undef amide enol stereo
    # TODO: test forgive undef phospho linker stereo
    # TODO: test forgive undef C=NH stereo
    # TODO: test forgive undef phospho stereo
    # Potentially better OE stereo check: OEFlipper — Toolkits - - Python
    # https: // docs.eyesopen.com / toolkits / python / omegatk / OEConfGenFunctions / OEFlipper.html

    @requires_rdkit
    def test_chemical_environment_matches_RDKit(self):
        """Test chemical environment matches"""
        # Create chiral molecule
        from simtk.openmm.app import element

        toolkit_wrapper = RDKitToolkitWrapper()
        molecule = Molecule()
        atom_C = molecule.add_atom(
            element.carbon.atomic_number, 0, False, stereochemistry="R", name="C"
        )
        atom_H = molecule.add_atom(element.hydrogen.atomic_number, 0, False, name="H")
        atom_Cl = molecule.add_atom(element.chlorine.atomic_number, 0, False, name="Cl")
        atom_Br = molecule.add_atom(element.bromine.atomic_number, 0, False, name="Br")
        atom_F = molecule.add_atom(element.fluorine.atomic_number, 0, False, name="F")
        molecule.add_bond(atom_C, atom_H, 1, False)
        molecule.add_bond(atom_C, atom_Cl, 1, False)
        molecule.add_bond(atom_C, atom_Br, 1, False)
        molecule.add_bond(atom_C, atom_F, 1, False)
        # Test known cases
        matches = molecule.chemical_environment_matches(
            "[#6:1]", toolkit_registry=toolkit_wrapper
        )
        assert (
            len(matches) == 1
        )  # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 1  # it should have one tagged atom
        assert set(matches[0]) == set([atom_C])
        matches = molecule.chemical_environment_matches(
            "[#6:1]~[#1:2]", toolkit_registry=toolkit_wrapper
        )
        assert (
            len(matches) == 1
        )  # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 2  # it should have two tagged atoms
        assert set(matches[0]) == set([atom_C, atom_H])
        matches = molecule.chemical_environment_matches(
            "[Cl:1]-[C:2]-[H:3]", toolkit_registry=toolkit_wrapper
        )
        assert (
            len(matches) == 1
        )  # there should be a unique match, so one atom tuple is returned
        assert len(matches[0]) == 3  # it should have three tagged atoms
        assert set(matches[0]) == set([atom_Cl, atom_C, atom_H])
        matches = molecule.chemical_environment_matches(
            "[#6:1]~[*:2]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 4  # there should be four matches
        for match in matches:
            assert len(match) == 2  # each match should have two tagged atoms
        # Test searching for stereo-specific SMARTS
        matches = molecule.chemical_environment_matches(
            "[#6@:1](-[F:2])(-[Cl])(-[Br])(-[H])", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 1  # there should be one match
        for match in matches:
            assert len(match) == 2  # each match should have two tagged atoms
        matches = molecule.chemical_environment_matches(
            "[#6@@:1](-[F:2])(-[Cl])(-[Br])(-[H])", toolkit_registry=toolkit_wrapper
        )
        assert (
            len(matches) == 0
        )  # this is the wrong stereochemistry, so there shouldn't be any matches

    @requires_rdkit
    @requires_openeye
    @pytest.mark.slow
    @pytest.mark.parametrize(
        ("toolkit", "method"),
        [
            ("openeye", "mmff94"),
            ("openeye", "am1bcc"),
            ("openeye", "am1-mulliken"),
            ("openeye", "gasteiger"),
            ("openeye", "am1bccnosymspt"),
            ("openeye", "am1elf10"),
            ("openeye", "am1bccelf10"),
            ("ambertools", "am1bcc"),
            ("ambertools", "gasteiger"),
            ("ambertools", "am1-mulliken"),
        ],
    )
    def test_assign_partial_charges(self, toolkit, method):
        """Test computation/retrieval of partial charges"""
        # TODO: Test only one molecule for speed?
        # TODO: Do we need to deepcopy each molecule, or is setUp called separately for each test method?
        import numpy as np
        from simtk import unit

        # Do not modify original molecules.
        # molecules = copy.deepcopy(mini_drug_bank())
        # In principle, testing for charge assignment over a wide set of molecules is important, but
        # I think that's covered in test_toolkits.py. Here, we should only be concerned with testing the
        # Molecule API, and therefore a single molecule should be sufficient
        molecule = Molecule.from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

        if toolkit == "openeye":
            toolkit_registry = ToolkitRegistry(
                toolkit_precedence=[OpenEyeToolkitWrapper]
            )
        elif toolkit == "ambertools":
            toolkit_registry = ToolkitRegistry(
                toolkit_precedence=[AmberToolsToolkitWrapper]
            )

        molecule.assign_partial_charges(
            partial_charge_method=method,
            toolkit_registry=toolkit_registry,
        )
        initial_charges = molecule._partial_charges

        # Make sure everything isn't 0s
        assert (abs(initial_charges / unit.elementary_charge) > 0.01).any()
        # Check total charge
        charges_sum_unitless = initial_charges.sum() / unit.elementary_charge
        total_charge_unitless = molecule.total_charge / unit.elementary_charge
        # if abs(charges_sum_unitless - total_charge_unitless) > 0.0001:
        # print(
        #     "molecule {}    charge_sum {}     molecule.total_charge {}".format(
        #         molecule.name, charges_sum_unitless, total_charge_unitless
        #     )
        # )
        np.allclose(charges_sum_unitless, total_charge_unitless, atol=0.002)

        # Call should be faster second time due to caching
        # TODO: Implement caching
        molecule.assign_partial_charges(
            partial_charge_method=method, toolkit_registry=toolkit_registry
        )
        recomputed_charges = molecule._partial_charges
        assert np.allclose(initial_charges, recomputed_charges, atol=0.002)

    @pytest.mark.parametrize("toolkit", ["openeye", "rdkit"])
    def test_apply_elf_conformer_selection(self, toolkit):
        """Test applying the ELF10 method."""

        if toolkit == "openeye":
            pytest.importorskip("openeye")
            toolkit_registry = ToolkitRegistry(
                toolkit_precedence=[OpenEyeToolkitWrapper]
            )
        elif toolkit == "rdkit":
            pytest.importorskip("rdkit")
            toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])

        molecule = Molecule.from_file(
            get_data_file_path(os.path.join("molecules", "z_3_hydroxy_propenal.sdf")),
            "SDF",
        )

        initial_conformers = [
            # Add a conformer with an internal H-bond.
            np.array(
                [
                    [0.5477, 0.3297, -0.0621],
                    [-0.1168, -0.7881, 0.2329],
                    [-1.4803, -0.8771, 0.1667],
                    [-0.2158, 1.5206, -0.4772],
                    [-1.4382, 1.5111, -0.5580],
                    [1.6274, 0.3962, -0.0089],
                    [0.3388, -1.7170, 0.5467],
                    [-1.8612, -0.0347, -0.1160],
                    [0.3747, 2.4222, -0.7115],
                ]
            )
            * unit.angstrom,
            # Add a conformer without an internal H-bond.
            np.array(
                [
                    [0.5477, 0.3297, -0.0621],
                    [-0.1168, -0.7881, 0.2329],
                    [-1.4803, -0.8771, 0.1667],
                    [-0.2158, 1.5206, -0.4772],
                    [0.3353, 2.5772, -0.7614],
                    [1.6274, 0.3962, -0.0089],
                    [0.3388, -1.7170, 0.5467],
                    [-1.7743, -1.7634, 0.4166],
                    [-1.3122, 1.4082, -0.5180],
                ]
            )
            * unit.angstrom,
        ]

        molecule._conformers = [*initial_conformers]

        # Apply ELF10
        molecule.apply_elf_conformer_selection(toolkit_registry=toolkit_registry)
        elf10_conformers = molecule.conformers

        assert len(elf10_conformers) == 1

        assert np.allclose(
            elf10_conformers[0].value_in_unit(unit.angstrom),
            initial_conformers[1].value_in_unit(unit.angstrom),
        )

    @requires_openeye
    def test_assign_fractional_bond_orders(self):
        """Test assignment of fractional bond orders"""
        # TODO: Test only one molecule for speed?
        # TODO: Do we need to deepcopy each molecule, or is setUp called separately for each test method?

        # Do not modify the original molecules.
        molecules = copy.deepcopy(mini_drug_bank())

        toolkits_to_bondorder_method = {
            (OpenEyeToolkitWrapper,): ["am1-wiberg", "pm3-wiberg"]
        }
        # Don't test AmberTools here since it takes too long
        # (AmberToolsToolkitWrapper, RDKitToolkitWrapper):['am1-wiberg']}
        for toolkits in list(toolkits_to_bondorder_method.keys()):
            toolkit_registry = ToolkitRegistry(toolkit_precedence=toolkits)
            for bond_order_model in toolkits_to_bondorder_method[toolkits]:
                for molecule in molecules[
                    :5
                ]:  # Just test first five molecules for speed
                    molecule.generate_conformers(toolkit_registry=toolkit_registry)
                    molecule.assign_fractional_bond_orders(
                        bond_order_model=bond_order_model,
                        toolkit_registry=toolkit_registry,
                        use_conformers=molecule.conformers,
                    )
                    fbo1 = [bond.fractional_bond_order for bond in molecule.bonds]
                    # TODO: Now that the assign_fractional_bond_orders function takes more kwargs,
                    #       how can we meaningfully cache its results?
                    # # Call should be faster the second time due to caching
                    # molecule.assign_fractional_bond_orders(bond_order_model=bond_order_model,
                    #                                        toolkit_registry=toolkit_registry)
                    # fbo2 = [bond.fractional_bond_order for bond in molecule.bonds]
                    # np.testing.assert_allclose(fbo1, fbo2, atol=1.e-4)

    def test_get_bond_between(self):
        """Test Molecule.get_bond_between"""
        mol = Molecule.from_smiles("C#C")

        # Dynamically get atoms in case from_smiles produces different atom order
        hydrogens = [a for a in mol.atoms if a.atomic_number == 1]
        carbons = [a for a in mol.atoms if a.atomic_number == 6]

        bond_from_atoms = mol.get_bond_between(carbons[0], carbons[1])
        bond_from_idx = mol.get_bond_between(
            carbons[0].molecule_atom_index, carbons[1].molecule_atom_index
        )
        assert bond_from_idx == bond_from_atoms

        with pytest.raises(NotBondedError):
            mol.get_bond_between(hydrogens[0], hydrogens[1])

    @pytest.mark.parametrize(
        ("smiles", "n_rings"),
        [
            ("CCO", 0),
            ("c1ccccc1", 1),
            ("C1CC2CCC1C2", 2),  # This should probably be 3?
            ("c1ccc(cc1)c2ccccc2", 2),
            ("c1ccc2ccccc2c1", 2),
            ("c1ccc2cc3ccccc3cc2c1", 3),
            ("C1C2C(CCC1)CC5C3C2CCC7C3C4C(CCC6C4C5CCC6)CC7", 7),
        ],
    )
    @requires_rdkit
    def test_molecule_rings(self, smiles, n_rings):
        """Test the Molecule.rings property"""
        assert (
            n_rings == Molecule.from_smiles(smiles, allow_undefined_stereo=True).n_rings
        )

    @pytest.mark.parametrize(
        ("smiles", "n_atom_rings", "n_bond_rings"),
        [
            ("c1ccc2ccccc2c1", 10, 11),
            ("c1ccc(cc1)c2ccccc2", 12, 12),
            ("Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C", 30, 30),
        ],
    )
    @requires_rdkit
    def test_is_in_ring(self, smiles, n_atom_rings, n_bond_rings):
        """Test Atom.is_in_ring and Bond.is_in_ring"""
        mol = Molecule.from_smiles(smiles)

        assert len([atom for atom in mol.atoms if atom.is_in_ring]) == n_atom_rings
        assert len([bond for bond in mol.bonds if bond.is_in_ring]) == n_bond_rings

    @requires_pkg("ipython")
    @requires_rdkit
    def test_visualize_rdkit(self):
        """Test that the visualize method returns an expected object when using RDKit to generate a 2-D representation"""
        import rdkit

        mol = Molecule().from_smiles("CCO")

        assert isinstance(mol.visualize(backend="rdkit"), rdkit.Chem.rdchem.Mol)

    @pytest.mark.skipif(
        has_pkg("rdkit"),
        reason="Test requires that RDKit is not installed",
    )
    def test_visualize_fallback(self):
        """Test falling back from RDKit to OpenEye if RDKit is specified but not installed"""
        mol = Molecule().from_smiles("CCO")
        with pytest.warns(UserWarning):
            mol.visualize(backend="rdkit")

    @requires_pkg("nglview")
    def test_visualize_nglview(self):
        """Test that the visualize method returns an NGLview widget. Note that
        nglview is not explicitly a requirement in the test environment, but
        is likely to get pulled in with other dependencies."""
        try:
            import nglview
        except ModuleNotFoundError:
            pass

        # Start with a molecule without conformers
        mol = Molecule().from_smiles("CCO")

        with pytest.raises(ValueError):
            mol.visualize(backend="nglview")

        # Add conformers
        mol.generate_conformers()

        # Ensure an NGLView widget is returned
        assert isinstance(mol.visualize(backend="nglview"), nglview.NGLWidget)

    @requires_pkg("ipython")
    @requires_openeye
    def test_visualize_openeye(self):
        """Test that the visualize method returns an expected object when using OpenEye to generate a 2-D representation"""
        import IPython

        mol = Molecule().from_smiles("CCO")

        assert isinstance(mol.visualize(backend="openeye"), IPython.core.display.Image)


class MyMol(FrozenMolecule):
    """
    Lightweight FrozenMolecule subclass for molecule-subclass tests below
    """


class TestMoleculeSubclass:
    """
    Test that the FrozenMolecule class is subclass-able, by ensuring that Molecule.from_X constructors
    return objects of the correct types
    """

    def test_molecule_subclass_from_smiles(self):
        """Ensure that the right type of object is returned when running MyMol.from_smiles"""
        mol = MyMol.from_smiles("CCO")
        assert isinstance(mol, MyMol)

    def test_molecule_subclass_from_inchi(self):
        """Ensure that the right type of object is returned when running MyMol.from_inchi"""
        mol = MyMol.from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
        assert isinstance(mol, MyMol)

    @requires_openeye
    def test_molecule_subclass_from_iupac(self):
        """Ensure that the right type of object is returned when running MyMol.from_iupac"""
        mol = MyMol.from_iupac("benzene")
        assert isinstance(mol, MyMol)

    def test_molecule_subclass_from_file(self):
        """Ensure that the right type of object is returned when running MyMol.from_file"""
        mol = MyMol.from_file(get_data_file_path("molecules/ethanol.sdf"))
        assert isinstance(mol, MyMol)

    def test_molecule_subclass_from_mapped_smiles(self):
        """Ensure that the right type of object is returned when running MyMol.from_mapped_smiles"""
        mol = MyMol.from_mapped_smiles("[H:1][C:2]([H:3])([H:4])([H:5])")
        assert isinstance(mol, MyMol)

    @requires_pkg("qcelemental")
    @requires_pkg("qcportal")
    def test_molecule_subclass_from_qcschema(self):
        """Ensure that the right type of object is returned when running MyMol.from_qcschema"""
        import qcportal as ptl

        client = ptl.FractalClient()
        ds = client.get_collection(
            "TorsionDriveDataset", "Fragment Stability Benchmark"
        )
        entry = ds.get_entry(
            "CC(=O)Nc1cc2c(cc1OC)nc[n:4][c:3]2[NH:2][c:1]3ccc(c(c3)Cl)F"
        )
        # now make the molecule from the record instance with and without the geometry
        mol = MyMol.from_qcschema(entry.dict(encoding="json"))
        assert isinstance(mol, MyMol)
        # Make from object, which will include geometry
        mol = MyMol.from_qcschema(entry, client)
        assert isinstance(mol, MyMol)

    def test_molecule_subclass_from_topology(self):
        """Ensure that the right type of object is returned when running MyMol.from_topology"""
        top = Molecule.from_smiles("CCO").to_topology()
        mol = MyMol.from_topology(top)
        assert isinstance(mol, MyMol)

    @requires_rdkit
    def test_molecule_subclass_from_pdb_and_smiles(self):
        """Ensure that the right type of object is returned when running MyMol.from_pdb_and_smiles"""
        mol = MyMol.from_pdb_and_smiles(
            get_data_file_path("molecules/toluene.pdb"), "Cc1ccccc1"
        )
        assert isinstance(mol, MyMol)

    def test_molecule_copy_constructor_from_other_subclass(self):
        """Ensure that the right type of object is returned when running the MyMol copy constructor"""
        normal_mol = MyMol.from_smiles("CCO")
        mol = MyMol(normal_mol)
        assert isinstance(mol, MyMol)

    def test_molecule_subclass_from_dict(self):
        """Ensure that the right type of object is returned when running MyMol.from_dict"""
        orig_mol = Molecule.from_smiles("CCO")
        mol = MyMol.from_dict(orig_mol.to_dict())
        assert isinstance(mol, MyMol)
