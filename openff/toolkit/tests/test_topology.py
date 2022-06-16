"""
Tests for Topology

"""

import itertools

import numpy as np
import pytest
from openff.units import unit
from openmm import app

from openff.toolkit.tests.create_molecules import (
    create_cyclohexane,
    create_ethanol,
    create_reversed_ethanol,
    cyx_hierarchy_added,
    dipeptide,
    dipeptide_hierarchy_added,
    dipeptide_residues_perceived,
    ethane_from_smiles,
    ethene_from_smiles,
    propane_from_smiles,
    toluene_from_sdf,
    topology_with_metadata,
)
from openff.toolkit.tests.utils import (
    get_data_file_path,
    requires_openeye,
    requires_pkg,
    requires_rdkit,
)
from openff.toolkit.topology import (
    Atom,
    ImproperDict,
    Molecule,
    TagSortedDict,
    Topology,
    ValenceDict,
)
from openff.toolkit.utils import (
    BASIC_CHEMINFORMATICS_TOOLKITS,
    OPENEYE_AVAILABLE,
    RDKIT_AVAILABLE,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
)
from openff.toolkit.utils.exceptions import (
    AtomNotInTopologyError,
    DuplicateUniqueMoleculeError,
    InvalidBoxVectorsError,
    InvalidPeriodicityError,
    MissingUniqueMoleculesError,
    MoleculeNotInTopologyError,
)


def assert_tuple_of_atoms_equal(
    atom_tuples1, atom_tuples2, transformed_dict_cls=ValenceDict
):
    """Check that two lists of atoms are the same.

    The function compares that the parent molecules are isomorphic and
    that the molecule index is the same.
    """
    assert len(atom_tuples1) == len(atom_tuples2)

    # They are atoms of isomorphic molecules. We assume here that all
    # atoms in the same list of tuples belong to the same molecule so
    # that we can perform the check only once.
    molecule1 = atom_tuples1[0][0]._molecule
    molecule2 = atom_tuples2[0][0]._molecule
    assert molecule1 == molecule2
    for atom_tuple in atom_tuples1:
        for a in atom_tuple:
            assert a._molecule is molecule1
    for atom_tuple in atom_tuples2:
        for a in atom_tuple:
            assert a._molecule is molecule2

    # All atoms are equal. Use ValenceDict for this
    atom_indices = []
    for atom_tuples in [atom_tuples1, atom_tuples2]:
        valence_dict = transformed_dict_cls()
        for atom_tuple in atom_tuples:
            key = tuple(a.molecule_atom_index for a in atom_tuple)
            valence_dict[key] = atom_tuple
        atom_indices.append(valence_dict)
    assert set(atom_indices[0]) == set(atom_indices[1])


# IF we've done our jobs right, it shouldn't matter which toolkit the tests for Topology run using (both's behaviors
# should be indistinguishable)
def test_cheminformatics_toolkit_is_installed():
    """Ensure that at least one supported cheminformatics toolkit is installed."""
    if not (RDKIT_AVAILABLE) and not (OPENEYE_AVAILABLE):
        msg = "No supported cheminformatics toolkits are installed. Please install a supported toolkit:\n"
        msg += str(BASIC_CHEMINFORMATICS_TOOLKITS)
        raise Exception(msg)


# TODO: Refactor this to pytest
class TestTopology:
    def test_empty(self):
        """Test creation of empty topology"""
        topology = Topology()
        assert topology.n_molecules == 0
        assert topology.n_atoms == 0
        assert topology.n_bonds == 0
        assert topology.box_vectors is None
        assert not topology.is_periodic
        assert len(topology.constrained_atom_pairs.items()) == 0

    def test_deprecated_api_points(self):
        """Ensure that some of the API deprecated circa v0.11.0 still exist."""
        from openff.toolkit.topology.topology import TopologyDeprecationWarning

        topology = Topology()

        for key in ["molecules", "atoms", "bonds", "particles"]:
            old_iterator = "topology_" + key
            old_counter = "n_topology_" + key
            with pytest.warns(
                TopologyDeprecationWarning,
                match=f"Topology.{old_counter} is deprecated. Use Topology.n_{key} instead.",
            ):
                assert getattr(topology, old_counter) == 0
            with pytest.warns(
                TopologyDeprecationWarning,
                match=f"Topology.{old_iterator} is deprecated. Use Topology.{key} instead.",
            ):
                assert len([*getattr(topology, old_iterator)]) == 0

        # Some particle-related methods are deprecated for reasons other than the `TopologyX` deprecation
        with pytest.warns(
            TopologyDeprecationWarning,
            match="Topology.n_particles is deprecated. Use Topology.n_atoms instead.",
        ):
            assert topology.n_particles == 0

        with pytest.warns(
            TopologyDeprecationWarning,
            match="Topology.particles is deprecated. Use Topology.atoms instead.",
        ):
            assert len([*topology.particles]) == 0

        topology = Molecule.from_smiles("O").to_topology()
        first_atom = [*topology.atoms][0]

        with pytest.warns(
            TopologyDeprecationWarning,
            match="Topology.particle_index is deprecated. Use Topology.atom_index instead.",
        ):
            assert topology.particle_index(first_atom) == 0

    def test_reinitialization_box_vectors(self):
        topology = Topology()
        assert Topology(topology).box_vectors is None

        topology.box_vectors = [1, 2, 3] * unit.nanometer
        topology_copy = Topology(topology)

        assert (topology.box_vectors == topology_copy.box_vectors).all()

    def test_box_vectors(self):
        """Test the getter and setter for box_vectors"""
        topology = Topology()
        good_box_vectors = unit.Quantity(np.eye(3) * 20, unit.angstrom)
        one_dim_vectors = unit.Quantity(np.ones(3) * 20, unit.angstrom)
        list_vectors = unit.Quantity([20, 20, 20], unit.angstrom)
        list_list_vectors = unit.Quantity(
            [[20, 0, 0], [0, 20, 0], [0, 0, 20]], unit.angstrom
        )
        bad_shape_vectors = unit.Quantity(np.ones(2) * 20, unit.angstrom)
        bad_units_vectors = unit.Quantity(np.ones(3) * 20, unit.year)
        bad_type_vectors = unit.Quantity(1.0, unit.nanometer)
        unitless_vectors = np.array([10, 20, 30])

        assert topology.box_vectors is None

        for bad_vectors in [
            bad_shape_vectors,
            bad_units_vectors,
            bad_type_vectors,
            unitless_vectors,
        ]:
            with pytest.raises(InvalidBoxVectorsError):
                topology.box_vectors = bad_vectors
            assert topology.box_vectors is None

        for good_vectors in [
            good_box_vectors,
            one_dim_vectors,
            list_vectors,
            list_list_vectors,
        ]:
            topology.box_vectors = good_vectors
            assert (topology.box_vectors == good_vectors * np.eye(3)).all()

    def test_is_periodic(self):
        """Test the getter and setter for is_periodic"""
        vacuum_top = Topology()
        assert vacuum_top.is_periodic is False

        with pytest.raises(InvalidPeriodicityError):
            vacuum_top.is_periodic = True

        solvent_box = Topology()
        solvent_box.box_vectors = np.eye(3) * 4 * unit.nanometer
        assert solvent_box.is_periodic is True

        with pytest.raises(InvalidPeriodicityError):
            solvent_box.is_periodic = False
        solvent_box.box_vectors = None
        assert solvent_box.is_periodic is False

    def test_from_smiles(self):
        """Test creation of a OpenFF Topology object from a SMILES string"""
        topology = Topology.from_molecules(ethane_from_smiles())

        assert topology.n_molecules == 1
        assert topology.n_atoms == 8
        assert topology.n_bonds == 7
        assert topology.box_vectors is None
        assert len(topology.constrained_atom_pairs.items()) == 0

        topology.add_molecule(ethane_from_smiles())
        assert topology.n_molecules == 2
        assert topology.n_atoms == 16
        assert topology.n_bonds == 14
        assert topology.box_vectors is None
        assert len(topology.constrained_atom_pairs.items()) == 0

    def test_from_smiles_unique_mols(self):
        """Test the addition of two different molecules to a topology"""
        topology = Topology.from_molecules(
            [ethane_from_smiles(), propane_from_smiles()]
        )
        assert topology.n_molecules == 2

    def test_n_atoms(self):
        """Test n_atoms function"""
        topology = Topology()
        assert topology.n_atoms == 0
        assert topology.n_bonds == 0
        topology.add_molecule(ethane_from_smiles())
        assert topology.n_atoms == 8
        assert topology.n_bonds == 7

    def test_get_atom(self):
        """Test Topology.atom function (atom lookup from index)"""
        topology = Topology()
        topology.add_molecule(ethane_from_smiles())
        with pytest.raises(Exception):
            topology.atom(-1)

        # Make sure we get 2 carbons and 8 hydrogens
        n_carbons = 0
        n_hydrogens = 0
        for index in range(8):
            if topology.atom(index).atomic_number == 6:
                n_carbons += 1
            if topology.atom(index).atomic_number == 1:
                n_hydrogens += 1
        assert n_carbons == 2
        assert n_hydrogens == 6

        with pytest.raises(Exception):
            topology.atom(8)

    def test_atom_index(self):
        topology = create_ethanol().to_topology()

        for index in range(topology.n_atoms):
            atom = topology.atom(index)
            assert topology.atom_index(atom) == index

        ghost_atom = Atom(atomic_number=1, formal_charge=0, is_aromatic=False)

        with pytest.raises(AtomNotInTopologyError):
            topology.atom_index(ghost_atom)

    def test_molecule_index(self):
        molecules = [Molecule.from_smiles("CCO"), Molecule.from_smiles("O")]

        topology = Topology.from_molecules(molecules)

        for index in range(topology.n_molecules):
            molecule = topology.molecule(index)
            assert topology.molecule_index(molecule) == index

        ghost_molecule = Molecule.from_smiles("N")

        with pytest.raises(MoleculeNotInTopologyError):
            topology.molecule_index(ghost_molecule)

    def test_atom_element_properties(self):
        """
        Test element-like getters of TopologyAtom atomic number. In 0.11.0, Atom.element
        was removed and replaced with Atom.atomic_number and Atom.symbol.
        """
        topology = Topology()
        topology.add_molecule(toluene_from_sdf())

        first_atom = topology.atom(0)
        eighth_atom = topology.atom(7)

        # These atoms are expected to be hydrogen and carbon, respectively
        assert first_atom.symbol == "C"
        assert first_atom.atomic_number == 6
        assert eighth_atom.symbol == "H"
        assert eighth_atom.atomic_number == 1

    def test_get_bond(self):
        """Test Topology.bond function (bond lookup from index)"""
        topology = Topology()
        topology.add_molecule(ethane_from_smiles())
        topology.add_molecule(ethene_from_smiles())
        with pytest.raises(Exception):
            topology.bond(-1)

        n_single_bonds = 0
        n_double_bonds = 0
        n_ch_bonds = 0
        n_cc_bonds = 0
        for index in range(12):  # 7 from ethane, 5 from ethene
            topology_bond = topology.bond(index)
            if topology_bond.bond_order == 1:
                n_single_bonds += 1
            if topology_bond.bond_order == 2:
                n_double_bonds += 1
            n_bond_carbons = 0
            n_bond_hydrogens = 0
            for atom in topology_bond.atoms:
                if atom.atomic_number == 6:
                    n_bond_carbons += 1
                if atom.atomic_number == 1:
                    n_bond_hydrogens += 1
            if n_bond_carbons == 2:
                n_cc_bonds += 1
            if n_bond_carbons == 1 and n_bond_hydrogens == 1:
                n_ch_bonds += 1

        assert n_single_bonds == 11
        assert n_double_bonds == 1
        assert n_cc_bonds == 2
        assert n_ch_bonds == 10

        with pytest.raises(Exception):
            topology_bond = topology.bond(12)

    def test_angles(self):
        """Topology.angles should return image angles of all topology molecules."""
        molecule1 = ethane_from_smiles()
        molecule2 = propane_from_smiles()

        # Create topology.
        topology = Topology()
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule2)

        # The topology should have the correct number of angles.
        topology_angles = list(topology.angles)
        assert len(topology_angles) == topology.n_angles
        assert topology.n_angles == 2 * molecule1.n_angles + molecule2.n_angles

        # Check that the topology angles are the correct ones.
        mol_angle_atoms1 = list(molecule1.angles)
        mol_angle_atoms2 = list(molecule2.angles)
        top_angle_atoms1 = [
            tuple(a for a in atoms) for atoms in topology_angles[: molecule1.n_angles]
        ]
        top_angle_atoms2 = [
            tuple(a for a in atoms)
            for atoms in topology_angles[molecule1.n_angles : 2 * molecule1.n_angles]
        ]
        top_angle_atoms3 = [
            tuple(a for a in atoms)
            for atoms in topology_angles[2 * molecule1.n_angles :]
        ]

        assert_tuple_of_atoms_equal(top_angle_atoms1, mol_angle_atoms1)
        assert_tuple_of_atoms_equal(top_angle_atoms2, mol_angle_atoms1)
        assert_tuple_of_atoms_equal(top_angle_atoms3, mol_angle_atoms2)

    def test_propers(self):
        """Topology.propers should return image propers torsions of all topology molecules."""
        molecule1 = ethane_from_smiles()
        molecule2 = propane_from_smiles()

        # Create topology.
        topology = Topology()
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule2)

        # The topology should have the correct number of propers.
        topology_propers = list(topology.propers)
        assert len(topology_propers) == topology.n_propers
        assert topology.n_propers == 2 * molecule1.n_propers + molecule2.n_propers

        # Check that the topology propers are the correct ones.
        mol_proper_atoms1 = list(molecule1.propers)
        mol_proper_atoms2 = list(molecule2.propers)
        top_proper_atoms1 = [
            tuple(a for a in atoms) for atoms in topology_propers[: molecule1.n_propers]
        ]
        top_proper_atoms2 = [
            tuple(a for a in atoms)
            for atoms in topology_propers[molecule1.n_propers : 2 * molecule1.n_propers]
        ]
        top_proper_atoms3 = [
            tuple(a for a in atoms)
            for atoms in topology_propers[2 * molecule1.n_propers :]
        ]

        assert_tuple_of_atoms_equal(top_proper_atoms1, mol_proper_atoms1)
        assert_tuple_of_atoms_equal(top_proper_atoms2, mol_proper_atoms1)
        assert_tuple_of_atoms_equal(top_proper_atoms3, mol_proper_atoms2)

    def test_impropers(self):
        """Topology.impropers should return image impropers torsions of all topology molecules."""
        molecule1 = ethane_from_smiles()
        molecule2 = propane_from_smiles()

        # Create topology.
        topology = Topology()
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule1)
        topology.add_molecule(molecule2)

        # The topology should have the correct number of impropers.
        topology_impropers = list(topology.impropers)
        assert len(topology_impropers) == topology.n_impropers
        assert topology.n_impropers == 2 * molecule1.n_impropers + molecule2.n_impropers

        # Check that the topology impropers are the correct ones.
        mol_improper_atoms1 = list(molecule1.impropers)
        mol_improper_atoms2 = list(molecule2.impropers)
        top_improper_atoms1 = [
            tuple(a for a in atoms)
            for atoms in topology_impropers[: molecule1.n_impropers]
        ]
        top_improper_atoms2 = [
            tuple(a for a in atoms)
            for atoms in topology_impropers[
                molecule1.n_impropers : 2 * molecule1.n_impropers
            ]
        ]
        top_improper_atoms3 = [
            tuple(a for a in atoms)
            for atoms in topology_impropers[2 * molecule1.n_impropers :]
        ]

        assert_tuple_of_atoms_equal(
            top_improper_atoms1, mol_improper_atoms1, transformed_dict_cls=ImproperDict
        )
        assert_tuple_of_atoms_equal(
            top_improper_atoms2, mol_improper_atoms1, transformed_dict_cls=ImproperDict
        )
        assert_tuple_of_atoms_equal(
            top_improper_atoms3, mol_improper_atoms2, transformed_dict_cls=ImproperDict
        )

    def test_pruned_impropers(self):
        """Test {smirnoff|amber}_impropers from the Topology API"""
        top = Topology.from_molecules(
            [Molecule.from_smiles(smi) for smi in ["N", "C=C"]]
        )

        assert len([*top.smirnoff_impropers]) == 18
        assert len([*top.amber_impropers]) == 18

        # Order not guaranteed, so cannot zip and compare directly
        for smirnoff_imp in top.smirnoff_impropers:
            # Convert SMIRNOFF-style improper into AMBER-style
            mod_imp = (
                smirnoff_imp[1],
                smirnoff_imp[0],
                smirnoff_imp[2],
                smirnoff_imp[3],
            )
            assert mod_imp in top.amber_impropers

    # test_two_of_same_molecule
    # test_two_different_molecules
    # test_get_molecule
    # test_is_bonded
    # TODO: Test serialization

    def test_from_openmm(self):
        """Test creation of an OpenFF Topology object from an OpenMM Topology and component molecules"""
        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )

        with pytest.raises(
            MissingUniqueMoleculesError, match="requires a list of Molecule objects"
        ):
            Topology.from_openmm(pdbfile.topology)

        molecules = [create_ethanol(), create_cyclohexane()]

        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        assert topology.n_molecules == 239

    def test_from_openmm_missing_reference(self):
        """Test creation of an OpenFF Topology object from an OpenMM Topology when missing a unique molecule"""
        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )

        molecules = [create_ethanol()]
        with pytest.raises(ValueError, match="No match found for molecule C6H12"):
            Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

    def test_from_openmm_missing_conect(self):
        """
        Test creation of an OpenFF Topology object from an OpenMM Topology
        when the origin PDB lacks CONECT records
        """
        pdbfile = app.PDBFile(
            get_data_file_path("systems/test_systems/1_ethanol_no_conect.pdb")
        )

        molecules = []
        molecules.append(Molecule.from_smiles("CCO"))
        with pytest.raises(
            ValueError,
            match="No match found for molecule C. This would be a "
            "very unusual molecule to try and parameterize, "
            "and it is likely that the data source it was "
            "read from does not contain connectivity "
            "information. If this molecule is coming from "
            "PDB, please ensure that the file contains CONECT "
            "records.",
        ):
            Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

    def test_to_from_openmm(self):
        """Test a round-trip OpenFF -> OpenMM -> OpenFF Topology."""
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        benzene = Molecule.from_smiles("c1ccccc1")
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])

        # Convert to OpenMM Topology.
        omm_topology = off_topology.to_openmm()

        # Check that bond orders are preserved.
        n_double_bonds = sum([b.order == 2 for b in omm_topology.bonds()])
        n_aromatic_bonds = sum([b.type is app.Aromatic for b in omm_topology.bonds()])
        assert n_double_bonds == 6
        assert n_aromatic_bonds == 12

        # Check that there is one residue and chain for each molecule.
        assert omm_topology.getNumResidues() == 3
        assert omm_topology.getNumChains() == 3

        # Convert back to OpenFF Topology.
        off_topology_copy = Topology.from_openmm(
            omm_topology, unique_molecules=[ethanol, benzene]
        )

        # The round-trip OpenFF Topology is identical to the original.
        # The reference molecules are the same.
        assert off_topology.n_molecules == off_topology_copy.n_molecules
        reference_molecules_copy = list(off_topology_copy.reference_molecules)
        for ref_mol_idx, ref_mol in enumerate(off_topology.reference_molecules):
            assert ref_mol == reference_molecules_copy[ref_mol_idx]

        # The number of topology molecules is the same.
        assert off_topology.n_molecules == off_topology_copy.n_molecules

        # Check atoms.
        assert off_topology.n_atoms == off_topology_copy.n_atoms
        for atom_idx, atom in enumerate(off_topology.atoms):
            atom_copy = off_topology_copy.atom(atom_idx)
            assert atom.atomic_number == atom_copy.atomic_number

        # Check bonds.
        for bond_idx, bond in enumerate(off_topology.bonds):
            # bond_copy = off_topology_copy.bond(bond_idx)
            bond_copy = off_topology_copy.get_bond_between(
                off_topology.atom_index(bond.atoms[0]),
                off_topology.atom_index(bond.atoms[1]),
            )
            bond_atoms = [a.atomic_number for a in bond.atoms]
            bond_atoms_copy = [a.atomic_number for a in bond_copy.atoms]
            assert bond_atoms == bond_atoms_copy
            assert bond.bond_order == bond_copy.bond_order
            assert bond.is_aromatic == bond_copy.is_aromatic

    def test_to_from_openmm_hierarchy_metadata(self):
        """
        Test roundtripping to/from ``OpenEyeToolkitWrapper`` for molecules with PDB hierarchy metadata
        """
        top = topology_with_metadata()
        omm_top = top.to_openmm()

        # Make a list of the unique molecules in the topology for the return from roundtripping to OpenMM
        unique_mols = []
        for mol in top.molecules:
            if mol.to_smiles() not in [m.to_smiles() for m in unique_mols]:
                unique_mols.append(mol)

        roundtrip_top = Topology.from_openmm(omm_top, unique_molecules=unique_mols)

        # Check OMM Atom
        for orig_atom, omm_atom in zip(top.atoms, omm_top.atoms()):
            if "residue_name" in orig_atom.metadata:
                assert orig_atom.metadata["residue_name"] == omm_atom.residue.name
            else:
                assert omm_atom.residue.name == "UNK"

            if "residue_number" in orig_atom.metadata:
                assert orig_atom.metadata["residue_number"] == int(omm_atom.residue.id)
            else:
                assert omm_atom.residue.id == "0"

            if "chain_id" in orig_atom.metadata:
                assert orig_atom.metadata["chain_id"] == omm_atom.residue.chain.id
            else:
                assert omm_atom.residue.chain.id == "X"

        # Check roundtripped OFFMol
        for orig_atom, roundtrip_atom in zip(top.atoms, roundtrip_top.atoms):
            if "residue_name" in orig_atom.metadata:
                original = orig_atom.metadata["residue_name"]
                roundtrip = roundtrip_atom.metadata["residue_name"]
                assert original == roundtrip
            else:
                assert roundtrip_atom.metadata["residue_name"] == "UNK"

            if "residue_number" in orig_atom.metadata:
                original = orig_atom.metadata["residue_number"]
                roundtrip = roundtrip_atom.metadata["residue_number"]
                assert original == roundtrip
            else:
                assert roundtrip_atom.metadata["residue_number"] == 0

            if "chain_id" in orig_atom.metadata:
                original = orig_atom.metadata["chain_id"]
                roundtrip = roundtrip_atom.metadata["chain_id"]
                assert original == roundtrip
            else:
                assert roundtrip_atom.metadata["chain_id"] == "X"

    @requires_pkg("mdtraj")
    def test_from_mdtraj(self):
        """Test construction of an OpenFF Topology from an MDTraj Topology object"""
        import mdtraj as md

        pdb_path = get_data_file_path(
            "systems/test_systems/1_cyclohexane_1_ethanol.pdb"
        )
        trj = md.load(pdb_path)

        with pytest.raises(
            MissingUniqueMoleculesError, match="requires a list of Molecule objects"
        ):
            Topology.from_mdtraj(trj.top)

        unique_molecules = [
            Molecule.from_smiles(mol_name) for mol_name in ["C1CCCCC1", "CCO"]
        ]
        top = Topology.from_mdtraj(trj.top, unique_molecules=unique_molecules)

        assert top.n_molecules == 2
        assert top.n_bonds == 26

    @requires_rdkit
    def test_to_file_units_check(self):
        """
        Checks that writing a PDB file with different coordinate representations results in the same output.
        - Angstrom "openff units" (default behavior if using Molecule.conformers[0])
        - nanometer "openff units"
        - unitless NumPy array
        - converted OpenMM quantity
        """
        from tempfile import NamedTemporaryFile

        from openff.units.openmm import to_openmm

        from openff.toolkit.topology import Molecule, Topology

        topology = Topology()
        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"), "CCO"
        )
        topology.add_molecule(mol)
        positions_angstrom = mol.conformers[0]

        def _check_file(topology, coordinates):
            # Write the molecule to PDB and ensure that the X coordinate of the first atom is 10.172
            count = 1
            with NamedTemporaryFile(suffix=".pdb") as iofile:
                topology.to_file(iofile.name, coordinates)
                data = open(iofile.name).readlines()
                for line in data:
                    if line.startswith("HETATM") and count == 1:
                        count = count + 1
                        coord = line.split()[-6]
            assert coord == "10.172"

        _check_file(topology, coordinates=positions_angstrom)
        _check_file(topology, coordinates=positions_angstrom.to(unit.nanometer))
        _check_file(topology, coordinates=positions_angstrom.m)
        _check_file(topology, coordinates=to_openmm(positions_angstrom))

        with pytest.raises(ValueError, match="Could not process.*list.*"):
            _check_file(topology, coordinates=positions_angstrom.m.tolist())

    @requires_rdkit
    def test_to_file_fileformat_lettercase(self):
        """
        Checks if fileformat specifier is indpendent of upper/lowercase
        """
        from tempfile import NamedTemporaryFile

        from openff.toolkit.topology import Molecule, Topology

        topology = Topology()
        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"), "CCO"
        )
        topology.add_molecule(mol)
        positions = mol.conformers[0]
        count = 1
        with NamedTemporaryFile(suffix=".pdb") as iofile:
            topology.to_file(iofile.name, positions, file_format="pDb")
            data = open(iofile.name).readlines()
            for line in data:
                if line.startswith("HETATM") and count == 1:
                    count = count + 1
                    coord = line.split()[-6]
        assert coord == "10.172"

    @requires_rdkit
    def test_to_file_fileformat_invalid(self):
        """
        Checks for invalid file format
        """
        from openff.toolkit.topology import Molecule, Topology

        topology = Topology()
        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"), "CCO"
        )
        topology.add_molecule(mol)
        positions = mol.conformers[0]
        fname = "ethanol_file.pdb"
        with pytest.raises(NotImplementedError):
            topology.to_file(fname, positions, file_format="AbC")

    def test_to_file_no_molecules(self):
        """
        Checks if Topology.to_file() writes a file with no topology and no coordinates
        """
        from tempfile import NamedTemporaryFile

        from openff.toolkit.topology import Topology

        topology = Topology()
        lines = []
        with NamedTemporaryFile(suffix=".pdb") as iofile:
            topology.to_file(iofile.name, [] * unit.nanometer)
            data = open(iofile.name).readlines()
            for line in data:
                lines.append(line.split())
        assert lines[1] == ["END"]

    @requires_rdkit
    def test_to_file_multi_molecule_different_order(self):
        """
        Checks for the following if Topology.to_write maintains the order of atoms
         for the same molecule with different indexing
        """
        from tempfile import NamedTemporaryFile

        from openff.toolkit.topology import Molecule, Topology

        topology = Topology()
        topology.add_molecule(create_ethanol())
        topology.add_molecule(create_reversed_ethanol())
        mol = Molecule.from_pdb_and_smiles(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"), "CCO"
        )
        positions = mol.conformers[0]
        # Make up coordinates for the second ethanol by translating the first by 10 angstroms (note that this will
        # still be a gibberish conformation, since the atom order in the second molecule is different)
        positions = np.concatenate([positions, positions + 10.0 * unit.angstrom])
        element_order = []

        with NamedTemporaryFile(suffix=".pdb") as iofile:
            topology.to_file(iofile.name, positions)
            data = open(iofile.name).readlines()
            for line in data:
                if line.startswith("HETATM"):
                    element_order.append(line.strip()[-1])
        assert element_order == [
            "C",
            "C",
            "O",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "O",
            "C",
            "C",
        ]

    @requires_openeye
    def test_from_openmm_duplicate_unique_mol(self):
        """
        Check that a DuplicateUniqueMoleculeError is raised if we try to pass in two indistinguishably unique mols
        """
        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )
        molecules = [
            Molecule.from_file(get_data_file_path(name))
            for name in (
                "molecules/ethanol.mol2",
                "molecules/ethanol_reordered.mol2",
                "molecules/cyclohexane.mol2",
            )
        ]
        with pytest.raises(DuplicateUniqueMoleculeError):
            Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

    @pytest.mark.skip
    def test_from_openmm_distinguish_using_stereochemistry(self):
        """Test creation of an OpenFF Topology object from an OpenMM topology with stereoisomers"""
        # From Jeff: The graph representation created from OMM molecules during the matching
        # process doesn't encode stereochemistry.
        raise NotImplementedError

    def test_to_openmm_assign_unique_atom_names(self):
        """
        Ensure that OFF topologies with no pre-existing atom names have unique
        atom names applied when being converted to openmm
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        benzene = Molecule.from_smiles("c1ccccc1")
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])
        omm_topology = off_topology.to_openmm()
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        # There should be 6 unique Cs, 6 unique Hs, and 1 unique O, for a total of 13 unique atom names
        assert len(atom_names) == 13

    def test_to_openmm_assign_some_unique_atom_names(self):
        """
        Ensure that OFF topologies with some pre-existing atom names have unique
        atom names applied to the other atoms when being converted to openmm
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        for atom in ethanol.atoms:
            atom.name = f"AT{atom.molecule_atom_index}"
        benzene = Molecule.from_smiles("c1ccccc1")
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])
        omm_topology = off_topology.to_openmm()
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        # There should be 9 "ATOM#"-labeled atoms, 6 unique Cs, and 6 unique Hs,
        # for a total of 21 unique atom names
        assert len(atom_names) == 21

    def test_to_openmm_assign_unique_atom_names_some_duplicates(self):
        """
        Ensure that OFF topologies where some molecules have invalid/duplicate
        atom names have unique atom names applied while the other molecules are unaffected.
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")

        # Assign duplicate atom names in ethanol (two AT0s)
        ethanol_atom_names_with_duplicates = [f"AT{i}" for i in range(ethanol.n_atoms)]
        ethanol_atom_names_with_duplicates[1] = "AT0"
        for atom, atom_name in zip(ethanol.atoms, ethanol_atom_names_with_duplicates):
            atom.name = atom_name

        # Assign unique atom names in benzene
        benzene = Molecule.from_smiles("c1ccccc1")
        benzene_atom_names = [f"AT{i}" for i in range(benzene.n_atoms)]
        for atom, atom_name in zip(benzene.atoms, benzene_atom_names):
            atom.name = atom_name

        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])
        omm_topology = off_topology.to_openmm()
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)

        # There should be  12 "AT#"-labeled atoms (from benzene), 2 unique Cs,
        # 1 unique O, and 6 unique Hs, for a total of 21 unique atom names
        assert len(atom_names) == 21

    def test_to_openmm_do_not_assign_unique_atom_names(self):
        """
        Test disabling unique atom name assignment in Topology.to_openmm
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        for atom in ethanol.atoms:
            atom.name = "eth_test"
        benzene = Molecule.from_smiles("c1ccccc1")
        benzene.atoms[0].name = "bzn_test"
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])
        omm_topology = off_topology.to_openmm(ensure_unique_atom_names=False)
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        # There should be 9 atom named "eth_test", 1 atom named "bzn_test",
        # and 12 atoms named "", for a total of 3 unique atom names
        assert len(atom_names) == 3

    def test_group_chemically_identical_molecules(self):
        """Test behavior and caching of Topology.group_chemically_identical_molecules"""
        top = Topology()

        # Test for correct behavior with empty topology
        assert top._identify_chemically_identical_molecules() == {}

        # Test for correct behavior with topology of one ethanol
        top.add_molecule(create_ethanol())
        groupings = top.identical_molecule_groups

        def assert_first_ethanol_is_grouped_correctly(groupings):
            assert groupings[0][0] == [0, {i: i for i in range(9)}]

        assert_first_ethanol_is_grouped_correctly(groupings)

        # Add an ethanol in reversed order
        top.add_molecule(create_reversed_ethanol())
        groupings = top.identical_molecule_groups

        def assert_reversed_ethanol_is_grouped_correctly(groupings):
            # Ensure that the second ethanol knows it's the same chemical species as the first ethanol
            assert groupings[0][1][0] == 1
            # Ensure that the second ethanol has the heavy atoms reversed
            assert groupings[0][1][1][0] == 8  # C
            assert groupings[0][1][1][1] == 7  # C
            assert groupings[0][1][1][2] == 6  # O
            # (we only check the hydroxyl H, since the other Hs have multiple valid matches)
            assert groupings[0][1][1][8] == 0  # HO

        assert_first_ethanol_is_grouped_correctly(groupings)
        assert_reversed_ethanol_is_grouped_correctly(groupings)

        # Add a cyclohexane, which should be unique from all the other molecules
        top.add_molecule(create_cyclohexane())

        def assert_cyclohexane_is_grouped_correctly(groupings):
            assert len(groupings[2]) == 1
            assert groupings[2][0] == [2, {i: i for i in range(18)}]

        groupings = top.identical_molecule_groups
        assert_first_ethanol_is_grouped_correctly(groupings)
        assert_reversed_ethanol_is_grouped_correctly(groupings)
        assert_cyclohexane_is_grouped_correctly(groupings)

        # Add a third ethanol, in the same order as the first
        top.add_molecule(create_ethanol())

        def assert_last_ethanol_is_grouped_correctly(groupings):
            # Ensure that the last ethanol knows it's the same chemical species as the first ethanol
            assert groupings[0][2][0] == 3
            # Ensure that the second ethanol has the heavy atoms matched correctly
            assert groupings[0][2][1][0] == 0
            assert groupings[0][2][1][1] == 1
            assert groupings[0][2][1][2] == 2
            # Ensure that the second ethanol has the hydroxyl hydrogen matched correctly
            assert groupings[0][2][1][8] == 8

        groupings = top.identical_molecule_groups
        assert_first_ethanol_is_grouped_correctly(groupings)
        assert_reversed_ethanol_is_grouped_correctly(groupings)
        assert_cyclohexane_is_grouped_correctly(groupings)
        assert_last_ethanol_is_grouped_correctly(groupings)

    @requires_openeye
    def test_chemical_environments_matches_OE(self):
        """Test Topology.chemical_environment_matches"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )
        # toolkit_wrapper = RDKitToolkitWrapper()
        molecules = [
            Molecule.from_file(get_data_file_path(name))
            for name in ("molecules/ethanol.mol2", "molecules/cyclohexane.mol2")
        ]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        # Test for substructure match
        matches = topology.chemical_environment_matches(
            "[C:1]-[C:2]-[O:3]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 143
        assert matches[0].topology_atom_indices == (1728, 1729, 1730)
        # Test for whole-molecule match
        matches = topology.chemical_environment_matches(
            "[H][C:1]([H])([H])-[C:2]([H])([H])-[O:3][H]",
            toolkit_registry=toolkit_wrapper,
        )
        assert (
            len(matches) == 1716
        )  # 143 * 12 (there are 12 possible hydrogen mappings)
        assert matches[0].topology_atom_indices == (1728, 1729, 1730)
        # Search for a substructure that isn't there
        matches = topology.chemical_environment_matches(
            "[C][C:1]-[C:2]-[O:3]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 0

    @requires_rdkit
    def test_chemical_environments_matches_RDK(self):
        """Test Topology.chemical_environment_matches"""
        toolkit_wrapper = RDKitToolkitWrapper()
        pdbfile = app.PDBFile(
            get_data_file_path("systems/packmol_boxes/cyclohexane_ethanol_0.4_0.6.pdb")
        )
        # toolkit_wrapper = RDKitToolkitWrapper()
        # molecules = [Molecule.from_file(get_data_file_path(name)) for name in ('molecules/ethanol.mol2',
        #                                                                      'molecules/cyclohexane.mol2')]
        molecules = []
        molecules.append(Molecule.from_smiles("CCO"))
        molecules.append(Molecule.from_smiles("C1CCCCC1"))
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
        # Count CCO matches
        matches = topology.chemical_environment_matches(
            "[C:1]-[C:2]-[O:3]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 143
        assert matches[0].topology_atom_indices == (1728, 1729, 1730)
        matches = topology.chemical_environment_matches(
            "[H][C:1]([H])([H])-[C:2]([H])([H])-[O:3][H]",
            toolkit_registry=toolkit_wrapper,
        )
        assert (
            len(matches) == 1716
        )  # 143 * 12 (there are 12 possible hydrogen mappings)
        assert matches[0].topology_atom_indices == (1728, 1729, 1730)
        # Search for a substructure that isn't there
        matches = topology.chemical_environment_matches(
            "[C][C:1]-[C:2]-[O:3]", toolkit_registry=toolkit_wrapper
        )
        assert len(matches) == 0

    def test_topology_hierarchy_iterators(
        self,
    ):
        top = Topology()
        # Ensure that an empty topology has no residues defined
        residues = list(top.hierarchy_iterator("residues"))
        assert residues == []
        # Ensure that a topology with no metadata has no residues defined
        top.add_molecule(dipeptide())
        residues = list(top.hierarchy_iterator("residues"))
        assert residues == []
        # Ensure that a topology with residues perceived has residues but not chains
        top.add_molecule(dipeptide_residues_perceived())
        residues = list(top.hierarchy_iterator("residues"))
        chains = list(top.hierarchy_iterator("chains"))
        assert chains == []
        assert [res.identifier for res in residues] == [
            ("None", 1, "ACE"),
            ("None", 2, "ALA"),
        ]
        # Ensure that adding molecules WITH hierarchy perceived DOES give the
        # topology residues and chains to iterate over
        top.add_molecule(dipeptide_hierarchy_added())
        top.add_molecule(cyx_hierarchy_added())
        residues = list(top.hierarchy_iterator("residues"))
        chains = list(top.hierarchy_iterator("chains"))
        assert [res.identifier for res in residues] == [
            ("None", 1, "ACE"),
            ("None", 2, "ALA"),
            ("None", 1, "ACE"),
            ("None", 2, "ALA"),
            ("None", 4, "ACE"),
            ("None", 5, "CYS"),
            ("None", 6, "NME"),
        ]
        assert len(chains) == 1
        # Haven't changed anything, so updating hierarchy schemes should give
        # same results
        top.molecule(3).update_hierarchy_schemes()
        residues = list(top.hierarchy_iterator("residues"))
        chains = list(top.hierarchy_iterator("chains"))
        assert [res.identifier for res in residues] == [
            ("None", 1, "ACE"),
            ("None", 2, "ALA"),
            ("None", 1, "ACE"),
            ("None", 2, "ALA"),
            ("None", 4, "ACE"),
            ("None", 5, "CYS"),
            ("None", 6, "NME"),
        ]
        assert len(chains) == 1


class TestAddTopology:
    def test_add_basic(self):
        topology1 = Molecule.from_smiles("O").to_topology()
        topology2 = Molecule.from_smiles("CO").to_topology()

        topology3 = topology1 + topology2

        assert topology3.n_atoms == 9
        assert topology3.n_bonds == 7
        assert topology3.n_molecules == 2

    def test_add_inplace(self):
        topology1 = Molecule.from_smiles("O").to_topology()
        topology2 = Molecule.from_smiles("CO").to_topology()

        topology3 = topology1 + topology2
        topology1 += topology2

        for attr in ["n_atoms", "n_bonds", "n_molecules"]:
            assert getattr(topology1, attr) == getattr(topology3, attr)

    def test_add_invalidate_cache(self):
        topology1 = Molecule.from_smiles("O").to_topology()
        topology2 = Molecule.from_smiles("CO").to_topology()

        topology1.add_constraint(0, 1, 1.01 * unit.angstrom)
        topology2.identical_molecule_groups

        assert topology1.constrained_atom_pairs[(0, 1)] == 1.01 * unit.angstrom
        assert len(topology2._cached_chemically_identical_molecules) == 1

        topology3 = topology1 + topology2

        assert topology3._cached_chemically_identical_molecules is None

        # Constrained atom pairs are intentionally not removed
        assert topology3.constrained_atom_pairs[(0, 1)] == 1.01 * unit.angstrom

    def test_add_with_constraints(self):
        # See https://github.com/openforcefield/openff-toolkit/pull/1194#discussion_r834768068
        methane = Molecule.from_mapped_smiles(
            "[H:2][C:1]([H:3])([H:4])[H:5]"
        ).to_topology()
        water = Molecule.from_mapped_smiles("[H:2][O:1][H:3]").to_topology()

        # Constrain one C-H bond in methane and H-H in water (0-indexed)
        methane.add_constraint(0, 1)
        water.add_constraint(1, 2, unit.Quantity(1.234, unit.angstrom))

        combined_topology = methane + water

        # Constraints are tracked i-j and j-i, so this dict is length 4
        assert len(combined_topology.constrained_atom_pairs) == 4

        # Methane is the first molecule, so (0, 1) is conserved. The atoms in the constrained
        # pair in water are offset by the number of atoms in methane (5)
        assert ((0, 1)) in combined_topology.constrained_atom_pairs
        assert ((1, 0)) in combined_topology.constrained_atom_pairs
        assert ((6, 7)) in combined_topology.constrained_atom_pairs
        assert ((7, 6)) in combined_topology.constrained_atom_pairs

        assert combined_topology.constrained_atom_pairs[(0, 1)]
        assert combined_topology.constrained_atom_pairs[(1, 0)]
        assert combined_topology.constrained_atom_pairs[(6, 7)] == unit.Quantity(
            1.234, unit.angstrom
        )
        assert combined_topology.constrained_atom_pairs[(7, 6)] == unit.Quantity(
            1.234, unit.angstrom
        )


class TestTopologySerialization:
    @pytest.fixture
    def oleic_acid(self):
        """Simple floppy molecule that can be assured to have multiple conformers"""
        return Molecule.from_smiles(r"CCCCCCCC\C=C/CCCCCCCC(O)=O")

    @pytest.mark.parametrize(("with_conformers"), [True, False])
    @pytest.mark.parametrize(("n_molecules"), [1, 2])
    @pytest.mark.parametrize(("format"), ["dict", "json"])
    def test_roundtrip(self, oleic_acid, with_conformers, n_molecules, format):

        if with_conformers:
            n_conformers = 2
            oleic_acid.generate_conformers(n_conformers=n_conformers)
        else:
            n_conformers = 0

        if format == "dict":
            roundtrip = Topology.from_dict(
                Topology.from_molecules(n_molecules * [oleic_acid]).to_dict()
            )
        elif format == "json":
            roundtrip = Topology.from_json(
                Topology.from_molecules(n_molecules * [oleic_acid]).to_json()
            )

        assert roundtrip.n_molecules == n_molecules
        assert roundtrip.n_atoms == oleic_acid.n_atoms * n_molecules
        assert [*roundtrip.molecules][0].n_conformers == n_conformers


@pytest.mark.parametrize(
    ("n_degrees", "num_pairs"),
    [
        (6, 0),
        (5, 3),
        (4, 14),
        (3, 28),
    ],
)
def test_nth_degree_neighbors(n_degrees, num_pairs):
    smiles = ["c1ccccc1", "N1ONON1"]
    topology = Topology.from_molecules([Molecule.from_smiles(smi) for smi in smiles])

    # See test_molecule.TestMolecule.test_nth_degree_neighbors_rings for values
    num_pairs_found = len([*topology.nth_degree_neighbors(n_degrees=n_degrees)])
    assert num_pairs_found == num_pairs


def _tagsorted_dict_init_ref_key(tsd):

    if tsd is None:
        tsd = TagSortedDict()

        ref_key = (0, 1, 2)
        tsd[ref_key] = 5
    else:
        ref_key = next(x for x in tsd)

    return tsd, ref_key


@pytest.mark.parametrize("tsd", [None, TagSortedDict({(0, 1, 2): 5})])
def test_tagsorted_dict_deduplication(tsd):
    """Test that all permutations of a key are present if one permutation is stored"""
    tsd, ref_key = _tagsorted_dict_init_ref_key(tsd)

    # only has one key, but all permutations match
    for key in itertools.permutations(ref_key):
        assert key in tsd

    assert len(tsd) == 1


@pytest.mark.parametrize("tsd", [None, TagSortedDict({(0, 1, 2): 5})])
def test_tagsorted_dict_permutation_equivalence(tsd):
    """Test that all permutations of a key would return the same result"""
    tsd, ref_key = _tagsorted_dict_init_ref_key(tsd)

    # lookups using any permutation will give the same return
    for key in itertools.permutations(ref_key):
        assert tsd[key] == 5

    assert len(tsd) == 1


@pytest.mark.parametrize("tsd", [None, TagSortedDict({(0, 1, 2): 5})])
def test_tagsorted_dict_key_transform(tsd):
    """Test that all key permutations transform to the same stored single key
    permutation"""
    tsd, ref_key = _tagsorted_dict_init_ref_key(tsd)

    # all permutations should resolve to ref_key
    for key in itertools.permutations(ref_key):
        tr_key = tsd.key_transform(key)
        assert tr_key == ref_key

    assert len(tsd) == 1


@pytest.mark.parametrize("tsd", [None, TagSortedDict({(0, 1, 2): 5})])
def test_tagsorted_dict_modify(tsd):
    """Test that modifying a key with another permutation replaces the original key"""
    tsd, ref_key = _tagsorted_dict_init_ref_key(tsd)

    # replace the ref_key since this is a permutation of it
    mod_key = (1, 0, 2)
    tsd[mod_key] = 6

    keys = list(tsd)
    assert ref_key not in keys

    for key in itertools.permutations(ref_key):
        assert key in tsd

    assert len(tsd) == 1  # should not be 2


@pytest.mark.parametrize("tsd", [TagSortedDict({(0, 1, 2): 5, (1, 2): 4})])
def test_tagsorted_dict_multiple_keys(tsd):
    """Test the use of multiple keys with similar values but different length"""
    tsd, ref_key = _tagsorted_dict_init_ref_key(tsd)

    # replace the ref_key since this is a permutation of it
    mod_key = (1, 0, 2)
    tsd[mod_key] = 6

    keys = list(tsd)
    assert ref_key not in keys  # ensure original key is gone

    for key in itertools.permutations(ref_key):
        assert key in tsd

    assert tsd[(2, 1)] == 4  # should be unchanged
    assert len(tsd) == 2  # should not be 3


@pytest.mark.parametrize("tsd", [TagSortedDict({(0, 1, 2): 5, (1, 2): 4})])
def test_tagsorted_dict_clear(tsd):
    """Test the clear method"""
    tsd, ref_key = _tagsorted_dict_init_ref_key(tsd)

    tsd.clear()

    assert len(tsd) == 0
