from copy import deepcopy

import numpy
import pytest
from openff.toolkit import Molecule, Topology
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest, needs_gmx
from openff.interchange.components.mdconfig import get_intermol_defaults
from openff.interchange.drivers.gromacs import _process, _run_gmx_energy
from openff.interchange.smirnoff._gromacs import _convert


@pytest.fixture()
def molecule1():
    molecule = Molecule.from_smiles(
        "[H][O][c]1[c]([H])[c]([O][H])[c]([H])[c]([O][H])[c]1[H]",
    )
    molecule.generate_conformers(n_conformers=1)
    molecule.name = "MOL1"

    return molecule


@pytest.fixture()
def molecule2():
    molecule = Molecule.from_smiles("C1=C(C=C(C=C1C(=O)O)C(=O)O)C(=O)O")
    molecule.generate_conformers(n_conformers=1)
    molecule.name = "MOL2"

    molecule.conformers[0] += numpy.array([5, 0, 0]) * unit.angstrom

    return molecule


@pytest.fixture()
def system1(molecule1, sage):
    box = 5 * numpy.eye(3) * unit.nanometer

    return _convert(Interchange.from_smirnoff(sage, [molecule1], box=box))


@pytest.fixture()
def system2(molecule2, sage):
    box = 5 * numpy.eye(3) * unit.nanometer

    return _convert(Interchange.from_smirnoff(sage, [molecule2], box=box))


@pytest.fixture()
def combined_system(molecule1, molecule2, sage):
    box = 5 * numpy.eye(3) * unit.nanometer

    return _convert(
        Interchange.from_smirnoff(
            sage,
            Topology.from_molecules([molecule1, molecule2]),
            box=box,
        ),
    )


class TestAddRemoveMoleculeType(_BaseTest):
    @needs_gmx
    @pytest.mark.parametrize("molecule_name", ["MOL1", "MOL2"])
    def test_remove_basic(self, combined_system, molecule_name):
        combined_system.remove_molecule_type(molecule_name)

        # Just a sanity check
        combined_system.to_files(prefix=molecule_name, decimal=8)

        get_intermol_defaults(periodic=True).write_mdp_file("tmp.mdp")

        _process(
            _run_gmx_energy(f"{molecule_name}.top", f"{molecule_name}.gro", "tmp.mdp"),
        )

    @pytest.mark.slow()
    @pytest.mark.parametrize("molecule_name", ["MOL1", "MOL2"])
    def test_add_existing_molecule_type(self, combined_system, molecule_name):
        with pytest.raises(
            ValueError,
            match=f"The molecule type {molecule_name} is already present in this system.",
        ):
            combined_system.add_molecule_type(
                combined_system.molecule_types[molecule_name],
                1,
            )

    @needs_gmx
    def test_different_force_field_different_energies(
        self,
        combined_system,
        system1,
        molecule2,
        sage,
        parsley,
    ):
        box = 5 * numpy.eye(3) * unit.nanometer

        parsley_system = deepcopy(system1)
        sage_system = deepcopy(system1)

        molecule2_parsley = _convert(
            Interchange.from_smirnoff(parsley, [molecule2], box=box),
        )
        molecule2_sage = _convert(Interchange.from_smirnoff(sage, [molecule2], box=box))

        parsley_system.add_molecule_type(molecule2_parsley.molecule_types["MOL2"], 1)
        sage_system.add_molecule_type(molecule2_sage.molecule_types["MOL2"], 1)

        parsley_system.positions = combined_system.positions
        sage_system.positions = combined_system.positions

        parsley_system.to_files(prefix="parsley", decimal=8)
        sage_system.to_files(prefix="sage", decimal=8)

        get_intermol_defaults(periodic=True).write_mdp_file("tmp.mdp")

        _parsley_energy = _process(
            _run_gmx_energy("parsley.top", "parsley.gro", "tmp.mdp"),
        )
        _sage_energy = _process(_run_gmx_energy("sage.top", "sage.gro", "tmp.mdp"))

        assert _parsley_energy != _sage_energy

    @needs_gmx
    def test_molecule_order_independent(self, system1, system2):
        positions1 = numpy.vstack([system1.positions, system2.positions])
        positions2 = numpy.vstack([system2.positions, system1.positions])

        system1.add_molecule_type(system2.molecule_types["MOL2"], 1)
        system1.positions = positions1

        system2.add_molecule_type(system1.molecule_types["MOL1"], 1)
        system2.positions = positions2

        system1.to_files(prefix="order1", decimal=8)
        system2.to_files(prefix="order2", decimal=8)

        get_intermol_defaults(periodic=True).write_mdp_file("tmp.mdp")

        _process(_run_gmx_energy("order1.top", "order1.gro", "tmp.mdp")).compare(
            _process(_run_gmx_energy("order2.top", "order2.gro", "tmp.mdp")),
        )

    @pytest.mark.slow()
    def test_clashing_atom_types(self, combined_system, system1, system2):
        with pytest.raises(
            ValueError,
            match="The molecule type MOL1 is already present in this system.",
        ):
            combined_system.add_molecule_type(system1.molecule_types["MOL1"], 1)

        with pytest.raises(
            ValueError,
            match="The molecule type MOL2 is already present in this system.",
        ):
            combined_system.add_molecule_type(system2.molecule_types["MOL2"], 1)

        with pytest.raises(
            ValueError,
            match="The molecule type MOL1 is already present in this system.",
        ):
            system1.add_molecule_type(system1.molecule_types["MOL1"], 1)

        with pytest.raises(
            ValueError,
            match="The molecule type MOL2 is already present in this system.",
        ):
            system2.add_molecule_type(system2.molecule_types["MOL2"], 1)


class TestToFiles(_BaseTest):
    @needs_gmx
    def test_identical_outputs(self, system1):
        system1.to_files(prefix="1", decimal=8)

        system1.to_top("2.top")
        system1.to_gro("2.gro", decimal=8)

        get_intermol_defaults(periodic=True).write_mdp_file("tmp.mdp")

        _process(_run_gmx_energy("1.top", "1.gro", "tmp.mdp")).compare(
            _process(_run_gmx_energy("2.top", "2.gro", "tmp.mdp")),
        )
