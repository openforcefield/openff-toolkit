import itertools

import numpy
import openmm
import pytest
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff.forcefield import ForceField
from openff.toolkit.typing.engines.smirnoff.parameters import (
    ElectrostaticsHandler,
    LibraryChargeHandler,
    VirtualSiteHandler,
    vdWHandler,
)
from openff.units import unit
from openff.units.openmm import to_openmm
from openmm import unit as openmm_unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest


def _get_interpolated_bond_k(bond_handler) -> float:
    for key in bond_handler.key_map:
        if key.bond_order is not None:
            topology_key = key
            break
    potential_key = bond_handler.key_map[topology_key]
    return bond_handler.potentials[potential_key].parameters["k"].m


class TestSMIRNOFFVirtualSites(_BaseTest):
    from openff.toolkit.tests.mocking import VirtualSiteMocking
    from openmm import unit as openmm_unit

    @classmethod
    def build_force_field(cls, v_site_handler: VirtualSiteHandler) -> ForceField:
        force_field = ForceField()

        force_field.get_parameter_handler("Bonds").add_parameter(
            {
                "smirks": "[*:1]~[*:2]",
                "k": 0.0 * unit.kilojoule_per_mole / unit.angstrom**2,
                "length": 0.0 * unit.angstrom,
            },
        )
        force_field.get_parameter_handler("Angles").add_parameter(
            {
                "smirks": "[*:1]~[*:2]~[*:3]",
                "k": 0.0 * unit.kilojoule_per_mole / unit.degree**2,
                "angle": 60.0 * unit.degree,
            },
        )
        force_field.get_parameter_handler("ProperTorsions").add_parameter(
            {
                "smirks": "[*:1]~[*:2]~[*:3]~[*:4]",
                "k": [0.0] * unit.kilojoule_per_mole,
                "phase": [0.0] * unit.degree,
                "periodicity": [2],
                "idivf": [1.0],
            },
        )
        force_field.get_parameter_handler("vdW").add_parameter(
            {
                "smirks": "[*:1]",
                "epsilon": 0.0 * unit.kilojoule_per_mole,
                "sigma": 1.0 * unit.angstrom,
            },
        )
        force_field.get_parameter_handler("LibraryCharges").add_parameter(
            {"smirks": "[*:1]", "charge": [0.0] * unit.elementary_charge},
        )
        force_field.get_parameter_handler("Electrostatics")

        force_field.register_parameter_handler(v_site_handler)

        return force_field

    @classmethod
    def generate_v_site_coordinates(
        cls,
        molecule: Molecule,
        input_conformer: unit.Quantity,
        parameter: VirtualSiteHandler.VirtualSiteType,
    ) -> unit.Quantity:
        # Compute the coordinates of the virtual site. Unfortunately OpenMM does not
        # seem to offer a more compact way to do this currently.
        handler = VirtualSiteHandler(version="0.3")
        handler.add_parameter(parameter=parameter)

        force_field = cls.build_force_field(handler)

        system = Interchange.from_smirnoff(
            force_field=force_field,
            topology=molecule.to_topology(),
        ).to_openmm(combine_nonbonded_forces=True)

        n_v_sites = sum(
            1 if system.isVirtualSite(i) else 0 for i in range(system.getNumParticles())
        )

        input_conformer = unit.Quantity(
            numpy.vstack(
                [
                    input_conformer.m_as(unit.angstrom),
                    numpy.zeros((n_v_sites, 3)),
                ],
            ),
            unit.angstrom,
        )

        context = openmm.Context(
            system,
            openmm.VerletIntegrator(1.0 * openmm_unit.femtosecond),
            openmm.Platform.getPlatformByName("Reference"),
        )

        context.setPositions(to_openmm(input_conformer))
        context.computeVirtualSites()

        output_conformer = context.getState(getPositions=True).getPositions(
            asNumpy=True,
        )

        return output_conformer[molecule.n_atoms :, :]

    @pytest.mark.skip(
        "Interchange does not allow this combination of non-bonded settings.",
    )
    @pytest.mark.parametrize(
        (
            "parameter",
            "smiles",
            "input_conformer",
            "atoms_to_shuffle",
            "expected_coordinates",
        ),
        [
            (
                VirtualSiteMocking.bond_charge_parameter("[Cl:1]-[C:2]"),
                "[Cl:1][C:2]([H:3])([H:4])[H:5]",
                VirtualSiteMocking.sp3_conformer(),
                (0, 1),
                unit.Quantity(numpy.array([[0.0, 3.0, 0.0]]), unit.angstrom),
            ),
            (
                VirtualSiteMocking.bond_charge_parameter("[C:1]#[C:2]"),
                "[H:1][C:2]#[C:3][H:4]",
                VirtualSiteMocking.sp1_conformer(),
                (2, 3),
                unit.Quantity(
                    numpy.array([[-3.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
                    unit.angstrom,
                ),
            ),
            (
                VirtualSiteMocking.monovalent_parameter("[O:1]=[C:2]-[H:3]"),
                "[O:1]=[C:2]([H:3])[H:4]",
                VirtualSiteMocking.sp2_conformer(),
                (0, 1, 2, 3),
                (
                    VirtualSiteMocking.sp2_conformer()[0]
                    + unit.Quantity(  # noqa
                        numpy.array(
                            [[1.0, numpy.sqrt(2), 1.0], [1.0, -numpy.sqrt(2), -1.0]],
                        ),
                        unit.angstrom,
                    )
                ),
            ),
            (
                VirtualSiteMocking.divalent_parameter(
                    "[H:2][O:1][H:3]",
                    match="once",
                    angle=0.0 * unit.degree,
                ),
                "[H:2][O:1][H:3]",
                VirtualSiteMocking.sp2_conformer()[1:, :],
                (0, 1, 2),
                numpy.array([[2.0, 0.0, 0.0]]) * unit.angstrom,
            ),
            (
                VirtualSiteMocking.divalent_parameter(
                    "[H:2][O:1][H:3]",
                    match="all_permutations",
                    angle=45.0 * unit.degree,
                ),
                "[H:2][O:1][H:3]",
                VirtualSiteMocking.sp2_conformer()[1:, :],
                (0, 1, 2),
                unit.Quantity(
                    numpy.array(
                        [
                            [numpy.sqrt(2), numpy.sqrt(2), 0.0],
                            [numpy.sqrt(2), -numpy.sqrt(2), 0.0],
                        ],
                    ),
                    unit.angstrom,
                ),
            ),
            (
                VirtualSiteMocking.trivalent_parameter("[N:1]([H:2])([H:3])[H:4]"),
                "[N:1]([H:2])([H:3])[H:4]",
                VirtualSiteMocking.sp3_conformer()[1:, :],
                (0, 1, 2, 3),
                numpy.array([[0.0, 2.0, 0.0]]) * unit.angstrom,
            ),
        ],
    )
    def test_v_site_geometry(
        self,
        parameter: VirtualSiteHandler.VirtualSiteType,
        smiles: str,
        input_conformer: unit.Quantity,
        atoms_to_shuffle: tuple[int, ...],
        expected_coordinates: unit.Quantity,
    ):
        """An integration test that virtual sites are placed correctly relative to the
        parent atoms"""

        for atom_permutation in itertools.permutations(atoms_to_shuffle):
            molecule = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
            molecule._conformers = [input_conformer]

            # We shuffle the ordering of the atoms involved in the v-site orientation
            # to ensure that the orientation remains invariant as expected.
            shuffled_atom_order = {i: i for i in range(molecule.n_atoms)}
            shuffled_atom_order.update(
                {
                    old_index: new_index
                    for old_index, new_index in zip(atoms_to_shuffle, atom_permutation)
                },
            )

            molecule = molecule.remap(shuffled_atom_order)

            output_coordinates = self.generate_v_site_coordinates(
                molecule,
                molecule.conformers[0],
                parameter,
            )

            assert output_coordinates.shape == expected_coordinates.shape

            def sort_coordinates(x: numpy.ndarray) -> numpy.ndarray:
                # Sort the rows by first, then second, then third columns as row
                # as the order of v-sites is not deterministic.
                return x[numpy.lexsort((x[:, 2], x[:, 1], x[:, 0])), :]

            found = sort_coordinates(
                output_coordinates.value_in_unit(openmm_unit.angstrom),
            )
            expected = sort_coordinates(expected_coordinates.m_as(unit.angstrom))

            assert numpy.allclose(found, expected), expected - found

    _E = unit.elementary_charge
    _A = unit.angstrom
    _KJ = unit.kilojoule_per_mole

    @pytest.mark.skip(
        "Interchange does not allow this combination of non-bonded settings.",
    )
    @pytest.mark.parametrize(
        ("topology", "parameters", "expected_parameters", "expected_n_v_sites"),
        [
            (
                Topology.from_molecules(
                    [
                        VirtualSiteMocking.chloromethane(reverse=False),
                        VirtualSiteMocking.chloromethane(reverse=True),
                    ],
                ),
                [VirtualSiteMocking.bond_charge_parameter("[Cl:1][C:2]")],
                [
                    # charge, sigma, epsilon
                    (0.1 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.2 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-0.3 * _E, 4.0 * _A, 3.0 * _KJ),
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.2 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.1 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-0.3 * _E, 4.0 * _A, 3.0 * _KJ),
                ],
                2,
            ),
            # Check a complex case of bond charge vsite assignment and overwriting
            (
                Topology.from_molecules(
                    [
                        VirtualSiteMocking.formaldehyde(reverse=False),
                        VirtualSiteMocking.formaldehyde(reverse=True),
                    ],
                ),
                [
                    VirtualSiteMocking.bond_charge_parameter("[O:1]=[*:2]"),
                    VirtualSiteMocking.bond_charge_parameter(
                        "[O:1]=[C:2]",
                        param_multiple=1.5,
                    ),
                    VirtualSiteMocking.bond_charge_parameter(
                        "[O:1]=[CX3:2]",
                        param_multiple=2.0,
                        name="EP2",
                    ),
                ],
                [
                    # charge, sigma, epsilon
                    (0.35 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.7 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-0.45 * _E, 6.0 * _A, 4.5 * _KJ),  # C=O vsite
                    (-0.6 * _E, 8.0 * _A, 6.0 * _KJ),  # CX3=O vsite
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.0 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.7 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.35 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-0.45 * _E, 6.0 * _A, 4.5 * _KJ),  # C=O vsite
                    (-0.6 * _E, 8.0 * _A, 6.0 * _KJ),  # CX3=O vsite
                ],
                4,
            ),
            (
                Topology.from_molecules(
                    [
                        VirtualSiteMocking.formaldehyde(reverse=False),
                        VirtualSiteMocking.formaldehyde(reverse=True),
                    ],
                ),
                [VirtualSiteMocking.monovalent_parameter("[O:1]=[C:2]-[H:3]")],
                [
                    # charge, sigma, epsilon
                    (0.2 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.4 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.3 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.3 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-0.6 * _E, 4.0 * _A, 5.0 * _KJ),
                    (-0.6 * _E, 4.0 * _A, 5.0 * _KJ),
                    (0.3 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.3 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.4 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.2 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-0.6 * _E, 4.0 * _A, 5.0 * _KJ),
                    (-0.6 * _E, 4.0 * _A, 5.0 * _KJ),
                ],
                4,
            ),
            (
                Topology.from_molecules(
                    [
                        VirtualSiteMocking.hypochlorous_acid(reverse=False),
                        VirtualSiteMocking.hypochlorous_acid(reverse=True),
                    ],
                ),
                [
                    VirtualSiteMocking.divalent_parameter(
                        "[H:2][O:1][Cl:3]",
                        match="all_permutations",
                    ),
                ],
                [
                    # charge, sigma, epsilon
                    (0.2 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.1 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.3 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-0.6 * _E, 4.0 * _A, 5.0 * _KJ),
                    (0.3 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.1 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.2 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-0.6 * _E, 4.0 * _A, 5.0 * _KJ),
                ],
                2,
            ),
            (
                Topology.from_molecules(
                    [
                        VirtualSiteMocking.fake_ammonia(reverse=False),
                        VirtualSiteMocking.fake_ammonia(reverse=True),
                    ],
                ),
                [VirtualSiteMocking.trivalent_parameter("[N:1]([Br:2])([Cl:3])[H:4]")],
                [
                    # charge, sigma, epsilon
                    (0.1 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.3 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.2 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.4 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-1.0 * _E, 5.0 * _A, 6.0 * _KJ),
                    (0.4 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.2 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.3 * _E, 10.0 * _A, 0.0 * _KJ),
                    (0.1 * _E, 10.0 * _A, 0.0 * _KJ),
                    (-1.0 * _E, 5.0 * _A, 6.0 * _KJ),
                ],
                2,
            ),
        ],
    )
    def test_create_force(
        self,
        topology: Topology,
        parameters: list[VirtualSiteHandler.VirtualSiteType],
        expected_parameters: list[
            tuple[openmm_unit.Quantity, openmm_unit.Quantity, openmm_unit.Quantity]
        ],
        expected_n_v_sites: int,
    ):
        expected_n_total = topology.n_atoms + expected_n_v_sites
        # sanity check the test input
        assert len(expected_parameters) == expected_n_total

        handler = VirtualSiteHandler(version="0.3")

        for parameter in parameters:
            handler.add_parameter(parameter=parameter)

        force_field = ForceField()

        force_field.register_parameter_handler(ElectrostaticsHandler(version=0.4))
        force_field.register_parameter_handler(LibraryChargeHandler(version=0.3))
        force_field.get_parameter_handler("LibraryCharges").add_parameter(
            {
                "smirks": "[*:1]",
                "charge": [0.0] * unit.elementary_charge,
            },
        )
        force_field.register_parameter_handler(vdWHandler(version=0.3))
        force_field.get_parameter_handler("vdW").add_parameter(
            {
                "smirks": "[*:1]",
                "epsilon": 0.0 * unit.kilojoule_per_mole,
                "sigma": 1.0 * unit.nanometer,
            },
        )
        force_field.register_parameter_handler(handler)

        system: openmm.System = Interchange.from_smirnoff(
            force_field,
            topology,
        ).to_openmm(
            combine_nonbonded_forces=True,
        )

        assert system.getNumParticles() == expected_n_total

        # TODO: Explicitly ensure virtual sites are collated between moleucles
        #       This is implicitly tested in the construction of the parameter arrays

        assert system.getNumForces() == 1
        force: openmm.NonbondedForce = next(iter(system.getForces()))

        total_charge = 0.0 * openmm_unit.elementary_charge

        for i, (expected_charge, expected_sigma, expected_epsilon) in enumerate(
            expected_parameters,
        ):
            charge, sigma, epsilon = force.getParticleParameters(i)

            # Make sure v-sites are massless.
            assert (
                numpy.isclose(system.getParticleMass(i)._value, 0.0)
            ) == system.isVirtualSite(i)

            assert numpy.isclose(
                expected_charge.m_as(unit.elementary_charge),
                charge.value_in_unit(openmm_unit.elementary_charge),
            )
            assert numpy.isclose(
                expected_sigma.m_as(unit.angstrom),
                sigma.value_in_unit(openmm_unit.angstrom),
            )
            assert numpy.isclose(
                expected_epsilon.m_as(unit.kilojoule_per_mole),
                epsilon.value_in_unit(openmm_unit.kilojoule_per_mole),
            )

            total_charge += charge

        expected_total_charge = sum(
            molecule.total_charge.m_as(unit.elementary_charge)
            for molecule in topology.molecules
        )

        assert numpy.isclose(
            expected_total_charge,
            total_charge.value_in_unit(openmm_unit.elementary_charge),
        )
