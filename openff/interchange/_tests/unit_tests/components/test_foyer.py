from pathlib import Path

import numpy as np
import parmed as pmd
import pytest
from openff.toolkit.topology import Molecule, Topology
from openff.units import unit
from openff.utilities.testing import has_package, skip_if_missing

from openff.interchange import Interchange
from openff.interchange._tests import (
    HAS_GROMACS,
    _BaseTest,
    get_test_files_dir_path,
    needs_gmx,
)
from openff.interchange.components.potentials import Potential
from openff.interchange.constants import kj_mol
from openff.interchange.drivers import get_openmm_energies
from openff.interchange.models import PotentialKey, TopologyKey

if has_package("foyer"):
    import foyer

    from openff.interchange.foyer._valence import _RBTorsionHandler

if HAS_GROMACS:
    from openff.interchange.drivers.gromacs import (
        _get_mdp_file,
        _process,
        _run_gmx_energy,
        get_gromacs_energies,
    )


@skip_if_missing("foyer")
class TestFoyer(_BaseTest):
    @pytest.fixture(scope="session")
    def oplsaa(self):
        return foyer.forcefields.load_OPLSAA()

    @pytest.fixture(scope="session")
    def oplsaa_interchange_ethanol(self, oplsaa):
        molecule = Molecule.from_file(
            get_test_files_dir_path("foyer_test_molecules") / "ethanol.sdf",
        )
        molecule.name = "ETH"

        top = molecule.to_topology()

        interchange = Interchange.from_foyer(topology=top, force_field=oplsaa)
        interchange.positions = molecule.conformers[0].m_as(unit.nanometer)
        interchange.box = [4, 4, 4]
        return interchange

    @pytest.fixture(scope="session")
    def get_interchanges(self, oplsaa):
        def interchanges_from_path(molecule_path):
            molecule_or_molecules = Molecule.from_file(molecule_path)
            mol_name = Path(molecule_path).stem[0:4].upper()

            if isinstance(molecule_or_molecules, list):
                for idx, mol in enumerate(molecule_or_molecules):
                    mol.name = f"{mol_name}{idx}"
            else:
                molecule_or_molecules.name = mol_name
                molecule_or_molecules = [molecule_or_molecules]

            topology = Topology.from_molecules(molecule_or_molecules)
            openff_interchange = Interchange.from_foyer(oplsaa, topology)

            if isinstance(molecule_or_molecules, list):
                openff_interchange.positions = np.vstack(
                    tuple(
                        molecule.conformers[0].m_as(unit.nanometer)
                        for molecule in molecule_or_molecules
                    ),
                )
            else:
                openff_interchange.positions = molecule_or_molecules.conformers[0].m_as(
                    unit.nanometer,
                )

            openff_interchange.box = [4, 4, 4]

            parmed_struct = pmd.openmm.load_topology(topology.to_openmm())
            parmed_struct.positions = openff_interchange.positions.m_as(unit.angstrom)
            parmed_struct.box = [40, 40, 40, 90, 90, 90]

            return openff_interchange, parmed_struct

        return interchanges_from_path

    def test_handlers_exist(self, oplsaa_interchange_ethanol):
        for _, handler in oplsaa_interchange_ethanol.collections.items():
            assert handler

        assert oplsaa_interchange_ethanol["vdW"].scale_14 == 0.5
        assert oplsaa_interchange_ethanol["Electrostatics"].scale_14 == 0.5

    @needs_gmx
    @pytest.mark.slow()
    @pytest.mark.skip(
        reason="Exporting to OpenMM with geometric mixing rules not yet implemented",
    )
    def test_ethanol_energies(self, oplsaa_interchange_ethanol):
        from openff.interchange.drivers import get_gromacs_energies

        gmx_energies = get_gromacs_energies(oplsaa_interchange_ethanol)
        omm_energies = get_openmm_energies(oplsaa_interchange_ethanol)

        gmx_energies.compare(
            omm_energies,
            {
                "vdW": 12.0 * unit.kilojoule / unit.mole,
                "Electrostatics": 12.0 * unit.kilojoule / unit.mole,
            },
        )

    @needs_gmx
    @pytest.mark.parametrize(
        argnames="molecule_path",
        argvalues=[*get_test_files_dir_path("foyer_test_molecules").glob("*.sdf")],
    )
    @pytest.mark.slow()
    def test_interchange_energies(self, molecule_path, get_interchanges, oplsaa):
        if "ethanol" in str(molecule_path) or "adamantane" in str(molecule_path):
            pytest.skip("Foyer/ParmEd bug with this molecule")
        openff_interchange, pmd_structure = get_interchanges(molecule_path)
        parameterized_pmd_structure = oplsaa.apply(pmd_structure)

        openff_energy = get_gromacs_energies(
            openff_interchange,
            mdp="cutoff_hbonds",
            round_positions=3,
        )

        parameterized_pmd_structure.save("from_foyer.gro")
        parameterized_pmd_structure.save("from_foyer.top")

        through_foyer = _process(
            _run_gmx_energy(
                top_file="from_foyer.top",
                gro_file="from_foyer.gro",
                mdp_file=_get_mdp_file("cutoff_hbonds"),
            ),
        )

        # TODO: Revisit after https://github.com/mosdef-hub/foyer/issues/431
        openff_energy.compare(
            through_foyer,
            {
                "Bond": 1e-3 * unit.kilojoule / unit.mole,
                "Angle": 1e-3 * unit.kilojoule / unit.mole,
                "Torsion": 1e-3 * unit.kilojoule / unit.mole,
                "vdW": 5 * unit.kilojoule / unit.mole,
                "Electrostatics": 1 * unit.kilojoule / unit.mole,
            },
        )


@skip_if_missing("foyer")
class TestRBTorsions(TestFoyer):
    @pytest.fixture()
    def ethanol_with_rb_torsions(self, sage):
        mol = Molecule.from_smiles("CC")
        mol.name = "ETH"
        mol.generate_conformers(n_conformers=1)
        top = mol.to_topology()

        out = Interchange.from_smirnoff(sage, top)
        out.box = [4, 4, 4]
        out.positions = mol.conformers[0]
        out.positions = np.round(out.positions, 2)

        rb_torsions = _RBTorsionHandler()
        smirks = "[#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]"
        pot_key = PotentialKey(id=smirks)
        for proper in top.propers:
            top_key = TopologyKey(atom_indices=tuple(top.atom_index(a) for a in proper))
            rb_torsions.key_map.update({top_key: pot_key})

        # Values from HC-CT-CT-HC RB torsion
        # https://github.com/mosdef-hub/foyer/blob/7816bf53a127502520a18d76c81510f96adfdbed/foyer/forcefields/xml/oplsaa.xml#L2585
        pot = Potential(
            parameters={
                "c0": 0.6276 * kj_mol,
                "c1": 1.8828 * kj_mol,
                "c2": 0.0 * kj_mol,
                "c3": -2.5104 * kj_mol,
                "c4": 0.0 * kj_mol,
                "c5": 0.0 * kj_mol,
            },
        )

        rb_torsions.potentials.update({pot_key: pot})

        out.collections.update({"RBTorsions": rb_torsions})
        out.collections.pop("ProperTorsions")

        return out

    @needs_gmx
    @pytest.mark.slow()
    @pytest.mark.skip(reason="Something is broken with RBTorsions in OpenMM export")
    def test_rb_torsions(self, ethanol_with_rb_torsions):
        omm = get_openmm_energies(ethanol_with_rb_torsions, round_positions=3).energies[
            "Torsion"
        ]
        gmx = get_gromacs_energies(
            ethanol_with_rb_torsions,
            round_positions=3,
        ).energies["Torsion"]

        assert (gmx - omm).m_as(kj_mol) < 1e-6

    @pytest.mark.slow()
    @skip_if_missing("foyer")
    @skip_if_missing("mbuild")
    @needs_gmx
    def test_rb_torsions_vs_foyer(self, oplsaa, ethanol_with_rb_torsions):
        # Given that these force constants are copied from Foyer's OPLS-AA file,
        # compare to processing through the current MoSDeF pipeline
        import mbuild

        comp = mbuild.load("CC", smiles=True)
        comp.xyz = ethanol_with_rb_torsions.positions.m_as(unit.nanometer)

        from_foyer = oplsaa.apply(comp)
        from_foyer.box = [40, 40, 40, 90, 90, 90]
        from_foyer.save("from_foyer.top")
        from_foyer.save("from_foyer.gro")

        rb_torsion_energy_from_foyer = _process(
            _run_gmx_energy(
                top_file="from_foyer.top",
                gro_file="from_foyer.gro",
                mdp_file=_get_mdp_file("default"),
            ),
        ).energies["Torsion"]

        # GROMACS vs. OpenMM was already compared, so just use one
        omm = get_gromacs_energies(
            ethanol_with_rb_torsions,
            round_positions=3,
        ).energies["Torsion"]

        assert (omm - rb_torsion_energy_from_foyer).m_as(kj_mol) < 1e-6
