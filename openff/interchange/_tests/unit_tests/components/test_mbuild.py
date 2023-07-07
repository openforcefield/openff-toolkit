import numpy as np
from openff.toolkit.topology import Molecule, Topology
from openff.units import unit
from openff.utilities.testing import skip_if_missing

from openff.interchange.components.mbuild import offmol_to_compound, offtop_to_compound


@skip_if_missing("mbuild")
class TestMBuildConversions:
    def test_basic_mol_to_compound(self):
        """Test basic behavior of conversion to mBuild Compound"""
        offmol = Molecule.from_smiles("CCO")
        offmol.generate_conformers(n_conformers=1)

        comp = offmol_to_compound(offmol)

        assert comp.n_particles == offmol.n_atoms
        assert comp.n_bonds == offmol.n_bonds

        np.testing.assert_equal(
            offmol.conformers[0].m_as(unit.nanometer),
            comp.xyz,
        )

    def test_mbuild_conversion_generate_conformers(self):
        """Test that a single conformer is automatically generated"""
        offmol = Molecule.from_smiles("CCO")

        comp = offmol_to_compound(offmol)

        assert comp.n_particles == offmol.n_atoms
        assert comp.n_bonds == offmol.n_bonds

        offmol.generate_conformers(n_conformers=1)
        expected_conf = offmol.conformers[0]

        np.testing.assert_equal(
            expected_conf.m_as(unit.nanometer),
            comp.xyz,
        )

    def test_mbuild_conversion_first_conformer_used(self):
        """Test that only the first conformer in an OFFMol is used"""
        offmol = Molecule.from_smiles("C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=CC=C3")
        offmol.generate_conformers(n_conformers=3, rms_cutoff=0.0 * unit.angstrom)

        comp = offmol_to_compound(offmol)

        np.testing.assert_equal(
            offmol.conformers[0].m_as(unit.nanometer),
            comp.xyz,
        )

        with np.testing.assert_raises(AssertionError):
            np.testing.assert_equal(
                offmol.conformers[1].m_as(unit.nanometer),
                comp.xyz,
            )

        with np.testing.assert_raises(AssertionError):
            np.testing.assert_equal(
                offmol.conformers[2].m_as(unit.nanometer),
                comp.xyz,
            )

    def test_mbuild_conversion_element_names(self):
        """Test that the generated Compound has particle names that can be
        interpreted as elements"""
        offmol = Molecule.from_smiles(
            "CSC(CO[P]([O-])([O-])=O)c1cc(Br)c(F)c(Cl)n1",
            allow_undefined_stereo=True,
        )
        comp = offmol_to_compound(offmol)

        known_elements = {"C", "H", "O", "N", "Cl", "Br", "F", "S", "P"}

        for particle in comp.particles():
            assert particle.name in known_elements

    def test_multi_mol_topology_to_compound(self):
        """Test the basic behavior of a (multi-mol) OFFTop"""
        ethanol = Molecule.from_smiles("CCO")
        ethanol.name = "ETH"
        methane = Molecule.from_smiles("C")
        methane.name = "MET"

        top = Topology.from_molecules([ethanol, ethanol, methane, methane])

        comp = offtop_to_compound(top)

        assert comp.n_particles == 28  # 2 * (9 + 5)
        assert comp.n_bonds == 24  # 2 * (8 + 4)
        assert len(comp.children) == 4  # 4 "molecules"
