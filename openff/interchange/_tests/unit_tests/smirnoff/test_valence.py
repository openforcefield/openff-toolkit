import numpy
import openmm
import pytest
from openff.toolkit import ForceField, Molecule, Topology
from openff.toolkit.tests.test_forcefield import (
    create_cyclohexane,
    create_ethanol,
    create_reversed_ethanol,
    create_water,
)
from openff.toolkit.tests.utils import requires_openeye
from openff.toolkit.typing.engines.smirnoff.parameters import (
    AngleHandler,
    BondHandler,
    ImproperTorsionHandler,
)
from openff.units import unit
from pydantic import ValidationError

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest
from openff.interchange.constants import kcal_mol_a2, kcal_mol_rad2
from openff.interchange.exceptions import DuplicateMoleculeError
from openff.interchange.models import AngleKey, BondKey, ImproperTorsionKey
from openff.interchange.smirnoff._create import _create_interchange
from openff.interchange.smirnoff._valence import (
    SMIRNOFFAngleCollection,
    SMIRNOFFBondCollection,
    SMIRNOFFConstraintCollection,
    SMIRNOFFImproperTorsionCollection,
    _check_molecule_uniqueness,
)


class TestSMIRNOFFValenceCollections(_BaseTest):
    def test_bond_collection(self):
        bond_handler = BondHandler(version=0.3)
        bond_handler.fractional_bondorder_method = "AM1-Wiberg"
        bond_parameter = BondHandler.BondType(
            smirks="[*:1]~[*:2]",
            k=1.5 * unit.kilocalorie_per_mole / unit.angstrom**2,
            length=1.5 * unit.angstrom,
            id="b1000",
        )
        bond_handler.add_parameter(bond_parameter.to_dict())

        forcefield = ForceField()
        forcefield.register_parameter_handler(bond_handler)
        bond_potentials = SMIRNOFFBondCollection.create(
            parameter_handler=forcefield["Bonds"],
            topology=Molecule.from_smiles("O").to_topology(),
        )

        top_key = BondKey(atom_indices=(0, 1))
        pot_key = bond_potentials.key_map[top_key]
        assert pot_key.associated_handler == "Bonds"
        pot = bond_potentials.potentials[pot_key]

        assert pot.parameters["k"].to(kcal_mol_a2).magnitude == pytest.approx(1.5)

    def test_angle_collection(self):
        angle_handler = AngleHandler(version=0.3)
        angle_parameter = AngleHandler.AngleType(
            smirks="[*:1]~[*:2]~[*:3]",
            k=2.5 * unit.kilocalorie_per_mole / unit.radian**2,
            angle=100 * unit.degree,
            id="b1000",
        )
        angle_handler.add_parameter(angle_parameter.to_dict())

        forcefield = ForceField()
        forcefield.register_parameter_handler(angle_handler)
        angle_potentials = SMIRNOFFAngleCollection.create(
            parameter_handler=forcefield["Angles"],
            topology=Molecule.from_smiles("CCC").to_topology(),
        )

        top_key = AngleKey(atom_indices=(0, 1, 2))
        pot_key = angle_potentials.key_map[top_key]
        assert pot_key.associated_handler == "Angles"
        pot = angle_potentials.potentials[pot_key]

        assert pot.parameters["k"].to(kcal_mol_rad2).magnitude == pytest.approx(2.5)

    def test_store_improper_torsion_matches(self, formaldehyde):
        parameter_handler = ImproperTorsionHandler(version=0.3)
        parameter_handler.add_parameter(
            parameter=ImproperTorsionHandler.ImproperTorsionType(
                smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]",
                periodicity1=2,
                phase1=180.0 * unit.degree,
                k1=1.1 * unit.kilocalorie_per_mole,
            ),
        )

        collection = SMIRNOFFImproperTorsionCollection()
        collection.store_matches(parameter_handler, formaldehyde.to_topology())

        assert len(collection.key_map) == 3

        for indices in [
            (0, 1, 2, 3),
            (0, 2, 3, 1),
            (0, 3, 1, 2),
        ]:
            key = ImproperTorsionKey(atom_indices=indices, mult=0)
            assert key in collection.key_map

    def test_store_nonstandard_improper_idivf(self, acetaldehyde):
        handler = ImproperTorsionHandler(version=0.3)
        handler.add_parameter(
            {
                "smirks": "[*:1]~[#6:2](~[#8:3])~[*:4]",
                "periodicity1": 2,
                "phase1": 180.0 * unit.degree,
                "k1": 1.1 * unit.kilocalorie_per_mole,
                "idivf1": 1.234 * unit.dimensionless,
                "id": "i1",
            },
        )

        collection = SMIRNOFFImproperTorsionCollection.create(
            parameter_handler=handler,
            topology=acetaldehyde.to_topology(),
        )

        handler = ImproperTorsionHandler(version=0.3)
        handler.default_idivf = 5.555
        handler.add_parameter(
            {
                "smirks": "[*:1]~[#6:2](~[#8:3])~[*:4]",
                "periodicity1": 2,
                "phase1": 180.0 * unit.degree,
                "k1": 1.1 * unit.kilocalorie_per_mole,
                "id": "i1",
            },
        )

        collection = SMIRNOFFImproperTorsionCollection.create(
            parameter_handler=handler,
            topology=acetaldehyde.to_topology(),
        )

        assert [*collection.potentials.values()][0].parameters[
            "idivf"
        ] == 5.555 * unit.dimensionless


class TestBondCollection(_BaseTest):
    def test_upconvert_warning(self, parsley, ethanol):
        from packaging.version import Version

        assert parsley["Bonds"].version == Version("0.3")
        with pytest.warns(
            UserWarning,
            match="Automatically up-converting BondHandler from version 0.3 to 0.4.",
        ):
            _create_interchange(parsley, [ethanol])


class TestConstraintCollection(_BaseTest):
    @pytest.mark.parametrize(
        ("mol", "n_constraints"),
        [
            ("C", 4),
            ("CC", 6),
        ],
    )
    def test_num_constraints(self, sage, mol, n_constraints):
        bond_handler = sage["Bonds"]
        constraint_handler = sage["Constraints"]

        topology = Molecule.from_smiles(mol).to_topology()

        constraints = SMIRNOFFConstraintCollection.create(
            parameter_handler=[bond_handler, constraint_handler],
            topology=topology,
            bonds=SMIRNOFFBondCollection.create(bond_handler, topology),
        )

        assert len(constraints.key_map) == n_constraints

    def test_constraints_with_distance(self, tip3p, water):
        topology = water.to_topology()
        topology.box_vectors = [4, 4, 4] * unit.nanometer

        constraints = SMIRNOFFConstraintCollection.create(
            parameter_handler=tip3p["Constraints"],
            topology=topology,
        )

        assert len(constraints.key_map) == 3
        assert len(constraints.potentials) == 2


class TestBondOrderInterpolation(_BaseTest):
    @pytest.mark.slow()
    def test_input_bond_orders_ignored(self):
        """Test that conformers existing in the topology are not considered in the bond order interpolation
        part of the parametrization process"""

        mol = create_ethanol()
        mol.assign_fractional_bond_orders(bond_order_model="am1-wiberg")
        mod_mol = Molecule(mol)
        for bond in mod_mol.bonds:
            bond.fractional_bond_order += 0.1

        top = Topology.from_molecules(mol)
        mod_top = Topology.from_molecules(mod_mol)

        forcefield = ForceField(
            "openff-2.0.0.offxml",
            self.xml_ff_bo_bonds,
        )

        bonds = SMIRNOFFBondCollection.create(
            parameter_handler=forcefield["Bonds"],
            topology=top,
        )
        bonds_mod = SMIRNOFFBondCollection.create(
            parameter_handler=forcefield["Bonds"],
            topology=mod_top,
        )

        for pot_key1, pot_key2 in zip(
            bonds.key_map.values(),
            bonds_mod.key_map.values(),
        ):
            k1 = bonds.potentials[pot_key1].parameters["k"].m_as(kcal_mol_a2)
            k2 = bonds_mod.potentials[pot_key2].parameters["k"].m_as(kcal_mol_a2)
            assert k1 == pytest.approx(k2, rel=1e-5), (k1, k2)

    def test_input_conformers_ignored(self):
        """Test that conformers existing in the topology are not considered in the bond order interpolation
        part of the parametrization process"""
        from openff.toolkit.tests.test_forcefield import create_ethanol

        mol = create_ethanol()
        mol.assign_fractional_bond_orders(bond_order_model="am1-wiberg")
        mod_mol = Molecule(mol)
        mod_mol.generate_conformers()
        tmp = mod_mol._conformers[0][0][0]
        mod_mol._conformers[0][0][0] = mod_mol._conformers[0][1][0]
        mod_mol._conformers[0][1][0] = tmp

        top = Topology.from_molecules(mol)
        mod_top = Topology.from_molecules(mod_mol)

        forcefield = ForceField(
            "openff-2.0.0.offxml",
            self.xml_ff_bo_bonds,
        )

        bonds = SMIRNOFFBondCollection.create(
            parameter_handler=forcefield["Bonds"],
            topology=top,
        )
        bonds_mod = SMIRNOFFBondCollection.create(
            parameter_handler=forcefield["Bonds"],
            topology=mod_top,
        )

        for key1, key2 in zip(bonds.potentials, bonds_mod.potentials):
            k1 = bonds.potentials[key1].parameters["k"].m_as(kcal_mol_a2)
            k2 = bonds_mod.potentials[key2].parameters["k"].m_as(kcal_mol_a2)
            assert k1 == pytest.approx(k2, rel=1e-5), (k1, k2)

    def test_fractional_bondorder_invalid_interpolation_method(self, ethanol):
        """
        Ensure that requesting an invalid interpolation method leads to a
        FractionalBondOrderInterpolationMethodUnsupportedError
        """
        forcefield = ForceField(
            "openff-2.0.0.offxml",
            self.xml_ff_bo_bonds,
        )
        forcefield["Bonds"]._fractional_bondorder_interpolation = "invalid method name"

        # TODO: Make this a more descriptive custom exception
        with pytest.raises(ValidationError):
            Interchange.from_smirnoff(forcefield, [ethanol])


class TestParameterInterpolation(_BaseTest):
    xml_ff_bo = """<?xml version='1.0' encoding='ASCII'?>
    <SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
      <Bonds version="0.3" fractional_bondorder_method="AM1-Wiberg"
        fractional_bondorder_interpolation="linear">
        <Bond
          smirks="[#6X4:1]~[#8X2:2]"
          id="bbo1"
          k_bondorder1="101.0 * kilocalories_per_mole/angstrom**2"
          k_bondorder2="123.0 * kilocalories_per_mole/angstrom**2"
          length_bondorder1="1.4 * angstrom"
          length_bondorder2="1.3 * angstrom"
          />
      </Bonds>
      <ProperTorsions version="0.3" potential="k*(1+cos(periodicity*theta-phase))">
        <Proper smirks="[*:1]~[#6X3:2]~[#6X3:3]~[*:4]" id="tbo1" periodicity1="2" phase1="0.0 * degree"
        k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
        <Proper smirks="[*:1]~[#6X4:2]~[#8X2:3]~[*:4]" id="tbo2" periodicity1="2" phase1="0.0 * degree"
        k1_bondorder1="1.00*kilocalories_per_mole" k1_bondorder2="1.80*kilocalories_per_mole" idivf1="1.0"/>
      </ProperTorsions>
    </SMIRNOFF>
    """

    @pytest.mark.xfail(reason="Not yet implemented using input bond orders")
    def test_bond_order_interpolation(self, ethanol):
        forcefield = ForceField(
            "openff-2.0.0.offxml",
            self.xml_ff_bo,
        )

        ethanol.generate_conformers(n_conformers=1)

        ethanol.bonds[1].fractional_bond_order = 1.5

        top = ethanol.to_topology()

        out = Interchange.from_smirnoff(forcefield, ethanol.to_topology())

        top_key = BondKey(
            atom_indices=(1, 2),
            bond_order=top.get_bond_between(1, 2).bond.fractional_bond_order,
        )

        found_k = out["Bonds"].potentials[out["Bonds"].key_map[top_key]].parameters["k"]
        assert found_k == 300 * kcal_mol_a2

    @pytest.mark.slow()
    @pytest.mark.xfail(reason="Not yet implemented using input bond orders")
    def test_bond_order_interpolation_similar_bonds(self):
        """Test that key mappings do not get confused when two bonds having similar SMIRKS matches
        have different bond orders"""
        forcefield = ForceField(
            "openff-2.0.0.offxml",
            self.xml_ff_bo,
        )

        # TODO: Construct manually to avoid relying on atom ordering
        mol = Molecule.from_smiles("C(CCO)O")
        mol.generate_conformers(n_conformers=1)

        mol.bonds[2].fractional_bond_order = 1.5
        mol.bonds[3].fractional_bond_order = 1.2

        top = mol.to_topology()

        out = Interchange.from_smirnoff(forcefield, top)

        bond1_top_key = BondKey(
            atom_indices=(2, 3),
            bond_order=top.get_bond_between(2, 3).bond.fractional_bond_order,
        )
        bond1_pot_key = out["Bonds"].key_map[bond1_top_key]

        bond2_top_key = BondKey(
            atom_indices=(0, 4),
            bond_order=top.get_bond_between(0, 4).bond.fractional_bond_order,
        )
        bond2_pot_key = out["Bonds"].key_map[bond2_top_key]

        assert numpy.allclose(
            out["Bonds"].potentials[bond1_pot_key].parameters["k"],
            300.0 * unit.Unit("kilocalories / mol / angstrom ** 2"),
        )

        assert numpy.allclose(
            out["Bonds"].potentials[bond2_pot_key].parameters["k"],
            180.0 * unit.Unit("kilocalories / mol / angstrom ** 2"),
        )

    @requires_openeye
    @pytest.mark.parametrize(
        (
            "get_molecule",
            "k_torsion_interpolated",
            "k_bond_interpolated",
            "length_bond_interpolated",
            "central_atoms",
        ),
        [
            (create_ethanol, 4.16586914, 42208.5402, 0.140054167256, (1, 2)),
            (create_reversed_ethanol, 4.16564555, 42207.9252, 0.14005483525, (7, 6)),
        ],
    )
    def test_fractional_bondorder_from_molecule(
        self,
        get_molecule,
        k_torsion_interpolated,
        k_bond_interpolated,
        length_bond_interpolated,
        central_atoms,
    ):
        """Copied from the toolkit with modified reference constants.
        Force constant computed by interpolating (k1, k2) = (101, 123) kcal/A**2/mol
        with bond order 1.00093035 (AmberTools 21.4, Python 3.8, macOS):
            101 + (123 - 101) * (0.00093035) = 101.0204677 kcal/A**2/mol
            = 42266.9637 kJ/nm**2/mol

        Same process with bond length (1.4, 1.3) A gives 0.1399906965 nm
        Same process with torsion k (1.0, 1.8) kcal/mol gives 4.18711406752 kJ/mol

        Using OpenEye (openeye-toolkits 2021.1.1, Python 3.8, macOS):
            bond order 0.9945832743790813
            bond k = 42208.5402 kJ/nm**2/mol
            bond length = 0.14005416725620918 nm
            torsion k = 4.16586914 kilojoules kJ/mol

        ... except OpenEye has a different fractional bond order for reversed ethanol
            bond order 0.9945164749654242
            bond k = 42207.9252 kJ/nm**2/mol
            bond length = 0.14005483525034576 nm
            torsion k = 4.16564555 kJ/mol

        """
        mol = get_molecule()
        forcefield = ForceField(
            "openff-2.0.0.offxml",
            self.xml_ff_bo,
        )
        topology = Topology.from_molecules(mol)

        out = Interchange.from_smirnoff(forcefield, topology)
        out.box = unit.Quantity(4 * numpy.eye(3), unit.nanometer)
        omm_system = out.to_openmm(combine_nonbonded_forces=True)

        # Verify that the assigned bond parameters were correctly interpolated
        off_bond_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.HarmonicBondForce)
        ][0]

        for idx in range(off_bond_force.getNumBonds()):
            params = off_bond_force.getBondParameters(idx)

            atom1, atom2 = params[0], params[1]
            atom1_mol, atom2_mol = central_atoms

            if ((atom1 == atom1_mol) and (atom2 == atom2_mol)) or (
                (atom1 == atom2_mol) and (atom2 == atom1_mol)
            ):
                k = params[-1]
                length = params[-2]
                numpy.testing.assert_allclose(
                    k / k.unit,
                    k_bond_interpolated,
                    atol=0,
                    rtol=2e-6,
                )
                numpy.testing.assert_allclose(
                    length / length.unit,
                    length_bond_interpolated,
                    atol=0,
                    rtol=2e-6,
                )

        # Verify that the assigned torsion parameters were correctly interpolated
        off_torsion_force = [
            force
            for force in omm_system.getForces()
            if isinstance(force, openmm.PeriodicTorsionForce)
        ][0]

        for idx in range(off_torsion_force.getNumTorsions()):
            params = off_torsion_force.getTorsionParameters(idx)

            atom2, atom3 = params[1], params[2]
            atom2_mol, atom3_mol = central_atoms

            if ((atom2 == atom2_mol) and (atom3 == atom3_mol)) or (
                (atom2 == atom3_mol) and (atom3 == atom2_mol)
            ):
                k = params[-1]
                numpy.testing.assert_allclose(
                    k / k.unit,
                    k_torsion_interpolated,
                    atol=0,
                    rtol=2e-6,
                )


def test_get_uses_interpolation():
    handler_no_interpolation = BondHandler(version=0.4)
    handler_all_interpolation = BondHandler(version=0.4)
    handler_some_interpolation = BondHandler(version=0.4)
    handler_partial_interpolation = BondHandler(version=0.4)

    for handler in [
        handler_no_interpolation,
        handler_some_interpolation,
    ]:
        handler.add_parameter(
            {
                "smirks": "[#1:1]-[*:2]",
                "k": 400 * kcal_mol_a2,
                "length": 1 * unit.angstrom,
            },
        )

    for handler in [
        handler_all_interpolation,
        handler_some_interpolation,
    ]:
        handler.add_parameter(
            {
                "smirks": "[#6:1]-[#6:2]",
                "k_bondorder1": 400 * kcal_mol_a2,
                "length_bondorder1": 1 * unit.angstrom,
                "k_bondorder2": 500 * kcal_mol_a2,
                "length_bondorder2": 1.2 * unit.angstrom,
            },
        )

    handler_partial_interpolation.add_parameter(
        {
            "smirks": "[#6:1]-[#7:2]",
            "k": 400 * kcal_mol_a2,
            "length_bondorder1": 1 * unit.angstrom,
            "length_bondorder2": 0.5 * unit.angstrom,
        },
    )

    assert not SMIRNOFFBondCollection()._get_uses_interpolation(
        handler_no_interpolation,
    )
    assert SMIRNOFFBondCollection()._get_uses_interpolation(handler_some_interpolation)
    assert SMIRNOFFBondCollection()._get_uses_interpolation(handler_all_interpolation)
    assert SMIRNOFFBondCollection()._get_uses_interpolation(
        handler_partial_interpolation,
    )


def test_check_molecule_uniqueness():
    _check_molecule_uniqueness(None)
    _check_molecule_uniqueness(list())
    _check_molecule_uniqueness([create_ethanol()])
    _check_molecule_uniqueness([create_cyclohexane()])
    _check_molecule_uniqueness([create_ethanol(), create_cyclohexane(), create_water()])

    with pytest.raises(DuplicateMoleculeError, match="Duplicate molecules"):
        _check_molecule_uniqueness(2 * [create_ethanol()])

    with pytest.raises(DuplicateMoleculeError, match="Duplicate molecules"):
        _check_molecule_uniqueness([create_ethanol(), create_reversed_ethanol()])
