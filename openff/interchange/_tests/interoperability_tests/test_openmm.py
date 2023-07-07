import math

import numpy
import openmm
import openmm.app
import openmm.unit
import pytest
from openff.toolkit.tests.test_forcefield import create_ethanol
from openff.toolkit.tests.utils import get_data_file_path
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField, VirtualSiteHandler
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest, get_test_file_path
from openff.interchange._tests.unit_tests.plugins.test_smirnoff_plugins import (
    TestDoubleExponential,
)
from openff.interchange.drivers.openmm import get_openmm_energies
from openff.interchange.exceptions import (
    MissingPositionsError,
    PluginCompatibilityError,
    UnsupportedCutoffMethodError,
    UnsupportedExportError,
)
from openff.interchange.interop.openmm import (
    from_openmm,
    to_openmm_positions,
    to_openmm_topology,
)

# WISHLIST: Add tests for reaction-field if implemented

nonbonded_methods = [
    {
        "vdw_method": "cutoff",
        "electrostatics_periodic": "PME",
        "periodic": True,
        "result": openmm.NonbondedForce.PME,
    },
    {
        "vdw_method": "cutoff",
        "electrostatics_periodic": "PME",
        "periodic": False,
        "result": openmm.NonbondedForce.NoCutoff,
    },
    {
        "vdw_method": "PME",
        "electrostatics_periodic": "PME",
        "periodic": True,
        "result": openmm.NonbondedForce.LJPME,
    },
    {
        "vdw_method": "PME",
        "electrostatics_periodic": "PME",
        "periodic": False,
        "result": UnsupportedCutoffMethodError,
    },
]


def _get_num_virtual_sites(openmm_topology: openmm.app.Topology) -> int:
    return sum(atom.element is None for atom in openmm_topology.atoms())


def _compare_openmm_topologies(top1: openmm.app.Topology, top2: openmm.app.Topology):
    """
    In lieu of first-class serializaiton in OpenMM (https://github.com/openmm/openmm/issues/1543),
    do some quick heuristics to roughly compare two OpenMM Topology objects.
    """
    for method_name in [
        "getNumAtoms",
        "getNumBonds",
        "getNumChains",
        "getNumResidues",
    ]:
        assert getattr(top1, method_name)() == getattr(top2, method_name)()

    assert (top1.getPeriodicBoxVectors() == top2.getPeriodicBoxVectors()).all()


class TestOpenMM(_BaseTest):
    @pytest.mark.parametrize("inputs", nonbonded_methods)
    def test_openmm_nonbonded_methods(self, inputs, sage):
        """See test_nonbonded_method_resolution in openff/toolkit/tests/test_forcefield.py"""
        vdw_method = inputs["vdw_method"]
        electrostatics_method = inputs["electrostatics_periodic"]
        periodic = inputs["periodic"]
        result = inputs["result"]

        molecules = [create_ethanol()]

        pdbfile = openmm.app.PDBFile(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"),
        )
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        if not periodic:
            topology.box_vectors = None

        sage.get_parameter_handler("vdW", {}).method = vdw_method
        sage.get_parameter_handler(
            "Electrostatics",
            {},
        ).periodic_potential = electrostatics_method
        interchange = Interchange.from_smirnoff(
            force_field=sage,
            topology=topology,
        )
        if type(result) == int:
            nonbonded_method = result
            # The method is validated and may raise an exception if it's not supported.
            sage.get_parameter_handler("vdW", {}).method = vdw_method
            sage.get_parameter_handler(
                "Electrostatics",
                {},
            ).periodic_potential = electrostatics_method
            interchange = Interchange.from_smirnoff(
                force_field=sage,
                topology=topology,
            )
            openmm_system = interchange.to_openmm(combine_nonbonded_forces=True)
            for force in openmm_system.getForces():
                if isinstance(force, openmm.NonbondedForce):
                    assert force.getNonbondedMethod() == nonbonded_method
                    break
            else:
                raise Exception
        elif issubclass(result, (BaseException, Exception)):
            exception = result
            with pytest.raises(exception):
                interchange.to_openmm(combine_nonbonded_forces=True)
        else:
            raise Exception

    @pytest.mark.skip(reason="Re-implement when SMIRNOFF supports more mixing rules")
    def test_unsupported_mixing_rule(self, sage):
        molecules = [create_ethanol()]
        pdbfile = openmm.app.PDBFile(
            get_data_file_path("systems/test_systems/1_ethanol.pdb"),
        )
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)

        interchange = Interchange.from_smirnoff(
            force_field=sage,
            topology=topology,
        )

        interchange["vdW"].mixing_rule = "geometric"

        with pytest.raises(UnsupportedExportError, match="default NonbondedForce"):
            interchange.to_openmm(combine_nonbonded_forces=True)

    @pytest.mark.xfail(reason="Broken because of splitting non-bonded forces")
    @pytest.mark.slow()
    @pytest.mark.parametrize("mol_smi", ["C", "CC", "CCO"])
    def test_openmm_roundtrip(self, sage, mol_smi):
        mol = Molecule.from_smiles(mol_smi)
        mol.generate_conformers(n_conformers=1)
        top = mol.to_topology()

        interchange = Interchange.from_smirnoff(sage, top)

        interchange.box = [4, 4, 4]
        interchange.positions = mol.conformers[0].value_in_unit(openmm.unit.nanometer)

        converted = from_openmm(
            topology=interchange.to_openmm_topology(),
            system=interchange.to_openmm(combine_nonbonded_forces=True),
        )

        converted.box = interchange.box
        converted.positions = interchange.positions

        get_openmm_energies(interchange).compare(
            get_openmm_energies(converted, combine_nonbonded_forces=True),
        )

    @pytest.mark.xfail(reason="Broken because of splitting non-bonded forces")
    @pytest.mark.slow()
    def test_combine_nonbonded_forces(self, sage):
        mol = Molecule.from_smiles("ClC#CCl")
        mol.name = "HPER"
        mol.generate_conformers(n_conformers=1)

        out = Interchange.from_smirnoff(force_field=sage, topology=mol.to_topology())
        out.box = [4, 4, 4]
        out.positions = mol.conformers[0]

        num_forces_combined = out.to_openmm(
            combine_nonbonded_forces=True,
        ).getNumForces()
        num_forces_uncombined = out.to_openmm(
            combine_nonbonded_forces=False,
        ).getNumForces()

        # The "new" forces are the split-off vdW forces, the 1-4 vdW, and the 1-4 electrostatics
        assert num_forces_combined + 3 == num_forces_uncombined

        separate = get_openmm_energies(out, combine_nonbonded_forces=False)
        combined = get_openmm_energies(out, combine_nonbonded_forces=True)

        assert (
            separate["vdW"] + separate["Electrostatics"] - combined["Nonbonded"]
        ).m < 0.001

    def test_openmm_no_angle_force_if_constrained(self):
        # Sage includes angle parameters for water and also TIP3P constraints
        tip3p = ForceField("openff-2.0.0.offxml")

        topology = Molecule.from_smiles("O").to_topology()
        topology.box_vectors = [4, 4, 4] * unit.nanometer

        interchange = Interchange.from_smirnoff(tip3p, topology)
        openmm_system = interchange.to_openmm(combine_nonbonded_forces=True)

        # The only angle in the system (H-O-H) includes bonds with constrained lengths
        # and a constrained angle, so by convention a force should NOT be added
        for force in openmm_system.getForces():
            if type(force) == openmm.HarmonicAngleForce:
                assert force.getNumAngles() == 0
                break
        else:
            raise Exception("No HarmonicAngleForce found")

    @pytest.mark.skip(reason="Rewrite as a plugin")
    def test_nonharmonic_angle(self, sage, ethanol_top):
        out = Interchange.from_smirnoff(sage, ethanol_top)
        out["Angles"].expression = "k/2*(cos(theta)-cos(angle))**2"

        system = out.to_openmm()

        def _is_custom_angle(force):
            return isinstance(force, openmm.CustomAngleForce)

        assert len([f for f in system.getForces() if _is_custom_angle(f)]) == 1

        for force in system.getForces():
            if _is_custom_angle(force):
                assert force.getEnergyFunction() == "k/2*(cos(theta)-cos(angle))^2"

    def test_openmm_no_valence_forces_with_no_handler(self, sage):
        ethanol = create_ethanol()

        original_system = Interchange.from_smirnoff(sage, [ethanol]).to_openmm(
            combine_nonbonded_forces=True,
        )
        assert original_system.getNumForces() == 4

        sage.deregister_parameter_handler("Constraints")
        sage.deregister_parameter_handler("Bonds")

        no_bonds = Interchange.from_smirnoff(sage, [ethanol]).to_openmm(
            combine_nonbonded_forces=True,
        )
        assert no_bonds.getNumForces() == 3

        sage.deregister_parameter_handler("Angles")

        no_angles = Interchange.from_smirnoff(sage, [ethanol]).to_openmm(
            combine_nonbonded_forces=True,
        )
        assert no_angles.getNumForces() == 2

    def test_openmm_only_electrostatics_no_vdw(self):
        force_field_only_charges = ForceField(get_test_file_path("no_vdw.offxml"))
        molecule = Molecule.from_smiles("[H][Cl]")

        system = Interchange.from_smirnoff(
            force_field_only_charges,
            [molecule],
        ).to_openmm(
            combine_nonbonded_forces=True,
        )

        assert system.getForce(0).getParticleParameters(0)[0]._value == 1.0
        assert system.getForce(0).getParticleParameters(1)[0]._value == -1.0

    def test_nonstandard_cutoffs_match(self, sage):
        """Test that multiple nonbonded forces use the same cutoff."""
        topology = Molecule.from_smiles("C").to_topology()
        topology.box_vectors = unit.Quantity([4, 4, 4], unit.nanometer)

        cutoff = unit.Quantity(1.555, unit.nanometer)

        sage["vdW"].cutoff = cutoff

        interchange = Interchange.from_smirnoff(
            force_field=sage,
            topology=topology,
        )

        system = interchange.to_openmm(combine_nonbonded_forces=False)

        # For now, just make sure all non-bonded forces use the vdW handler's cutoff
        for force in system.getForces():
            if type(force) in (openmm.NonbondedForce, openmm.CustomNonbondedForce):
                assert force.getCutoffDistance().value_in_unit(
                    openmm.unit.nanometer,
                ) == pytest.approx(cutoff.m_as(unit.nanometer))


class TestOpenMMSwitchingFunction(_BaseTest):
    def test_switching_function_applied(self, sage, basic_top):
        out = Interchange.from_smirnoff(force_field=sage, topology=basic_top).to_openmm(
            combine_nonbonded_forces=True,
        )

        found_force = False
        for force in out.getForces():
            if isinstance(force, openmm.NonbondedForce):
                found_force = True
                assert force.getUseSwitchingFunction()
                assert force.getSwitchingDistance().value_in_unit(
                    openmm.unit.angstrom,
                ) == pytest.approx(8), force.getSwitchingDistance()

        assert found_force, "NonbondedForce not found in system"

    def test_switching_function_not_applied(self, sage, basic_top):
        sage["vdW"].switch_width = 0.0 * unit.angstrom

        out = Interchange.from_smirnoff(force_field=sage, topology=basic_top).to_openmm(
            combine_nonbonded_forces=True,
        )

        found_force = False
        for force in out.getForces():
            if isinstance(force, openmm.NonbondedForce):
                found_force = True
                assert not force.getUseSwitchingFunction()
                assert force.getSwitchingDistance() == -1 * openmm.unit.nanometer

        assert found_force, "NonbondedForce not found in system"

    def test_switching_function_nonstandard(self, sage, basic_top):
        sage["vdW"].switch_width = 0.12345 * unit.angstrom

        out = Interchange.from_smirnoff(force_field=sage, topology=basic_top).to_openmm(
            combine_nonbonded_forces=True,
        )

        found_force = False
        for force in out.getForces():
            if isinstance(force, openmm.NonbondedForce):
                found_force = True
                assert force.getUseSwitchingFunction()
                assert (
                    force.getSwitchingDistance() - (9 - 0.12345) * openmm.unit.angstrom
                ) < 1e-10 * openmm.unit.angstrom

        assert found_force, "NonbondedForce not found in system"


class TestOpenMMWithPlugins(TestDoubleExponential):
    pytest.importorskip("deforcefields")

    def test_combine_compatibility(self, de_force_field):
        out = Interchange.from_smirnoff(
            force_field=de_force_field,
            topology=[Molecule.from_smiles("CO")],
        )

        with pytest.raises(
            PluginCompatibilityError,
            match="failed a compatibility check",
        ) as exception:
            out.to_openmm(combine_nonbonded_forces=True)

        assert isinstance(exception.value.__cause__, AssertionError)

    def test_double_exponential_create_simulation(self, de_force_field):
        from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper

        molecule = Molecule.from_smiles("CCO")
        molecule.generate_conformers(n_conformers=1)
        topology = molecule.to_topology()
        topology.box_vectors = unit.Quantity([4, 4, 4], unit.nanometer)

        out = Interchange.from_smirnoff(
            de_force_field,
            topology,
        )

        system = out.to_openmm(combine_nonbonded_forces=False)

        simulation = openmm.app.Simulation(
            to_openmm_topology(out),
            system,
            openmm.LangevinIntegrator(300, 1, 0.002),
            openmm.Platform.getPlatformByName("CPU"),
        )

        simulation.context.setPositions(
            to_openmm_positions(out, include_virtual_sites=False),
        )
        simulation.context.setPeriodicBoxVectors(*out.box.to_openmm())

        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().in_units_of(openmm.unit.kilojoule_per_mole)

        if OpenEyeToolkitWrapper.is_available():
            expected_energy = 13.591709748611304
        else:
            expected_energy = 37.9516622967221

        # Different operating systems report different energies around 0.001 kJ/mol,
        # locally testing this should enable something like 1e-6 kJ/mol
        assert abs(energy._value - expected_energy) < 3e-3


@pytest.mark.slow()
class TestOpenMMVirtualSites(_BaseTest):
    @pytest.fixture()
    def sage_with_sigma_hole(self, sage):
        """Fixture that loads an SMIRNOFF XML with a C-Cl sigma hole."""
        # TODO: Move this into BaseTest to that GROMACS and others can access it
        virtual_site_handler = VirtualSiteHandler(version=0.3)

        sigma_type = VirtualSiteHandler.VirtualSiteType(
            name="EP",
            smirks="[#6:1]-[#17:2]",
            distance=1.4 * unit.angstrom,
            type="BondCharge",
            match="once",
            charge_increment1=0.1 * unit.elementary_charge,
            charge_increment2=0.2 * unit.elementary_charge,
        )

        virtual_site_handler.add_parameter(parameter=sigma_type)
        sage.register_parameter_handler(virtual_site_handler)

        return sage

    @pytest.fixture()
    def sage_with_monovalent_lone_pair(self, sage):
        """Fixture that loads an SMIRNOFF XML for argon"""
        virtual_site_handler = VirtualSiteHandler(version=0.3)

        carbonyl_type = VirtualSiteHandler.VirtualSiteMonovalentLonePairType(
            name="EP",
            smirks="[O:1]=[C:2]-[C:3]",
            distance=0.3 * unit.angstrom,
            type="MonovalentLonePair",
            match="once",
            outOfPlaneAngle=0.0 * unit.degree,
            inPlaneAngle=120.0 * unit.degree,
            charge_increment1=0.05 * unit.elementary_charge,
            charge_increment2=0.1 * unit.elementary_charge,
            charge_increment3=0.15 * unit.elementary_charge,
        )

        virtual_site_handler.add_parameter(parameter=carbonyl_type)
        sage.register_parameter_handler(virtual_site_handler)

        return sage

    def test_valence_term_paticle_index_offsets(self):
        # Use a questionable version of TIP5P that includes angle parameters, since that's what's being tested
        tip5p_offxml = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <LibraryCharges version="0.3">
        <LibraryCharge
            name="tip5p"
            smirks="[#1:1]-[#8X2H2+0:2]-[#1:3]"
            charge1="0.*elementary_charge"
            charge2="0.*elementary_charge"
            charge3="0.*elementary_charge"/>
    </LibraryCharges>
    <vdW
        version="0.3"
        potential="Lennard-Jones-12-6"
        combining_rules="Lorentz-Berthelot"
        scale12="0.0"
        scale13="0.0"
        scale14="0.5"
        scale15="1.0"
        switch_width="0.0 * angstrom"
        cutoff="9.0 * angstrom" method="cutoff">
            <Atom
                smirks="[#1:1]-[#8X2H2+0]-[#1]"
                epsilon="0.*mole**-1*kilojoule"
                sigma="1.0 * nanometer"/>
            <Atom
                smirks="[#1]-[#8X2H2+0:1]-[#1]"
                epsilon="0.66944*mole**-1*kilojoule"
                sigma="0.312*nanometer"/>
    </vdW>
    <Bonds
        version="0.4"
        potential="harmonic"
        fractional_bondorder_method="AM1-Wiberg"
        fractional_bondorder_interpolation="linear">
        <Bond
            smirks="[#1:1]-[#8X2H2+0:2]-[#1]"
            length="0.9572*angstrom"
            k="462750.4*nanometer**-2*mole**-1*kilojoule"/>
    </Bonds>
    <Angles version="0.3" potential="harmonic">
        <Angle
            smirks="[#1:1]-[#8X2H2+0:2]-[#1:3]"
            angle="1.82421813418*radian"
            k="836.8*mole**-1*radian**-2*kilojoule"
            id="a1"/>
    </Angles>
    <VirtualSites version="0.3">
        <VirtualSite
            type="DivalentLonePair"
            name="EP"
            smirks="[#1:2]-[#8X2H2+0:1]-[#1:3]"
            distance="0.70 * angstrom"
            charge_increment1="0.0*elementary_charge"
            charge_increment2="0.1205*elementary_charge"
            charge_increment3="0.1205*elementary_charge"
            sigma="10.0*angstrom"
            epsilon="0.0*kilocalories_per_mole"
            outOfPlaneAngle="54.71384225*degree"
            match="all_permutations" >
        </VirtualSite>
    </VirtualSites>
    <Electrostatics
        version="0.3"
        method="PME"
        scale12="0.0"
        scale13="0.0"
        scale14="0.833333"
        scale15="1.0"
        switch_width="0.0 * angstrom"
        cutoff="9.0 * angstrom"/>
</SMIRNOFF>
"""
        tip5p = ForceField(tip5p_offxml)
        water = Molecule.from_mapped_smiles("[H:2][O:1][H:3]")

        out = Interchange.from_smirnoff(tip5p, [water, water]).to_openmm(
            combine_nonbonded_forces=True,
        )

        assert out.getNumForces() == 3

        for force in out.getForces():
            if isinstance(force, openmm.HarmonicAngleForce):
                p1, p2, p3, _, _ = force.getAngleParameters(1)

                assert p1 == 6
                assert p2 == 5
                assert p3 == 7


class TestOpenMMVirtualSiteExclusions(_BaseTest):
    def test_tip5p_num_exceptions(self):
        tip5p = ForceField(get_test_file_path("tip5p.offxml"))
        water = Molecule.from_smiles("O")
        water.generate_conformers(n_conformers=1)

        out = Interchange.from_smirnoff(tip5p, [water]).to_openmm(
            combine_nonbonded_forces=True,
        )

        # In a TIP5P water    expected exceptions include (total 10)
        #
        # V(3)  V(4)          Oxygen to hydrogens and particles (4)
        #    \ /                - (0, 1), (0, 2), (0, 3), (0, 4)
        #     O(0)            Hyrogens to virtual particles (4)
        #    / \                - (1, 3), (1, 4), (2, 3), (2, 4)
        # H(1)  H(2)          Hydrogens and virtual particles to each other (2)
        #                       - (1, 2), (3, 4)

        for force in out.getForces():
            if isinstance(force, openmm.NonbondedForce):
                assert force.getNumExceptions() == 10

    def test_dichloroethane_exceptions(self, sage):
        """Test a case in which a parent's 1-4 exceptions must be 'imported'."""
        from openff.toolkit.tests.mocking import VirtualSiteMocking

        # This molecule has heavy atoms with indices (1-indexed) CL1, C2, C3, Cl4,
        # resulting in 1-4 interactions between the Cl-Cl pair and some Cl-H pairs
        dichloroethane = Molecule.from_mapped_smiles(
            "[Cl:1][C:2]([H:5])([H:6])[C:3]([H:7])([H:8])[Cl:4]",
        )

        # This parameter pulls 0.1 and 0.2e from Cl (parent) and C, respectively, and has
        # LJ parameters of 4 A, 3 kJ/mol
        parameter = VirtualSiteMocking.bond_charge_parameter("[Cl:1]-[C:2]")

        handler = VirtualSiteHandler(version="0.3")
        handler.add_parameter(parameter=parameter)

        sage.register_parameter_handler(handler)

        system = Interchange.from_smirnoff(sage, [dichloroethane]).to_openmm(
            combine_nonbonded_forces=True,
        )

        assert system.isVirtualSite(8)
        assert system.isVirtualSite(9)

        non_bonded_force = [
            f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)
        ][0]

        for exception_index in range(non_bonded_force.getNumExceptions()):
            p1, p2, q, sigma, epsilon = non_bonded_force.getExceptionParameters(
                exception_index,
            )
            if p2 == 8:
                # Parent Cl, adjacent C and its bonded H, and the 1-3 C
                if p1 in (0, 1, 2, 4, 5):
                    assert q._value == epsilon._value == 0.0
                # 1-4 Cl or 1-4 Hs
                if p1 in (3, 6, 7):
                    for value in (q, sigma, epsilon):
                        assert value._value != 0, (q, sigma, epsilon)
            if p2 == 9:
                if p1 in (3, 1, 2, 6, 7):
                    assert q._value == epsilon._value == 0.0
                if p1 in (0, 4, 5):
                    for value in (q, sigma, epsilon):
                        assert value._value != 0, (q, sigma, epsilon)


class TestToOpenMMTopology(_BaseTest):
    def test_num_virtual_sites(self):
        tip4p = ForceField("openff-2.0.0.offxml", get_test_file_path("tip4p.offxml"))
        water = Molecule.from_smiles("O")

        out = Interchange.from_smirnoff(tip4p, [water])

        assert _get_num_virtual_sites(to_openmm_topology(out)) == 1

        # TODO: Monkeypatch Topology.to_openmm() and emit a warning when it seems
        #       to be used while virtual sites are present in a handler
        assert _get_num_virtual_sites(out.topology.to_openmm()) == 0

    def test_interchange_method(self):
        """
        Ensure similar-ish behavior between `to_openmm_topology` as a standalone function
        and as the wrapped method of the same name on the `Interchange` class.
        """
        tip4p = ForceField("openff-2.0.0.offxml", get_test_file_path("tip4p.offxml"))
        topology = Molecule.from_smiles("O").to_topology()
        topology.box_vectors = unit.Quantity([4, 4, 4], unit.nanometer)

        out = Interchange.from_smirnoff(tip4p, topology)

        _compare_openmm_topologies(out.to_openmm_topology(), to_openmm_topology(out))

    @pytest.mark.parametrize("ensure_unique_atom_names", [True, "residues", "chains"])
    def test_assign_unique_atom_names(self, ensure_unique_atom_names):
        """
        Ensure that OFF topologies with no pre-existing atom names have unique
        atom names applied when being converted to openmm
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        benzene = Molecule.from_smiles("c1ccccc1")
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])

        # This test uses molecules with no hierarchy schemes, so the parametrized
        # ensure_unique_atom_names values should behave identically.
        assert not any(
            [mol._hierarchy_schemes for mol in off_topology.molecules],
        ), "Test assumes no hierarchy schemes"

        sage = ForceField("openff-2.0.0.offxml")
        interchange = Interchange.from_smirnoff(sage, off_topology)

        omm_topology = interchange.to_openmm_topology(
            ensure_unique_atom_names=ensure_unique_atom_names,
        )
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        # There should be 6 unique Cs, 6 unique Hs, and 1 unique O, for a total of 13 unique atom names
        assert len(atom_names) == 13

    @pytest.mark.parametrize("ensure_unique_atom_names", [True, "residues", "chains"])
    def test_assign_some_unique_atom_names(self, ensure_unique_atom_names):
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

        # This test uses molecules with no hierarchy schemes, so the parametrized
        # ensure_unique_atom_names values should behave identically.
        assert not any(
            [mol._hierarchy_schemes for mol in off_topology.molecules],
        ), "Test assumes no hierarchy schemes"

        sage = ForceField("openff-2.0.0.offxml")
        interchange = Interchange.from_smirnoff(sage, off_topology)

        omm_topology = interchange.to_openmm_topology(
            ensure_unique_atom_names=ensure_unique_atom_names,
        )
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        # There should be 9 "ATOM#"-labeled atoms, 6 unique Cs, and 6 unique Hs,
        # for a total of 21 unique atom names
        assert len(atom_names) == 21

    @pytest.mark.parametrize("ensure_unique_atom_names", [True, "residues", "chains"])
    def test_assign_unique_atom_names_some_duplicates(self, ensure_unique_atom_names):
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

        # This test uses molecules with no hierarchy schemes, so the parametrized
        # ensure_unique_atom_names values should behave identically.
        assert not any(
            [mol._hierarchy_schemes for mol in off_topology.molecules],
        ), "Test assumes no hierarchy schemes"

        sage = ForceField("openff-2.0.0.offxml")
        interchange = Interchange.from_smirnoff(sage, off_topology)

        omm_topology = interchange.to_openmm_topology(
            ensure_unique_atom_names=ensure_unique_atom_names,
        )
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)

        # There should be  12 "AT#"-labeled atoms (from benzene), 2 unique Cs,
        # 1 unique O, and 6 unique Hs, for a total of 21 unique atom names
        assert len(atom_names) == 21

    def test_do_not_assign_unique_atom_names(self):
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

        sage = ForceField("openff-2.0.0.offxml")
        interchange = Interchange.from_smirnoff(sage, off_topology)

        omm_topology = interchange.to_openmm_topology(ensure_unique_atom_names=False)
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        # There should be 9 atom named "eth_test", 1 atom named "bzn_test",
        # and 12 atoms named "", for a total of 3 unique atom names
        assert len(atom_names) == 3

    @pytest.mark.parametrize("explicit_arg", [True, False])
    def test_preserve_per_residue_unique_atom_names(self, explicit_arg):
        """
        Test that to_openmm preserves atom names that are unique per-residue by default
        """
        # Create a topology from a capped dialanine
        peptide = Molecule.from_polymer_pdb(
            get_data_file_path("proteins/MainChain_ALA_ALA.pdb"),
        )
        off_topology = Topology.from_molecules([peptide])

        # Assert the test's assumptions
        _ace, ala1, ala2, _nme = off_topology.hierarchy_iterator("residues")
        assert [a.name for a in ala1.atoms] == [
            a.name for a in ala2.atoms
        ], "Test assumes both alanines have same atom names"

        for res in off_topology.hierarchy_iterator("residues"):
            res_atomnames = [atom.name for atom in res.atoms]
            assert len(set(res_atomnames)) == len(
                res_atomnames,
            ), f"Test assumes atom names are already unique per-residue in {res}"

        # Record the initial atom names
        init_atomnames = [str(atom.name) for atom in off_topology.atoms]

        sage = ForceField("openff-2.0.0.offxml")
        interchange = Interchange.from_smirnoff(sage, off_topology)

        # Perform the test
        if explicit_arg:
            omm_topology = interchange.to_openmm_topology(
                ensure_unique_atom_names="residues",
            )
        else:
            omm_topology = interchange.to_openmm_topology()

        # Check that the atom names were preserved
        final_atomnames = [str(atom.name) for atom in omm_topology.atoms()]
        assert final_atomnames == init_atomnames

    @pytest.mark.parametrize("explicit_arg", [True, False])
    def test_generate_per_residue_unique_atom_names(self, explicit_arg):
        """
        Test that to_openmm generates atom names that are unique per-residue
        """
        # Create a topology from a capped dialanine
        peptide = Molecule.from_polymer_pdb(
            get_data_file_path("proteins/MainChain_ALA_ALA.pdb"),
        )
        off_topology = Topology.from_molecules([peptide])

        # Remove atom names from some residues, make others have duplicate atom names
        ace, ala1, ala2, nme = off_topology.hierarchy_iterator("residues")
        for atom in ace.atoms:
            atom._name = None
        for atom in ala1.atoms:
            atom.name = ""
        for atom in ala2.atoms:
            atom.name = "ATX2"
        for atom in nme.atoms:
            if atom.name == "H2":
                atom.name = "H1"
                break

        # Assert assumptions
        for res in off_topology.hierarchy_iterator("residues"):
            res_atomnames = [atom.name for atom in res.atoms]
            assert len(set(res_atomnames)) != len(
                res_atomnames,
            ), f"Test assumes atom names are not unique per-residue in {res}"
        assert off_topology.n_atoms == 32, "Test assumes topology has 32 atoms"

        sage = ForceField("openff-2.0.0.offxml")
        interchange = Interchange.from_smirnoff(sage, off_topology)

        # Perform the test
        if explicit_arg:
            omm_topology = interchange.to_openmm_topology(
                ensure_unique_atom_names="residues",
            )
        else:
            omm_topology = interchange.to_openmm_topology()

        # Check that the atom names are now unique per-residue but not per-molecule
        for res in omm_topology.residues():
            res_atomnames = [atom.name for atom in res.atoms()]
            assert len(set(res_atomnames)) == len(
                res_atomnames,
            ), f"Final atom names are not unique in residue {res}"

        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        assert (
            len(atom_names) < 32
        ), "There should be duplicate atom names in this output topology"

    @pytest.mark.parametrize("ensure_unique_atom_names", ["chains", True])
    def test_generate_per_molecule_unique_atom_names_with_residues(
        self,
        ensure_unique_atom_names,
    ):
        """
        Test that to_openmm can generate atom names that are unique per-molecule
        when the topology has residues
        """
        # Create a topology from a capped dialanine
        peptide = Molecule.from_polymer_pdb(
            get_data_file_path("proteins/MainChain_ALA_ALA.pdb"),
        )
        off_topology = Topology.from_molecules([peptide])

        # Remove atom names from some residues, make others have duplicate atom names
        ace, ala1, ala2, nme = off_topology.hierarchy_iterator("residues")
        for atom in ace.atoms:
            atom._name = None
        for atom in ala1.atoms:
            atom.name = ""
        for atom in ala2.atoms:
            atom.name = "ATX2"
        for atom in nme.atoms:
            if atom.name == "H2":
                atom.name = "H1"
                break

        # Assert assumptions
        for res in off_topology.hierarchy_iterator("residues"):
            res_atomnames = [atom.name for atom in res.atoms]
            assert len(set(res_atomnames)) != len(
                res_atomnames,
            ), f"Test assumes atom names are not unique per-residue in {res}"
        assert off_topology.n_atoms == 32, "Test assumes topology has 32 atoms"

        sage = ForceField("openff-2.0.0.offxml")
        interchange = Interchange.from_smirnoff(sage, off_topology)

        # Perform the test
        omm_topology = interchange.to_openmm_topology(
            ensure_unique_atom_names=ensure_unique_atom_names,
        )

        # Check that the atom names are now unique across the topology (of 1 molecule)
        atom_names = set()
        for atom in omm_topology.atoms():
            atom_names.add(atom.name)
        assert (
            len(atom_names) == 32
        ), "There should not be duplicate atom names in this output topology"

    @pytest.mark.parametrize(
        "ensure_unique_atom_names",
        [True, "residues", "chains", False],
    )
    def test_to_openmm_copies_molecules(self, ensure_unique_atom_names):
        """
        Check that generating new atom names doesn't affect the input topology
        """
        # Create OpenFF topology with 1 ethanol and 2 benzenes.
        ethanol = Molecule.from_smiles("CCO")
        for atom in ethanol.atoms:
            atom.name = f"AT{atom.molecule_atom_index}"
        benzene = Molecule.from_smiles("c1ccccc1")
        off_topology = Topology.from_molecules(molecules=[ethanol, benzene, benzene])

        # This test uses molecules with no hierarchy schemes, so the parametrized
        # ensure_unique_atom_names values should behave identically (except False).
        assert not any(
            [mol._hierarchy_schemes for mol in off_topology.molecules],
        ), "Test assumes no hierarchy schemes"

        sage = ForceField("openff-2.0.0.offxml")
        interchange = Interchange.from_smirnoff(sage, off_topology)

        # Record the initial atom names to compare to later
        init_atomnames = [str(atom.name) for atom in interchange.topology.atoms]

        omm_topology = interchange.to_openmm_topology(
            ensure_unique_atom_names=ensure_unique_atom_names,
        )

        # Get the atom names back from the initial molecules after calling to_openmm
        final_atomnames_mols = [
            atom.name for atom in [*ethanol.atoms, *benzene.atoms, *benzene.atoms]
        ]
        # Get the atom names back from the initial topology after calling to_openmm
        final_atomnames_offtop = [atom.name for atom in off_topology.atoms]
        # Get the atom names back from the new OpenMM topology
        final_atomnames_ommtop = [atom.name for atom in omm_topology.atoms()]

        # Check the appropriate properties!
        assert (
            init_atomnames == final_atomnames_mols
        ), "Molecules' atom names were changed"
        assert (
            init_atomnames == final_atomnames_offtop
        ), "Topology's atom names were changed"
        if ensure_unique_atom_names:
            assert (
                init_atomnames != final_atomnames_ommtop
            ), "New atom names should've been generated but weren't"


class TestToOpenMMPositions(_BaseTest):
    def test_missing_positions(self):
        with pytest.raises(
            MissingPositionsError,
            match=r"are required.*\.positions=None",
        ):
            to_openmm_positions(Interchange())

    @pytest.mark.parametrize("include_virtual_sites", [True, False])
    def test_positions_basic(self, include_virtual_sites):
        force_field = ForceField(
            "openff-2.0.0.offxml",
            get_test_file_path(
                "tip4p.offxml" if include_virtual_sites else "tip3p.offxml",
            ),
        )
        water = Molecule.from_smiles("O")
        water.generate_conformers(n_conformers=1)

        out = Interchange.from_smirnoff(force_field, [water])

        positions = to_openmm_positions(
            out,
            include_virtual_sites=include_virtual_sites,
        )

        assert isinstance(positions, openmm.unit.Quantity)
        assert positions.shape == (4, 3) if include_virtual_sites else (3, 3)

        numpy.testing.assert_allclose(
            positions.value_in_unit(openmm.unit.angstrom)[:3],
            water.conformers[0].m_as(unit.angstrom),
        )

    @pytest.mark.parametrize("include_virtual_sites", [True, False])
    def test_given_positions(self, include_virtual_sites):
        """Test issue #616"""
        force_field = ForceField(
            "openff-2.0.0.offxml",
            get_test_file_path(
                "tip4p.offxml" if include_virtual_sites else "tip3p.offxml",
            ),
        )

        water = Molecule.from_smiles("O")
        water.generate_conformers(n_conformers=1)

        topology = Topology.from_molecules([water, water])
        out = Interchange.from_smirnoff(force_field, topology)

        # Approximate conformer position with a duplicate 5 A away in x
        out.positions = unit.Quantity(
            numpy.array(
                [
                    [0.85, 1.17, 0.84],
                    [1.51, 0.47, 0.75],
                    [0.0, 0.71, 0.76],
                    [5.85, 1.17, 0.84],
                    [6.51, 0.47, 0.75],
                    [5.0, 0.71, 0.76],
                ],
            ),
            unit.angstrom,
        )

        positions = to_openmm_positions(
            out,
            include_virtual_sites=include_virtual_sites,
        )

        assert isinstance(positions, openmm.unit.Quantity)

        # Number of particles per molecule
        n = 3 + int(include_virtual_sites)

        assert numpy.allclose(
            (positions[n:][:3] - positions[:n][:3]).value_in_unit(openmm.unit.angstrom),
            numpy.array([[5, 0, 0], [5, 0, 0], [5, 0, 0]]),
        )


class TestOpenMMToPDB(_BaseTest):
    def test_to_pdb(self, sage):
        import mdtraj as md

        molecule = Molecule.from_smiles("O")

        out = Interchange.from_smirnoff(sage, molecule.to_topology())

        with pytest.raises(MissingPositionsError):
            out.to_pdb("file_should_not_exist.pdb")

        molecule.generate_conformers(n_conformers=1)
        out.positions = molecule.conformers[0]

        out.to_pdb("out.pdb")

        md.load("out.pdb")

        with pytest.raises(UnsupportedExportError):
            out.to_pdb("file_should_not_exist.pdb", writer="magik")


class TestBuckingham:
    def test_water_with_virtual_sites(self):
        force_field = ForceField(
            get_test_file_path("buckingham_virtual_sites.offxml"),
            load_plugins=True,
        )

        water = Molecule.from_mapped_smiles("[H:2][O:1][H:3]")
        water.generate_conformers(n_conformers=1)
        topology = water.to_topology()

        interchange = Interchange.from_smirnoff(
            force_field=force_field,
            topology=topology,
            box=[4, 4, 4],
        )

        with pytest.raises(PluginCompatibilityError):
            interchange.to_openmm(combine_nonbonded_forces=True)

        system = interchange.to_openmm(combine_nonbonded_forces=False)

        assert system.getNumForces() == 4

        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                electrostatics = force
                continue
            elif isinstance(force, openmm.CustomNonbondedForce):
                vdw = force
                continue
            elif isinstance(force, openmm.CustomBondForce):
                if "qq" in force.getEnergyFunction():
                    electrostatics14 = force
                    continue

        assert system.getNumParticles() == 4

        masses = [15.99943, 1.007947, 1.007947, 0.0]

        for particle_index in range(system.getNumParticles()):
            assert system.isVirtualSite(particle_index) == (particle_index == 3)
            assert system.getParticleMass(particle_index)._value == pytest.approx(
                masses[particle_index],
            )

        charges = openmm.unit.Quantity(
            [0.0, 0.53254, 0.53254, -1.06508],
            openmm.unit.elementary_charge,
        )

        for index, charge in enumerate(charges):
            assert electrostatics.getParticleParameters(index)[0] == charge

        for index in range(vdw.getNumParticles()):
            parameters = vdw.getParticleParameters(index)
            for p in parameters:
                assert (p == 0) == (index > 0)

        assert vdw.getParticleParameters(0) == (1600000.0, 42.0, 0.003)

        # This test should be replaced with one that uses a more complex
        # system than a single water molecule and look at vdw14 force
        assert electrostatics14.getNumBonds() == 0

        with pytest.raises(PluginCompatibilityError):
            get_openmm_energies(interchange, combine_nonbonded_forces=True)

        with pytest.warns(
            UserWarning,
            match="energies from split forces with virtual sites",
        ):
            assert not math.isnan(
                get_openmm_energies(
                    interchange,
                    combine_nonbonded_forces=False,
                ).total_energy.m,
            )


class TestGBSA(_BaseTest):
    def test_create_gbsa(self):
        force_field = ForceField(
            "openff-2.0.0.offxml",
            get_test_file_path("gbsa.offxml"),
        )

        molecule = Molecule.from_smiles("CCO")
        molecule.generate_conformers(n_conformers=1)

        interchange = Interchange.from_smirnoff(
            force_field=force_field,
            topology=molecule.to_topology(),
            box=[4, 4, 4] * unit.nanometer,
        )

        assert get_openmm_energies(interchange).total_energy is not None

    def test_cannot_split_nonbonded_forces(self):
        force_field = ForceField(
            "openff-2.0.0.offxml",
            get_test_file_path("gbsa.offxml"),
        )

        force_field["Electrostatics"]
        molecule = Molecule.from_smiles("CCO")
        molecule.generate_conformers(n_conformers=1)

        with pytest.raises(UnsupportedExportError, match="exactly one"):
            Interchange.from_smirnoff(
                force_field=force_field,
                topology=molecule.to_topology(),
                box=[4, 4, 4] * unit.nanometer,
            ).to_openmm(combine_nonbonded_forces=False)

    def test_no_cutoff(self):
        force_field = ForceField(
            "openff-2.0.0.offxml",
            get_test_file_path("gbsa.offxml"),
        )

        force_field["Electrostatics"]
        molecule = Molecule.from_smiles("CCO")
        molecule.generate_conformers(n_conformers=1)

        system = Interchange.from_smirnoff(
            force_field=force_field,
            topology=molecule.to_topology(),
            box=None,
        ).to_openmm(combine_nonbonded_forces=True)

        for force in system.getForces():
            if isinstance(force, openmm.CustomGBForce):
                assert force.getNonbondedMethod() == openmm.CustomGBForce.NoCutoff
                # This should be set to OpenMM's default, though not used
                assert force.getCutoffDistance() == 1.0 * openmm.unit.nanometer
                break

        else:
            raise Exception
