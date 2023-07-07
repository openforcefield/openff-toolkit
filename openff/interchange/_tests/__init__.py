"""Assorted utilities used in tests."""
import pathlib
import sys
from collections import defaultdict
from typing import DefaultDict, Optional

import numpy as np
import openmm
import pytest
from openff.toolkit.tests.create_molecules import create_ammonia, create_ethanol
from openff.toolkit.tests.utils import get_data_file_path
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.units import unit
from openmm import unit as openmm_unit

from openff.interchange import Interchange
from openff.interchange.drivers.gromacs import _find_gromacs_executable
from openff.interchange.drivers.lammps import _find_lammps_executable

if sys.version_info >= (3, 10):
    from importlib import resources
else:
    import importlib_resources as resources


def get_test_file_path(test_file: str) -> pathlib.Path:
    """Given a filename in the collection of data files, return its full path."""
    test_dir_path = get_test_files_dir_path()
    test_file_path = test_dir_path / test_file

    if test_file_path.is_file():
        return test_file_path
    else:
        raise FileNotFoundError(
            f"could not file file {test_file} in path {test_file_path}",
        )


def get_test_files_dir_path(dirname: Optional[str] = None) -> pathlib.Path:
    """Given a directory with a collection of test data files, return its full path."""
    dir_path = resources.files("openff.interchange._tests.data")

    if dirname:
        test_dir: pathlib.PosixPath = dir_path / dirname  # type: ignore[assignment]
    else:
        test_dir = dir_path  # type: ignore[assignment]

    if test_dir.is_dir():
        return test_dir
    else:
        raise NotADirectoryError(
            f"Provided directory {dirname} doesn't exist in {dir_path}",
        )


class _BaseTest:
    @pytest.fixture(autouse=True)
    def _initdir(self, tmpdir):
        tmpdir.chdir()

    # TODO: group fixtures up as dicts, i.e. argon['forcefield'], argon['topology'], ...
    @pytest.fixture()
    def argon_ff(self):
        """Fixture that loads an SMIRNOFF XML for argon."""
        return ForceField(get_test_file_path("argon.offxml"))

    @pytest.fixture()
    def argon(self):
        """Fixture that builds a simple arogon topology."""
        argon = Molecule()
        argon.add_atom(
            atomic_number=18,
            formal_charge=0,
            is_aromatic=False,
        )
        return argon

    @pytest.fixture()
    def ammonia_ff(self):
        """Fixture that loads an SMIRNOFF XML for ammonia."""
        return ForceField(get_test_file_path("ammonia.offxml"))

    @pytest.fixture()
    def ammonia(self):
        return create_ammonia()

    @pytest.fixture()
    def ethanol(self):
        return create_ethanol()

    @pytest.fixture()
    def basic_top(self):
        top = Molecule.from_smiles("C").to_topology()
        top.box_vectors = unit.Quantity([5, 5, 5], unit.nanometer)
        return top

    @pytest.fixture()
    def ammonia_top(self, ammonia):
        """Fixture that builds a simple ammonia topology."""
        return Topology.from_molecules(4 * [ammonia])

    @pytest.fixture()
    def ethanol_top(self, ethanol):
        """Fixture that builds a simple four ethanol topology."""
        return Topology.from_molecules(4 * [ethanol])

    @pytest.fixture()
    def sage(self):
        return ForceField("openff-2.0.0.offxml")

    @pytest.fixture()
    def sage_unconstrained(self):
        return ForceField("openff_unconstrained-2.0.0.offxml")

    @pytest.fixture()
    def mainchain_ala(self):
        molecule = Molecule.from_file(get_data_file_path("proteins/MainChain_ALA.sdf"))
        molecule._add_default_hierarchy_schemes()
        molecule.perceive_residues()
        molecule.perceive_hierarchy()

        return molecule

    @pytest.fixture()
    def mainchain_arg(self):
        molecule = Molecule.from_file(get_data_file_path("proteins/MainChain_ARG.sdf"))
        molecule._add_default_hierarchy_schemes()
        molecule.perceive_residues()
        molecule.perceive_hierarchy()

        return molecule

    @pytest.fixture()
    def two_peptides(self, mainchain_ala, mainchain_arg):
        return Topology.from_molecules([mainchain_ala, mainchain_arg])

    @pytest.fixture()
    def tip3p_xml(self):
        # Modified (added Electrostatics tag) from below link
        # https://github.com/openforcefield/openff-toolkit/blob/0.10.2/openff/toolkit/data/test_forcefields/tip3p.offxml
        return """<?xml version="1.0" encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <Constraints version="0.3">
    <Constraint smirks="[#1:1]-[*:2]" id="c1"></Constraint>
    <Constraint smirks="[#1:1]-[#8X2H2+0:2]-[#1]"
        id="c-tip3p-H-O" distance="0.9572 * angstrom"></Constraint>
    <Constraint smirks="[#1:1]-[#8X2H2+0]-[#1:2]"
        id="c-tip3p-H-O-H" distance="1.5139006545247014*angstrom"></Constraint>
  </Constraints>
  <vdW
    version="0.3"
    potential="Lennard-Jones-12-6"
    combining_rules="Lorentz-Berthelot"
    scale12="0.0"
    scale13="0.0"
    scale14="0.5"
    scale15="1"
    switch_width="1.0*angstroms"
    cutoff="9.0*angstroms" method="cutoff"
  >
    <Atom
      smirks="[#1]-[#8X2H2+0:1]-[#1]"
      id="n1"
      sigma="0.31507524065751241*nanometers"
      epsilon="0.635968*kilojoules_per_mole"
    />
    <Atom
      smirks="[#1:1]-[#8X2H2+0]-[#1]"
      id="n2"
      sigma="1*nanometers"
      epsilon="0*kilojoules_per_mole" />
  </vdW>
  <Electrostatics
    version="0.3"
    method="PME"
    scale12="0.0"
    scale13="0.0"
    scale14="0.833333"
    scale15="1.0"
    switch_width="0.0*angstrom"
    cutoff="9.0*angstrom"
  ></Electrostatics>
  <LibraryCharges version="0.3">
    <LibraryCharge
        name="TIP3P"
        smirks="[#1:1]-[#8X2H2+0:2]-[#1:3]"
        charge1="0.417*elementary_charge"
        charge2="-0.834*elementary_charge"
        charge3="0.417*elementary_charge"/>
  </LibraryCharges>
</SMIRNOFF>"""

    @pytest.fixture()
    def tip3p(self, tip3p_xml):
        return ForceField(tip3p_xml)

    @pytest.fixture()
    def tip3p_missing_electrostatics_xml(self):
        # Stripped from below, follow link for details
        # https://github.com/openforcefield/openff-toolkit/blob/0.10.2/openff/toolkit/data/test_forcefields/tip3p.offxml
        return """<?xml version="1.0" encoding='ASCII'?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  <vdW
    version="0.3"
    potential="Lennard-Jones-12-6"
    combining_rules="Lorentz-Berthelot"
    scale12="0.0"
    scale13="0.0"
    scale14="0.5"
    scale15="1"
    switch_width="1.0*angstroms"
    cutoff="9.0*angstroms" method="cutoff"
  >
    <Atom
      smirks="[#1]-[#8X2H2+0:1]-[#1]"
      id="n1"
      sigma="0.31507524065751241*nanometers"
      epsilon="0.635968*kilojoules_per_mole"
    />
    <Atom
      smirks="[#1:1]-[#8X2H2+0]-[#1]"
      id="n2"
      sigma="1*nanometers"
      epsilon="0*kilojoules_per_mole" />
  </vdW>
  <LibraryCharges version="0.3">
    <LibraryCharge
        name="TIP3P"
        smirks="[#1:1]-[#8X2H2+0:2]-[#1:3]"
        charge1="0.417*elementary_charge"
        charge2="-0.834*elementary_charge"
        charge3="0.417*elementary_charge"/>
  </LibraryCharges>
</SMIRNOFF>"""

    xml_ff_bo_bonds = """<?xml version='1.0' encoding='ASCII'?>
    <SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
      <Bonds version="0.3" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
        <Bond smirks="[#6:1]~[#8:2]" id="bbo1"
            k_bondorder1="100.0*kilocalories_per_mole/angstrom**2"
            k_bondorder2="1000.0*kilocalories_per_mole/angstrom**2"
            length_bondorder1="1.5*angstrom"
            length_bondorder2="1.0*angstrom"/>
      </Bonds>
    </SMIRNOFF>
    """

    @pytest.fixture()
    def methane(self):
        return Molecule.from_smiles("C")

    @pytest.fixture()
    def parsley(self):
        return ForceField("openff-1.0.0.offxml")

    @pytest.fixture()
    def hydrogen_cyanide(self):
        return Molecule.from_mapped_smiles("[H:1][C:2]#[N:3]")

    @pytest.fixture()
    def hydrogen_cyanide_reversed(self):
        return Molecule.from_mapped_smiles("[H:3][C:2]#[N:1]")

    @pytest.fixture()
    def hexane_diol(self):
        molecule = Molecule.from_smiles("OCCCCCCO")
        molecule.assign_partial_charges(partial_charge_method="gasteiger")
        molecule.partial_charges.m
        return molecule

    @pytest.fixture()
    def hydrogen_chloride(self):
        return Molecule.from_mapped_smiles("[Cl:1][H:2]")

    @pytest.fixture()
    def formaldehyde(self):
        return Molecule.from_mapped_smiles("[H:3][C:1]([H:4])=[O:2]")

    @pytest.fixture()
    def acetaldehyde(self):
        return Molecule.from_mapped_smiles(
            "[C:1]([C:2](=[O:3])[H:7])([H:4])([H:5])[H:6]",
        )

    @pytest.fixture()
    def water(self):
        return Molecule.from_mapped_smiles("[H:2][O:1][H:3]")


HAS_GROMACS = _find_gromacs_executable() is not None
HAS_LAMMPS = _find_lammps_executable() is not None

needs_gmx = pytest.mark.skipif(not HAS_GROMACS, reason="Needs GROMACS")
needs_not_gmx = pytest.mark.skipif(
    HAS_GROMACS,
    reason="Needs GROMACS to NOT be installed",
)
needs_lmp = pytest.mark.skipif(not HAS_LAMMPS, reason="Needs LAMMPS")
needs_not_lmp = pytest.mark.skipif(
    HAS_LAMMPS,
    reason="Needs LAMMPS to NOT be installed",
)

kj_nm2_mol = openmm_unit.kilojoule_per_mole / openmm_unit.nanometer**2
kj_rad2_mol = openmm_unit.kilojoule_per_mole / openmm_unit.radian**2


def _get_charges_from_openmm_system(omm_sys: openmm.System):
    for force in omm_sys.getForces():
        if type(force) == openmm.NonbondedForce:
            break
    for idx in range(omm_sys.getNumParticles()):
        param = force.getParticleParameters(idx)
        yield param[0].value_in_unit(openmm_unit.elementary_charge)


def _get_sigma_from_nonbonded_force(
    n_particles: int,
    nonbond_force: openmm.NonbondedForce,
):
    for idx in range(n_particles):
        param = nonbond_force.getParticleParameters(idx)
        yield param[1].value_in_unit(openmm_unit.nanometer)


def _get_epsilon_from_nonbonded_force(
    n_particles: int,
    nonbond_force: openmm.NonbondedForce,
):
    for idx in range(n_particles):
        param = nonbond_force.getParticleParameters(idx)
        yield param[2].value_in_unit(openmm_unit.kilojoule_per_mole)


def _get_lj_params_from_openmm_system(omm_sys: openmm.System):
    for force in omm_sys.getForces():
        if type(force) == openmm.NonbondedForce:
            break
    n_particles = omm_sys.getNumParticles()
    sigmas = np.asarray([*_get_sigma_from_nonbonded_force(n_particles, force)])
    epsilons = np.asarray([*_get_epsilon_from_nonbonded_force(n_particles, force)])

    return sigmas, epsilons


def _get_charges_from_openff_interchange(interchange: Interchange):
    charges_ = [*interchange["Electrostatics"].charges.values()]
    charges = np.asarray([charge.magnitude for charge in charges_])
    return charges


def _create_torsion_dict(torsion_force) -> dict[tuple[int], list[tuple]]:
    torsions: DefaultDict = defaultdict(list)

    for i in range(torsion_force.getNumTorsions()):
        p1, p2, p3, p4, periodicity, phase, k = torsion_force.getTorsionParameters(i)
        key = (p1, p2, p3, p4)
        torsions[key]
        torsions[key].append((periodicity, phase, k))

    return torsions


def _create_bond_dict(bond_force):
    bonds = dict()

    for i in range(bond_force.getNumBonds()):
        p1, p2, length, k = bond_force.getBondParameters(i)
        key = (p1, p2)
        bonds[key] = (length, k)

    return bonds


def _create_angle_dict(angle_force):
    angles = dict()

    for i in range(angle_force.getNumAngles()):
        p1, p2, p3, theta, k = angle_force.getAngleParameters(i)
        key = (p1, p2, p3)
        angles[key] = (theta, k)

    return angles


def _compare_individual_torsions(x, y):
    assert x[0] == y[0]
    assert x[1] == y[1]
    assert (x[2] - y[2]) < 1e-15 * openmm_unit.kilojoule_per_mole


def _compare_torsion_forces(force1, force2):
    sorted1 = _create_torsion_dict(torsion_force=force1)
    sorted2 = _create_torsion_dict(torsion_force=force2)

    assert sum(len(v) for v in sorted1.values()) == force1.getNumTorsions()
    assert sum(len(v) for v in sorted2.values()) == force2.getNumTorsions()
    assert len(sorted1) == len(sorted2)

    for key in sorted1:
        for i in range(len(sorted1[key])):
            _compare_individual_torsions(sorted1[key][i], sorted2[key][i])


def _compare_bond_forces(force1, force2):
    assert force1.getNumBonds() == force2.getNumBonds()

    bonds1 = _create_bond_dict(force1)
    bonds2 = _create_bond_dict(force2)

    for key in bonds1:
        length_diff = bonds2[key][0] - bonds1[key][0]
        assert (
            abs(length_diff) < 1e-15 * openmm_unit.nanometer
        ), f"Bond lengths differ by {length_diff}"
        k_diff = bonds2[key][1] - bonds1[key][1]
        assert abs(k_diff) < 1e-9 * kj_nm2_mol, f"bond k differ by {k_diff}"


def _compare_angle_forces(force1, force2):
    assert force1.getNumAngles() == force2.getNumAngles()

    angles1 = _create_angle_dict(force1)
    angles2 = _create_angle_dict(force2)

    for key in angles1:
        angle_diff = angles2[key][0] - angles1[key][0]
        assert (
            abs(angle_diff) < 1e-15 * openmm_unit.radian
        ), f"angles differ by {angle_diff}"
        k_diff = angles2[key][1] - angles1[key][1]
        assert abs(k_diff) < 1e-10 * kj_rad2_mol, f"angle k differ by {k_diff}"


def _compare_nonbonded_settings(force1, force2):
    for attr in dir(force1):
        if not attr.startswith("get") or attr in [
            "getExceptionParameterOffset",
            "getExceptionParameters",
            "getGlobalParameterDefaultValue",
            "getGlobalParameterName",
            "getLJPMEParametersInContext",
            "getPMEParametersInContext",
            "getParticleParameterOffset",
            "getParticleParameters",
            "getForceGroup",
        ]:
            continue
        assert getattr(force1, attr)() == getattr(force2, attr)(), attr


def _compare_nonbonded_parameters(force1, force2):
    assert (
        force1.getNumParticles() == force2.getNumParticles()
    ), "found different number of particles"

    for i in range(force1.getNumParticles()):
        q1, sig1, eps1 = force1.getParticleParameters(i)
        q2, sig2, eps2 = force2.getParticleParameters(i)
        assert (
            abs(q2 - q1) < 1e-8 * openmm_unit.elementary_charge
        ), f"charge mismatch in particle {i}: {q1} vs {q2}"
        assert (
            abs(sig2 - sig1) < 1e-12 * openmm_unit.nanometer
        ), f"sigma mismatch in particle {i}: {sig1} vs {sig2}"
        assert (
            abs(eps2 - eps1) < 1e-12 * openmm_unit.kilojoule_per_mole
        ), f"epsilon mismatch in particle {i}: {eps1} vs {eps2}"


def _compare_exceptions(force1, force2):
    assert (
        force1.getNumExceptions() == force2.getNumExceptions()
    ), "found different number of exceptions"

    for i in range(force1.getNumExceptions()):
        _, _, q1, sig1, eps1 = force1.getExceptionParameters(i)
        _, _, q2, sig2, eps2 = force2.getExceptionParameters(i)
        assert (
            abs(q2 - q1) < 1e-12 * openmm_unit.elementary_charge**2
        ), f"charge mismatch in exception {i}"
        assert (
            abs(sig2 - sig1) < 1e-12 * openmm_unit.nanometer
        ), f"sigma mismatch in exception {i}"
        assert (
            abs(eps2 - eps1) < 1e-12 * openmm_unit.kilojoule_per_mole
        ), f"epsilon mismatch in exception {i}"


def _get_force(openmm_sys: openmm.System, force_type):
    forces = [f for f in openmm_sys.getForces() if type(f) == force_type]

    if len(forces) > 1:
        raise NotImplementedError("Not yet able to process duplicate forces types")
    return forces[0]
