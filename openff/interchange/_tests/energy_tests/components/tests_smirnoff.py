import numpy as np
import openmm
import pytest
from openff.toolkit.tests.utils import get_data_file_path
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff.forcefield import ForceField
from openff.units import unit
from openmm import unit as openmm_unit

from openff.interchange import Interchange
from openff.interchange._tests.unit_tests.smirnoff.test_valence import (
    TestBondOrderInterpolation,
)
from openff.interchange.drivers.openmm import _get_openmm_energies, get_openmm_energies


class TestBondOrderInterpolationEnergies(TestBondOrderInterpolation):
    @pytest.mark.slow()
    def test_basic_bond_order_interpolation_energies(self):
        forcefield = ForceField(
            "test_forcefields/test_forcefield.offxml",
            self.xml_ff_bo_bonds,
        )

        mol = Molecule.from_file(get_data_file_path("molecules/CID20742535_anion.sdf"))
        mol.generate_conformers(n_conformers=1)
        top = mol.to_topology()

        out = Interchange.from_smirnoff(forcefield, top)
        out.box = [4, 4, 4] * unit.nanometer
        out.positions = mol.conformers[0]

        interchange_bond_energy = get_openmm_energies(
            out,
            combine_nonbonded_forces=True,
        ).energies["Bond"]
        toolkit_bond_energy = _get_openmm_energies(
            forcefield.create_openmm_system(top),
            box_vectors=[[4, 0, 0], [0, 4, 0], [0, 0, 4]] * openmm_unit.nanometer,
            positions=mol.conformers[0],
        ).energies["Bond"]

        assert abs(interchange_bond_energy - toolkit_bond_energy).m < 1e-2

        new = out.to_openmm(combine_nonbonded_forces=True)
        ref = forcefield.create_openmm_system(top)

        new_k = []
        new_length = []
        for force in new.getForces():
            if type(force) == openmm.HarmonicBondForce:
                for i in range(force.getNumBonds()):
                    new_k.append(force.getBondParameters(i)[3]._value)
                    new_length.append(force.getBondParameters(i)[2]._value)

        ref_k = []
        ref_length = []
        for force in ref.getForces():
            if type(force) == openmm.HarmonicBondForce:
                for i in range(force.getNumBonds()):
                    ref_k.append(force.getBondParameters(i)[3]._value)
                    ref_length.append(force.getBondParameters(i)[2]._value)

        np.testing.assert_allclose(ref_k, new_k, rtol=3e-5)
