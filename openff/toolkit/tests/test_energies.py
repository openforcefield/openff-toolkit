import json

import numpy as np
import pytest

try:
    import openmm
    from openff.units.openmm import to_openmm
    from openmm import unit as openmm_unit
except ImportError:
    from simtk import unit as openmm_unit, openmm
    from openff.units.simtk import to_simtk as to_openmm

from openff.toolkit.tests.utils import get_data_file_path, requires_rdkit
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField


@requires_rdkit
@pytest.mark.parametrize("constrained", [True, False])
@pytest.mark.parametrize(
    "mol",
    [
        "ethanol.sdf",
        "methane_multiconformer.sdf",
        "CID20742535_anion.sdf",
    ],
)
def test_reference(constrained, mol):
    """Minimal regression test comparing molecule energies to energies computed
    by version 0.8.0 of the toolkit"""
    # TODO: Also test periodic vs. vacuum
    with open(
        get_data_file_path("reference_energies/reference_after_1007.json"), "r"
    ) as fi:
        reference = json.loads(fi.read())

    name = mol + "_"
    if not constrained:
        name += "un"
    name += "constrained"
    reference_energy = reference[name]

    omm_sys, positions, off_top = _build_system(
        mol=mol,
        constrained=constrained,
    )

    simulation = _build_simulation(omm_sys=omm_sys, off_top=off_top)
    derived_energy = _get_energy(simulation=simulation, positions=positions)

    np.testing.assert_almost_equal(
        actual=derived_energy.value_in_unit(openmm_unit.kilojoule_per_mole),
        desired=reference_energy,
        decimal=5,
    )


def generate_reference():
    """Function to generate reference files, if this script is called directly by a Python interpreter.
    This is ignored by pytest while running tests."""
    reference = dict()
    for mol in [
        "ethanol.sdf",
        "methane_multiconformer.sdf",
        "CID20742535_anion.sdf",
    ]:
        for constrained in [True, False]:
            omm_sys, positions, off_top = _build_system(
                mol=mol,
                constrained=constrained,
            )

            simulation = _build_simulation(omm_sys=omm_sys, off_top=off_top)
            energy = _get_energy(simulation=simulation, positions=positions)

            name = mol + "_"
            if not constrained:
                name += "un"
            name += "constrained"
            reference.update(
                {name: energy.value_in_unit(openmm_unit.kilojoule_per_mole)}
            )

    import openff.toolkit

    toolkit_version = openff.toolkit.__version__

    with open(f"reference_{toolkit_version}.json", "w") as json_out:
        json.dump(reference, json_out)


def _build_system(mol, constrained):
    if constrained:
        parsley = ForceField("openff-1.0.0.offxml")
    else:
        parsley = ForceField("openff_unconstrained-1.0.0.offxml")

    mol = Molecule.from_file(get_data_file_path("molecules/" + mol), file_format="sdf")

    if type(mol) == Molecule:
        off_top = mol.to_topology()
        positions = mol.conformers[0]
    elif type(mol) == list:
        # methane_multiconformer case is a list of two mols
        off_top = Topology()
        for mol_i in mol:
            off_top.add_molecule(mol_i)
        positions = np.vstack([mol[0].conformers[0], mol[1].conformers[0]])

    from openff.toolkit.utils.toolkits import (
        AmberToolsToolkitWrapper,
        RDKitToolkitWrapper,
        ToolkitRegistry,
    )

    toolkit_registry = ToolkitRegistry(
        toolkit_precedence=[RDKitToolkitWrapper, AmberToolsToolkitWrapper]
    )

    omm_sys = parsley.create_openmm_system(off_top, toolkit_registry=toolkit_registry)

    return omm_sys, positions, off_top


def _build_simulation(omm_sys, off_top):
    """Given an OpenMM System, initialize a barebones OpenMM Simulation."""
    # Use OpenMM to compute initial and minimized energy for all conformers
    integrator = openmm.VerletIntegrator(1 * openmm_unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName("Reference")
    omm_top = off_top.to_openmm()
    simulation = openmm.app.Simulation(omm_top, omm_sys, integrator, platform)

    return simulation


def _get_energy(simulation, positions):
    """Given an OpenMM simulation and position, return its energy"""
    simulation.context.setPositions(to_openmm(positions))
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy()

    return energy


if __name__ == "__main__":
    generate_reference()
