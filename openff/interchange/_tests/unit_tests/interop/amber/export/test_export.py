import openmm
import openmm.app
import parmed
import pytest
from openff.toolkit import ForceField, Molecule
from openff.toolkit.tests.utils import requires_openeye

from openff.interchange import Interchange
from openff.interchange._tests import get_test_file_path


@requires_openeye
@pytest.mark.parametrize(
    "molecule",
    [
        "lig_CHEMBL3265016-1.pdb",
        "c1ccc2ccccc2c1",
    ],
)
def test_atom_names_with_padding(molecule):
    # pytest processes fixtures before the decorator can be applied
    if molecule.endswith(".pdb"):
        molecule = Molecule(get_test_file_path(molecule).as_posix())
    else:
        molecule = Molecule.from_smiles(molecule)

    # Unclear if the toolkit will always load PDBs with padded whitespace in name
    Interchange.from_smirnoff(
        ForceField("openff-2.0.0.offxml"),
        molecule.to_topology(),
    ).to_prmtop("tmp.prmtop")

    # Loading with ParmEd striggers #679 if exclusions lists are wrong
    parmed.load_file("tmp.prmtop")


@pytest.mark.parametrize("molecule", ["C1=CN=CN1", "c1ccccc1", "c1ccc2ccccc2c1"])
def exclusions_in_rings(molecule):
    molecule.generate_conformers(n_conformers=1)
    topology = molecule.to_topology()
    topology.box_vectors = [4, 4, 4]

    sage_no_impropers = ForceField("openff-2.0.0.offxml")
    sage_no_impropers.deregister_parameter_handler("ImproperTorsions")

    interchange = sage_no_impropers.create_interchange(topology)

    interchange.to_prmtop("tmp.prmtop")

    # Use the OpenMM export as source of truth
    openmm_system = interchange.to_openmm()
    for force in openmm_system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            reference = force.getNumExceptions()

    loaded_system = openmm.app.AmberPrmtopFile("tmp.prmtop").createSystem()
    for force in loaded_system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            assert force.getNumExceptions() == reference
