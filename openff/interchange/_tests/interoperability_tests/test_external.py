import numpy as np
import pytest
from openff.toolkit.topology import Molecule, Topology
from openff.units import unit
from openff.units.openmm import from_openmm
from openmm import app
from openmm import unit as openmm_unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest, get_test_file_path


class TestFromOpenMM(_BaseTest):
    @pytest.mark.slow()
    def test_from_openmm_pdbfile(self, argon_ff):
        pdb_file_path = get_test_file_path("10-argons.pdb").as_posix()
        pdbfile = app.PDBFile(pdb_file_path)

        mol = Molecule.from_smiles("[#18]")
        tmp = Topology.from_openmm(pdbfile.topology, unique_molecules=[mol])

        out = Interchange.from_smirnoff(argon_ff, tmp)
        out.box = from_openmm(pdbfile.topology.getPeriodicBoxVectors())
        out.positions = from_openmm(pdbfile.getPositions())

        assert np.allclose(
            out.positions.to(unit.nanometer).magnitude,
            pdbfile.getPositions().value_in_unit(openmm_unit.nanometer),
        )
