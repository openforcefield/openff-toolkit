"""
Test the behavior of the drivers.all module
"""
import math
from shutil import which

import pandas
import pytest
from openff.toolkit.topology import Molecule

from openff.interchange import Interchange
from openff.interchange._tests import (
    _BaseTest,
    needs_gmx,
    needs_lmp,
    needs_not_gmx,
    needs_not_lmp,
)
from openff.interchange.drivers.all import get_all_energies, get_summary_data
from openff.interchange.drivers.gromacs import _find_gromacs_executable
from openff.interchange.drivers.lammps import _find_lammps_executable


@pytest.mark.slow()
class TestDriversAll(_BaseTest):
    @pytest.fixture()
    def basic_interchange(self, sage_unconstrained):
        molecule = Molecule.from_smiles("CCO")
        molecule.generate_conformers(n_conformers=1)
        molecule.name = "MOL"
        topology = molecule.to_topology()

        out = Interchange.from_smirnoff(sage_unconstrained, topology)
        out.positions = molecule.conformers[0]
        out.box = [4, 4, 4]

        return out

    @needs_gmx
    @needs_lmp
    def test_all_with_all(self, basic_interchange):
        summary = get_all_energies(basic_interchange)

        assert len(summary) == 4

    @needs_not_gmx
    @needs_not_lmp
    def test_all_with_minimum(self, basic_interchange):
        summary = get_all_energies(basic_interchange)

        assert len(summary) == 2

    def test_skipping(self, basic_interchange):
        summary = get_all_energies(basic_interchange)

        assert ("GROMACS" in summary) == (_find_gromacs_executable() is not None)
        assert ("Amber" in summary) == (which("sander") is not None)
        assert ("LAMMPS" in summary) == (_find_lammps_executable() is not None)

    # TODO: Also run all of this with h-bond constraints
    def test_summary_data(self, basic_interchange):
        summary = get_summary_data(basic_interchange)

        assert isinstance(summary, pandas.DataFrame)

        assert "OpenMM" in summary.index

        assert ("GROMACS" in summary.index) == (_find_gromacs_executable() is not None)
        assert ("Amber" in summary.index) == (which("sander") is not None)
        assert ("LAMMPS" in summary.index) == (_find_lammps_executable() is not None)

        # Check that (some of) the data is reasonable, this tolerance should be greatly reduced
        # See https://github.com/openforcefield/openff-interchange/issues/632
        for key in ["Bond", "Angle", "Torsion"]:
            assert summary.describe().loc["std", key] < 0.001

        # Check that (some of) the data did not NaN out
        for val in summary["Torsion"].to_dict().values():
            assert not math.isnan(val)
