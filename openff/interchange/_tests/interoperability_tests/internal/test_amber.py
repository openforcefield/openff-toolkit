import MDAnalysis
import mdtraj
import numpy as np
import parmed
import pytest
from openff.toolkit import ForceField, Molecule
from openff.toolkit.tests.utils import get_data_file_path
from openff.units import unit
from openmm import app

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest
from openff.interchange.drivers import get_amber_energies, get_openmm_energies


class TestAmber(_BaseTest):
    @pytest.mark.skip(reason="Need replacement route to reference positions")
    def test_inpcrd(self, sage):
        mol = Molecule.from_smiles(10 * "C")
        mol.name = "HPER"
        mol.generate_conformers(n_conformers=1)

        out = Interchange.from_smirnoff(force_field=sage, topology=mol.to_topology())
        out.box = [4, 4, 4]
        out.positions = mol.conformers[0]
        out.positions = unit.nanometer * np.round(out.positions.m_as(unit.nanometer), 5)

        out.to_inpcrd("internal.inpcrd")
        # This method no longer exists
        out._to_parmed().save("parmed.inpcrd")

        coords1 = parmed.load_file("internal.inpcrd").coordinates
        coords2 = parmed.load_file("parmed.inpcrd").coordinates

        np.testing.assert_equal(coords1, coords2)

    @pytest.mark.slow()
    @pytest.mark.parametrize(
        "smiles",
        [
            "C",
            "CC",  # Adds a proper torsion term(s)
            "C=O",  # Simplest molecule with any improper torsion
            "OC=O",  # Simplest molecule with a multi-term torsion
            "CCOC",  # This hits t86, which has a non-1.0 idivf
            "C1COC(=O)O1",  # This adds an improper, i2
        ],
    )
    def test_amber_energy(self, sage_unconstrained, smiles):
        """
        Basic test to see if the amber energy driver is functional.

        Note this test can only use the unconstrained version of Sage because sander applies SHAKE
        constraints in the single-point energy calculation, i.e. uses geometries with constraints
        applied, NOT what is in the coordinate file. See issue
        https://github.com/openforcefield/openff-interchange/issues/323
        """
        mol = Molecule.from_smiles(smiles)
        mol.generate_conformers(n_conformers=1)
        top = mol.to_topology()

        interchange = Interchange.from_smirnoff(sage_unconstrained, top)

        interchange.box = [5, 5, 5]
        interchange.positions = mol.conformers[0]

        omm_energies = get_openmm_energies(interchange, combine_nonbonded_forces=True)
        amb_energies = get_amber_energies(interchange)

        # TODO: More investigation into possible non-bonded energy differences and better reporting.
        #       03/02/2023 manually inspected some files and charges and vdW parameters are
        #       precisely identical. Passing box vectors to prmtop files might not always work.
        omm_energies.energies.pop("Nonbonded")
        amb_energies.energies.pop("vdW")
        amb_energies.energies.pop("Electrostatics")

        omm_energies.compare(amb_energies)


class TestPRMTOP(_BaseTest):
    def test_atom_names_pdb(self):
        peptide = Molecule.from_polymer_pdb(
            get_data_file_path("proteins/MainChain_ALA_ALA.pdb"),
        )
        ff14sb = ForceField("ff14sb_off_impropers_0.0.3.offxml")

        Interchange.from_smirnoff(ff14sb, peptide.to_topology()).to_prmtop(
            "atom_names.prmtop",
        )

        pdb_object = app.PDBFile(get_data_file_path("proteins/MainChain_ALA_ALA.pdb"))
        openmm_object = app.AmberPrmtopFile("atom_names.prmtop")
        mdanalysis_object = MDAnalysis.Universe("atom_names.prmtop")
        # This may not be useful, see
        # https://github.com/mdtraj/mdtraj/blob/6bb35a8d78a5758ff1f72b3af1fc21d2e38f1029/mdtraj/formats/prmtop.py#L47-L49
        mdtraj_object = mdtraj.load_prmtop("atom_names.prmtop")

        pdb_atom_names = [atom.name for atom in pdb_object.topology.atoms()]

        openmm_atom_names = [atom.name for atom in openmm_object.topology.atoms()]
        mdanalysis_atom_names = [atom.name for atom in mdanalysis_object.atoms]
        mdtraj_atom_names = [atom.name for atom in mdtraj_object.atoms]

        assert openmm_atom_names == pdb_atom_names
        assert mdanalysis_atom_names == pdb_atom_names
        assert mdtraj_atom_names == pdb_atom_names


class TestAmberResidues(_BaseTest):
    @pytest.mark.parametrize("patch_residue_name", [True, False])
    def test_single_residue_system_residue_name(
        self,
        tmp_path,
        sage,
        ethanol,
        patch_residue_name,
    ):
        if patch_residue_name:
            for atom in ethanol.atoms:
                atom.metadata["residue_name"] = "YUP"

            ethanol.add_default_hierarchy_schemes()

        Interchange.from_smirnoff(sage, [ethanol]).to_prmtop("test.prmtop")

        residue_names = [r.name for r in parmed.load_file("test.prmtop").residues]

        if patch_residue_name:
            assert residue_names == ["YUP"]
        else:
            assert residue_names == ["RES"]
