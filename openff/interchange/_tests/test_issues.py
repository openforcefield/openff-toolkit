"""Tests reproducing specific issues that are otherwise uncategorized."""
import parmed
from openff.toolkit import ForceField, Molecule


def test_issue_723():
    force_field = ForceField("openff-2.1.0.offxml")

    molecule = Molecule.from_smiles("C#N")
    molecule.generate_conformers(n_conformers=1)

    force_field.create_interchange(molecule.to_topology()).to_top("_x.top")

    parmed.load_file("_x.top")
