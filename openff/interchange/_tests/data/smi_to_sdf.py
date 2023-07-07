"""Script used to convert molecules.smi to molecules.sdf."""
from openff.toolkit.topology.molecule import Molecule
from rdkit import Chem

mols = Molecule.from_file("molecules.smi", allow_undefined_stereo=True)

with Chem.SDWriter("molecules.sdf") as w:
    for mol in mols:
        try:
            mol.generate_conformers(n_conformers=1)
        except:  # noqa
            continue
        mol.assign_partial_charges(partial_charge_method="am1bcc")
        rdmol = mol.to_rdkit()
        w.write(rdmol)
