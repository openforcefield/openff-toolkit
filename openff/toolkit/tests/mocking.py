import numpy
from openff.units import unit

from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import VirtualSiteHandler


class VirtualSiteMocking:
    @staticmethod
    def sp1_conformer() -> unit.Quantity:
        return unit.Quantity(
            numpy.array(
                [
                    [-2.0, +0.0, +0.0],
                    [-1.0, +0.0, +0.0],
                    [+1.0, +0.0, +0.0],
                    [+2.0, +0.0, +0.0],
                ]
            ),
            unit.angstrom,
        )

    @staticmethod
    def sp2_conformer() -> unit.Quantity:
        return unit.Quantity(
            numpy.array(
                [
                    [+1.0, +0.0, +0.0],
                    [+0.0, +0.0, +0.0],
                    [-1.0, +0.0, -1.0],
                    [-1.0, +0.0, +1.0],
                ]
            ),
            unit.angstrom,
        )

    @staticmethod
    def sp3_conformer() -> unit.Quantity:
        """Returns the conformer of a dummy tetrahedral molecule with the
        first atom positioned at (0, 1, 0), the second atom at (0, 0, 0) and
        the remaining atoms at the corners of the tetrahedron.
        """
        return unit.Quantity(
            numpy.array(
                [
                    [+0.0, +1.0, +0.0],
                    [+0.0, +0.0, +0.0],
                    [-numpy.sqrt(3) / 2.0, -0.5, +0.0],
                    [+numpy.sqrt(3) / 4.0, -0.5, +3.0 / 4.0],
                    [+numpy.sqrt(3) / 4.0, -0.5, -3.0 / 4.0],
                ]
            ),
            unit.angstrom,
        )

    @staticmethod
    def bond_charge_parameter(
        smirks: str, name: str = "EP", param_multiple: float = 1.0
    ) -> VirtualSiteHandler.VirtualSiteType:
        return VirtualSiteHandler.VirtualSiteType(
            type="BondCharge",
            smirks=smirks,
            name=name,
            charge_increment=unit.Quantity(
                [0.1 * param_multiple, 0.2 * param_multiple],
                unit.elementary_charge,
            ),
            sigma=4.0 * param_multiple * unit.angstrom,
            epsilon=3.0 * param_multiple * unit.kilojoule_per_mole,
            match="all_permutations",
            distance=2.0 * param_multiple * unit.angstrom,
        )

    @staticmethod
    def monovalent_parameter(
        smirks: str, name: str = "EP"
    ) -> VirtualSiteHandler.VirtualSiteType:
        return VirtualSiteHandler.VirtualSiteType(
            type="MonovalentLonePair",
            smirks=smirks,
            name=name,
            charge_increment=[0.1, 0.2, 0.3] * unit.elementary_charge,
            sigma=4.0 * unit.angstrom,
            epsilon=5.0 * unit.kilojoule_per_mole,
            match="all_permutations",
            distance=2.0 * unit.angstrom,
            inPlaneAngle=135.0 * unit.degree,
            outOfPlaneAngle=45.0 * unit.degree,
        )

    @staticmethod
    def divalent_parameter(
        smirks: str,
        match: str,
        name: str = "EP",
        angle: unit.Quantity = 0.0 * unit.degree,
    ) -> VirtualSiteHandler.VirtualSiteType:
        return VirtualSiteHandler.VirtualSiteType(
            type="DivalentLonePair",
            smirks=smirks,
            name=name,
            charge_increment=[0.1, 0.2, 0.3] * unit.elementary_charge,
            sigma=4.0 * unit.angstrom,
            epsilon=5.0 * unit.kilojoule_per_mole,
            match=match,
            distance=2.0 * unit.angstrom,
            outOfPlaneAngle=angle,
        )

    @staticmethod
    def trivalent_parameter(
        smirks: str, name: str = "EP"
    ) -> VirtualSiteHandler.VirtualSiteType:
        return VirtualSiteHandler.VirtualSiteType(
            type="TrivalentLonePair",
            smirks=smirks,
            name=name,
            charge_increment=[0.1, 0.2, 0.3, 0.4] * unit.elementary_charge,
            sigma=5.0 * unit.angstrom,
            epsilon=6.0 * unit.kilojoule_per_mole,
            match="once",
            distance=2.0 * unit.angstrom,
        )

    @staticmethod
    def molecule_from_smiles(smiles: str, reverse: bool) -> Molecule:
        molecule = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)

        if reverse:
            mapping = {i: molecule.n_atoms - i - 1 for i in range(molecule.n_atoms)}
            molecule = molecule.remap(mapping)

        return molecule

    @staticmethod
    def chloromethane(reverse: bool = False) -> Molecule:
        return VirtualSiteMocking.molecule_from_smiles(
            "[Cl:1][C:2]([H:3])([H:4])[H:5]", reverse
        )

    @staticmethod
    def formaldehyde(reverse: bool = False) -> Molecule:
        return VirtualSiteMocking.molecule_from_smiles(
            "[O:1]=[C:2]([H:3])[H:4]", reverse
        )

    @staticmethod
    def hypochlorous_acid(reverse: bool = False) -> Molecule:
        return VirtualSiteMocking.molecule_from_smiles("[H:1][O:2][Cl:3]", reverse)

    @staticmethod
    def fake_ammonia(reverse: bool = False) -> Molecule:
        return VirtualSiteMocking.molecule_from_smiles(
            "[N:1]([Cl:2])([Br:3])[H:4]", reverse
        )
