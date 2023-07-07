import networkx
from openff.toolkit import Topology
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange._experimental import experimental
from openff.interchange.common._nonbonded import ElectrostaticsCollection, vdWCollection
from openff.interchange.common._valence import (
    AngleCollection,
    BondCollection,
    ConstraintCollection,
    ImproperTorsionCollection,
    ProperTorsionCollection,
    RyckaertBellemansTorsionCollection,
)
from openff.interchange.components.potentials import Potential
from openff.interchange.interop.gromacs.models.models import (
    GROMACSSystem,
    PeriodicImproperDihedral,
    PeriodicProperDihedral,
    RyckaertBellemansDihedral,
)
from openff.interchange.models import (
    AngleKey,
    BondKey,
    ImproperTorsionKey,
    PotentialKey,
    ProperTorsionKey,
    TopologyKey,
)


@experimental
def to_interchange(
    system: GROMACSSystem,
) -> Interchange:
    """
    Convert a GROMACS system to an Interchange object.

    Parameters
    ----------
    system: GROMACSSystem
        The GROMACS system to convert.

    Returns
    -------
    interchange
        The converted Interchange object.

    """
    if (system.nonbonded_function, system.gen_pairs) != (1, True):
        raise NotImplementedError()

    if system.combination_rule == 2:
        mixing_rule = "lorentz-berthelot"
    elif system.combination_rule == 3:
        mixing_rule = "geometric"

    vdw = vdWCollection(scale14=system.vdw_14, mixing_rule=mixing_rule)
    electrostatics = ElectrostaticsCollection(version=0.4, scale_14=system.coul_14)

    constraints = ConstraintCollection()
    bonds = BondCollection()
    angles = AngleCollection()
    # TODO: Split out periodic and RB/other styles?
    periodic_propers = ProperTorsionCollection()
    rb_torsions = RyckaertBellemansTorsionCollection()
    impropers = ImproperTorsionCollection()

    vdw.potentials = {
        PotentialKey(id=f"{atom_type.name}"): Potential(
            parameters={"sigma": atom_type.sigma, "epsilon": atom_type.epsilon},
        )
        for atom_type in system.atom_types.values()
    }

    molecule_start_index = 0

    for molecule_name, molecule_type in system.molecule_types.items():
        for _ in range(system.molecules[molecule_name]):
            for atom in molecule_type.atoms:
                topology_atom_index = molecule_start_index + atom.index - 1
                topology_key = TopologyKey(
                    atom_indices=(topology_atom_index,),
                )

                vdw.key_map.update(
                    {topology_key: PotentialKey(id=f"{atom.atom_type}")},
                )

                # GROMACS does NOT necessarily tie partial charges to atom types, so need a new key for each atom
                electrostatics_key = PotentialKey(id=f"{topology_key.atom_indices[0]}")
                electrostatics.key_map.update(
                    {topology_key: electrostatics_key},
                )

                electrostatics.potentials.update(
                    {electrostatics_key: Potential(parameters={"charge": atom.charge})},
                )

            if len(molecule_type.settles) > 0:
                if [
                    system.atom_types[a.atom_type].atomic_number
                    for a in molecule_type.atoms
                ] != [8, 1, 1]:
                    raise NotImplementedError(
                        "Settles have only been implemented for water with OHH ordering.",
                    )

            for settle in molecule_type.settles:
                oxygen_hydrogen1 = BondKey(
                    atom_indices=(
                        1 + molecule_start_index - 1,
                        2 + molecule_start_index - 1,
                    ),
                )
                oxygen_hydrogen2 = BondKey(
                    atom_indices=(
                        1 + molecule_start_index - 1,
                        3 + molecule_start_index - 1,
                    ),
                )
                hydrogen_hydrogen = BondKey(
                    atom_indices=(
                        2 + molecule_start_index - 1,
                        3 + molecule_start_index - 1,
                    ),
                )

                # Could probably look up the key and potential instead of re-creating?
                oxygen_hydrogen_key = PotentialKey(
                    id="O-H-settles",
                    associated_handler="Constraints",
                )

                hydrogen_hydrogen_key = PotentialKey(
                    id="H-H-settles",
                    associated_handler="Constraints",
                )

                oxygen_hydrogen_constraint = Potential(
                    parameters={
                        "distance": settle.oxygen_hydrogen_distance,
                    },
                )
                hydrogen_hydrogen_constraint = Potential(
                    parameters={
                        "distance": settle.hydrogen_hydrogen_distance,
                    },
                )

                constraints.key_map.update({oxygen_hydrogen1: oxygen_hydrogen_key})
                constraints.key_map.update({oxygen_hydrogen2: oxygen_hydrogen_key})
                constraints.key_map.update({hydrogen_hydrogen: hydrogen_hydrogen_key})

                constraints.potentials.update(
                    {oxygen_hydrogen_key: oxygen_hydrogen_constraint},
                )
                constraints.potentials.update(
                    {hydrogen_hydrogen_key: hydrogen_hydrogen_constraint},
                )

            # TODO: Build constraints from settles
            for bond in molecule_type.bonds:
                topology_key = BondKey(
                    atom_indices=(
                        bond.atom1 + molecule_start_index - 1,
                        bond.atom2 + molecule_start_index - 1,
                    ),
                )

                potential_key = PotentialKey(
                    id="-".join(map(str, topology_key.atom_indices)),
                )

                potential = Potential(
                    parameters={
                        "k": bond.k,
                        "length": bond.length,
                    },
                )

                bonds.key_map.update({topology_key: potential_key})
                bonds.potentials.update({potential_key: potential})

            for angle in molecule_type.angles:
                topology_key = AngleKey(
                    atom_indices=(
                        angle.atom1 + molecule_start_index - 1,
                        angle.atom2 + molecule_start_index - 1,
                        angle.atom3 + molecule_start_index - 1,
                    ),
                )

                potential_key = PotentialKey(
                    id="-".join(map(str, topology_key.atom_indices)),
                )

                potential = Potential(
                    parameters={
                        "k": angle.k,
                        "angle": angle.angle,
                    },
                )

                angles.key_map.update({topology_key: potential_key})
                angles.potentials.update({potential_key: potential})

            for dihedral in molecule_type.dihedrals:
                if isinstance(dihedral, PeriodicProperDihedral):
                    key_type = ProperTorsionKey
                    collection = periodic_propers
                elif isinstance(dihedral, RyckaertBellemansDihedral):
                    key_type = ProperTorsionKey
                    collection = rb_torsions  # type: ignore[assignment]
                elif isinstance(dihedral, PeriodicImproperDihedral):
                    key_type = ImproperTorsionKey
                    collection = impropers  # type: ignore[assignment]
                else:
                    raise NotImplementedError(
                        f"Dihedral type {type(dihedral)} not implemented.",
                    )

                topology_key = key_type(
                    atom_indices=(
                        dihedral.atom1 + molecule_start_index - 1,
                        dihedral.atom2 + molecule_start_index - 1,
                        dihedral.atom3 + molecule_start_index - 1,
                        dihedral.atom4 + molecule_start_index - 1,
                    ),
                )

                potential_key = PotentialKey(
                    id="-".join(map(str, topology_key.atom_indices)),
                )

                if isinstance(
                    dihedral,
                    (PeriodicProperDihedral, PeriodicImproperDihedral),
                ):
                    potential = Potential(
                        parameters={
                            "periodicity": unit.Quantity(
                                dihedral.multiplicity,
                                unit.dimensionless,
                            ),
                            "phase": dihedral.phi,
                            "k": dihedral.k,
                            "idivf": 1 * unit.dimensionless,
                        },
                    )

                elif isinstance(dihedral, RyckaertBellemansDihedral):
                    potential = Potential(
                        parameters={
                            "c0": dihedral.c0,
                            "c1": dihedral.c1,
                            "c2": dihedral.c2,
                            "c3": dihedral.c3,
                            "c4": dihedral.c4,
                            "c5": dihedral.c5,
                        },
                    )

                else:
                    raise NotImplementedError()

                collection.key_map.update({topology_key: potential_key})
                collection.potentials.update({potential_key: potential})

            molecule_start_index += len(molecule_type.atoms)

    interchange = Interchange()

    interchange.collections.update(
        {
            "vdW": vdw,
            "Electrostatics": electrostatics,
            "Constraints": constraints,
            "Bonds": bonds,
            "Angles": angles,
            "ProperTorsions": periodic_propers,
            "RBTorsions": rb_torsions,
            "ImproperTorsions": impropers,
        },
    )

    interchange.topology = _convert_topology(system)
    interchange.positions = system.positions
    interchange.box = system.box

    return interchange


def _convert_topology(
    system: GROMACSSystem,
) -> Topology:
    from openff.toolkit.topology._mm_molecule import _SimpleMolecule

    topology = Topology()

    for molecule_name, molecule_type in system.molecule_types.items():
        graph = networkx.Graph()

        n_molecules = system.molecules[molecule_name]

        for atom in molecule_type.atoms:
            graph.add_node(
                atom.index - 1,
                atomic_number=system.atom_types[atom.atom_type].atomic_number,
            )

        if [
            system.atom_types[atom.atom_type].atomic_number
            for atom in molecule_type.atoms
        ] == [8, 1, 1]:
            graph.add_edge(0, 1)
            graph.add_edge(0, 2)
        else:
            for bond in molecule_type.bonds:
                graph.add_edge(
                    bond.atom1 - 1,
                    bond.atom2 - 1,
                )

        molecule = _SimpleMolecule._from_subgraph(graph)
        molecule.name = molecule_name

        for _ in range(n_molecules):
            topology.add_molecule(molecule)

    return topology
