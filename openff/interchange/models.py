"""Custom Pydantic models."""
import abc
from typing import Literal, Optional

from openff.models.models import DefaultModel
from pydantic import Field


class TopologyKey(DefaultModel, abc.ABC):
    """
    A unique identifier of a segment of a chemical topology.

    These refer to a single portion of a chemical graph, i.e. a single valence term,
    (a bond, angle, or dihedral) or a single atom. These target only the information in
    the chemical graph and do not store physics parameters. For example, a TopologyKey
    corresponding to a bond would store the indices of the two atoms that compose the
    bond, but not the force constant or equilibrium bond length as determined by the
    force field.

    Examples
    --------
    Create a TopologyKey identifying some speicfic angle

    .. code-block:: pycon

        >>> from openff.interchange.models import TopologyKey
        >>> this_angle = TopologyKey(atom_indices=(2, 1, 3))
        >>> this_angle
        TopologyKey with atom indices (2, 1, 3)

    Create a TopologyKey indentifying just one atom

    .. code-block:: pycon

        >>> this_atom = TopologyKey(atom_indices=(4,))
        >>> this_atom
        TopologyKey with atom indices (4,)

    """

    # TODO: Swith to `pydantic.contuple` once 1.10.3 or 2.0.0 is released
    atom_indices: tuple[int, ...] = Field(
        description="The indices of the atoms occupied by this interaction",
    )

    def __hash__(self) -> int:
        return hash(tuple(self.atom_indices))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__} with atom indices {self.atom_indices}"


class BondKey(TopologyKey):
    """
    A unique identifier of the atoms associated in a bond potential.
    """

    atom_indices: tuple[int, ...] = Field(
        description="The indices of the atoms occupied by this interaction",
    )

    bond_order: Optional[float] = Field(
        None,
        description=(
            "If this key represents as topology component subject to interpolation between "
            "multiple parameters(s), the bond order determining the coefficients of the wrapped "
            "potentials."
        ),
    )

    def __hash__(self) -> int:
        return hash((tuple(self.atom_indices), self.bond_order))

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__} with atom indices {self.atom_indices}"
            f"{'' if self.bond_order is None else ', bond order ' + str(self.bond_order)}"
        )


class AngleKey(TopologyKey):
    """
    A unique identifier of the atoms associated in an angle potential.
    """

    atom_indices: tuple[int, ...] = Field(
        description="The indices of the atoms occupied by this interaction",
    )


class ProperTorsionKey(TopologyKey):
    """
    A unique identifier of the atoms associated in a proper torsion potential.
    """

    atom_indices: tuple[int, ...] = Field(
        description="The indices of the atoms occupied by this interaction",
    )

    mult: Optional[int] = Field(
        None,
        description="The index of this duplicate interaction",
    )

    phase: Optional[float] = Field(
        None,
        description="If this key represents as topology component subject to interpolation between "
        "multiple parameters(s), the phase determining the coefficients of the wrapped "
        "potentials.",
    )

    bond_order: Optional[float] = Field(
        None,
        description=(
            "If this key represents as topology component subject to interpolation between "
            "multiple parameters(s), the bond order determining the coefficients of the wrapped "
            "potentials."
        ),
    )

    def __hash__(self) -> int:
        return hash((tuple(self.atom_indices), self.mult, self.bond_order, self.phase))

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__} with atom indices {self.atom_indices}"
            f"{'' if self.mult is None else ', mult ' + str(self.mult)}"
            f"{'' if self.bond_order is None else ', bond order ' + str(self.bond_order)}"
        )


class ImproperTorsionKey(ProperTorsionKey):
    """
    A unique identifier of the atoms associated in an improper torsion potential.

    The central atom is the second atom in the `atom_indices` tuple, or accessible via `get_central_atom_index`.
    """

    def get_central_atom_index(self) -> int:
        """Get the index of the central atom of this improper torsion."""
        return self.atom_indices[1]


class LibraryChargeTopologyKey(DefaultModel):
    """
    A unique identifier of the atoms associated with a library charge.
    """

    # TODO: Store all atoms associated with this charge?
    # TODO: Is there an upper bound on the number of atoms that can be associated with a LibraryChargeType?
    # TODO: Eventually rename this for coherence with `TopologyKey`
    this_atom_index: int

    @property
    def atom_indices(self) -> tuple[int, ...]:
        """Alias for `this_atom_index`."""
        return (self.this_atom_index,)

    def __hash__(self) -> int:
        return hash((self.this_atom_index,))


class SingleAtomChargeTopologyKey(LibraryChargeTopologyKey):
    """
    Shim class for storing the result of charge_from_molecules.
    """


class ChargeModelTopologyKey(DefaultModel):
    """Subclass of `TopologyKey` for use with charge models only."""

    this_atom_index: int
    partial_charge_method: str

    @property
    def atom_indices(self) -> tuple[int, ...]:
        """Alias for `this_atom_index`."""
        return (self.this_atom_index,)

    def __hash__(self) -> int:
        return hash((self.this_atom_index, self.partial_charge_method))


class ChargeIncrementTopologyKey(DefaultModel):
    """Subclass of `TopologyKey` for use with charge increments only."""

    # TODO: Eventually rename this for coherence with `TopologyKey`
    this_atom_index: int
    other_atom_indices: tuple[int, ...]

    @property
    def atom_indices(self) -> tuple[int, ...]:
        """Alias for `this_atom_index`."""
        return (self.this_atom_index,)

    def __hash__(self) -> int:
        return hash((self.this_atom_index, self.other_atom_indices))


class VirtualSiteKey(TopologyKey):
    """A unique identifier of a virtual site in the scope of a chemical topology."""

    # TODO: Overriding the attribute of a parent class is clumsy, but less grief than
    #       having this not inherit from `TopologyKey`. It might be useful to just have
    #       orientation_atom_indices point to the same thing.
    atom_indices: Optional[tuple[int]] = None  # type: ignore[assignment]

    orientation_atom_indices: tuple[int, ...] = Field(
        description="The indices of the 'orientation atoms' which are used to define the position "
        "of this virtual site. The first atom is the 'parent atom' which defines which atom the "
        "virtual site is 'attached' to.",
    )
    type: str = Field(description="The type of this virtual site parameter.")
    name: str = Field(description="The name of this virtual site parameter.")
    match: Literal["once", "all_permutations"] = Field(
        description="The `match` attribute of the associated virtual site type",
    )

    def __hash__(self) -> int:
        return hash(
            (
                self.orientation_atom_indices,
                self.name,
                self.type,
                self.match,
            ),
        )


class PotentialKey(DefaultModel):
    """
    A unique identifier of an instance of physical parameters as applied to a segment of a chemical topology.

    These refer to a single term in a force field as applied to a single segment of a chemical
    topology, i.e. a single atom or dihedral. For example, a PotentialKey corresponding to a
    bond would store the the force constant and the equilibrium bond length as determined by
    the force field. These keys to not have direct knowledge of where in a topology they have been
    applied.

    Examples
    --------
    Create a PotentialKey corresponding to the parameter with id `b55` in OpenFF "Parsley" 1.0.0

    .. code-block:: pycon

        >>> from openff.interchange.models import PotentialKey
        >>> from openff.toolkit.typing.engines.smirnoff import ForceField
        >>> parsley = ForceField("openff-1.0.0.offxml")
        >>> param = parsley["Bonds"].get_parameter({"id": "b55"})[0]
        >>> bond_55 = PotentialKey(id=param.smirks)
        >>> bond_55
        PotentialKey associated with handler 'None' with id '[#16X4,#16X3:1]-[#8X2:2]'

    Create a PotentialKey corresponding to the angle parameters in OPLS-AA defined
    between atom types opls_135, opls_135, and opls_140

    .. code-block:: pycon

        >>> oplsaa_angle = PotentialKey(id="opls_135-opls_135-opls_140")
        >>> oplsaa_angle
        PotentialKey associated with handler 'None' with id 'opls_135-opls_135-opls_140'

    """

    id: str = Field(
        ...,
        description="A unique identifier of this potential, i.e. a SMARTS pattern or an atom type",
    )
    mult: Optional[int] = Field(
        None,
        description="The index of this duplicate interaction",
    )
    associated_handler: Optional[str] = Field(
        None,
        description="The type of handler this potential key is associated with, "
        "i.e. 'Bonds', 'vdW', or 'LibraryCharges",
    )
    bond_order: Optional[float] = Field(
        None,
        description="If this is a key to a WrappedPotential interpolating multiple parameter(s), "
        "the bond order determining the coefficients of the wrapped potentials.",
    )

    def __hash__(self) -> int:
        return hash((self.id, self.mult, self.associated_handler, self.bond_order))

    def __repr__(self) -> str:
        return (
            f"PotentialKey associated with handler '{self.associated_handler}' with id '{self.id}'"
            f"{'' if self.mult is None else ', mult ' + str(self.mult)}"
            f"{'' if self.bond_order is None else ', bond order ' + str(self.bond_order)}"
        )
