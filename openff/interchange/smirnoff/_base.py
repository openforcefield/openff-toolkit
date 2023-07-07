import abc
import json
from typing import Optional, TypeVar, Union

from openff.models.models import DefaultModel
from openff.models.types import custom_quantity_encoder
from openff.toolkit.topology import Topology
from openff.toolkit.typing.engines.smirnoff.parameters import (
    AngleHandler,
    BondHandler,
    ImproperTorsionHandler,
    ParameterHandler,
    ProperTorsionHandler,
)
from openff.units import unit

from openff.interchange.components.potentials import Collection, Potential
from openff.interchange.exceptions import (
    InvalidParameterHandlerError,
    SMIRNOFFParameterAttributeNotImplementedError,
    UnassignedAngleError,
    UnassignedBondError,
    UnassignedTorsionError,
)
from openff.interchange.models import (
    LibraryChargeTopologyKey,
    PotentialKey,
    TopologyKey,
)

T = TypeVar("T", bound="SMIRNOFFCollection")
TP = TypeVar("TP", bound="ParameterHandler")


def _sanitize(o) -> Union[str, dict]:
    # `BaseModel.json()` assumes that all keys and values in dicts are JSON-serializable, which is a problem
    # for the mapping dicts `key_map` and `potentials`.
    if isinstance(o, dict):
        return {_sanitize(k): _sanitize(v) for k, v in o.items()}
    elif isinstance(o, DefaultModel):
        return o.json()
    elif isinstance(o, unit.Quantity):
        return custom_quantity_encoder(o)
    return o


def dump_collection(v, *, default):
    """Dump a SMIRNOFFCollection to JSON after converting to compatible types."""
    return json.dumps(_sanitize(v), default=default)


def collection_loader(data: str) -> dict:
    """Load a JSON blob dumped from a `Collection`."""
    tmp: dict[str, Optional[Union[int, bool, str]]] = {}

    for key, val in json.loads(data).items():
        if isinstance(val, (str, bool, type(None))):
            # These are stored as string but must be parsed into `Quantity`
            if key in ("cutoff", "switch_width"):
                tmp[key] = unit.Quantity(*json.loads(val).values())  # type: ignore[arg-type]
            else:
                tmp[key] = val
        elif isinstance(val, dict):
            if key == "key_map":
                key_map = {}

                for key_, val_ in val.items():
                    if "atom_indices" in key_:
                        topology_key: Union[
                            TopologyKey,
                            LibraryChargeTopologyKey,
                        ] = TopologyKey.parse_raw(key_)

                    else:
                        topology_key = LibraryChargeTopologyKey.parse_raw(key_)

                    potential_key = PotentialKey(**val_)

                    key_map[topology_key] = potential_key

                tmp[key] = key_map  # type: ignore[assignment]

            elif key == "potentials":
                potentials = {}

                for key_, val_ in val.items():
                    potential_key = PotentialKey.parse_raw(key_)
                    potential = Potential.parse_raw(json.dumps(val_))

                    potentials[potential_key] = potential

                tmp[key] = potentials  # type: ignore[assignment]

            else:
                raise NotImplementedError(f"Cannot parse {key} in this JSON.")

    return tmp


# Coped from the toolkit, see
# https://github.com/openforcefield/openff-toolkit/blob/0133414d3ab51e1af0996bcebe0cc1bdddc6431b/
# openff/toolkit/typing/engines/smirnoff/parameters.py#L2318
def _check_all_valence_terms_assigned(
    handler,
    assigned_terms,
    topology,
    valence_terms,
):
    """Check that all valence terms have been assigned."""
    if len(assigned_terms) == len(valence_terms):
        return

    # Convert the valence term to a valence dictionary to make sure
    # the order of atom indices doesn't matter for comparison.
    valence_terms_dict = assigned_terms.__class__()
    for atoms in valence_terms:
        atom_indices = (topology.atom_index(a) for a in atoms)
        valence_terms_dict[atom_indices] = atoms

    # Check that both valence dictionaries have the same keys (i.e. terms).
    assigned_terms_set = set(assigned_terms.keys())
    valence_terms_set = set(valence_terms_dict.keys())
    unassigned_terms = valence_terms_set.difference(assigned_terms_set)
    not_found_terms = assigned_terms_set.difference(valence_terms_set)

    # Raise an error if there are unassigned terms.
    err_msg = ""

    if len(unassigned_terms) > 0:
        unassigned_atom_tuples = []

        unassigned_str = ""
        for unassigned_tuple in unassigned_terms:
            unassigned_str += "\n- Topology indices " + str(unassigned_tuple)
            unassigned_str += ": names and elements "

            unassigned_atoms = []

            # Pull and add additional helpful info on missing terms
            for atom_idx in unassigned_tuple:
                atom = topology.atom(atom_idx)
                unassigned_atoms.append(atom)
                unassigned_str += f"({atom.name} {atom.symbol}), "
            unassigned_atom_tuples.append(tuple(unassigned_atoms))
        err_msg += (
            "{parameter_handler} was not able to find parameters for the following valence terms:\n"
            "{unassigned_str}"
        ).format(
            parameter_handler=handler.__class__.__name__,
            unassigned_str=unassigned_str,
        )
    if len(not_found_terms) > 0:
        if err_msg != "":
            err_msg += "\n"
        not_found_str = "\n- ".join([str(x) for x in not_found_terms])
        err_msg += (
            "{parameter_handler} assigned terms that were not found in the topology:\n"
            "- {not_found_str}"
        ).format(
            parameter_handler=handler.__class__.__name__,
            not_found_str=not_found_str,
        )
    if err_msg:
        err_msg += "\n"

        if isinstance(handler, BondHandler):
            exception_class = UnassignedBondError
        elif isinstance(handler, AngleHandler):
            exception_class = UnassignedAngleError
        elif isinstance(handler, (ProperTorsionHandler, ImproperTorsionHandler)):
            exception_class = UnassignedTorsionError
        else:
            raise RuntimeError(
                f"Could not find an exception class for handler {handler}",
            )

        exception = exception_class(err_msg)
        exception.unassigned_topology_atom_tuples = unassigned_atom_tuples
        exception.handler_class = handler.__class__
        raise exception


class SMIRNOFFCollection(Collection, abc.ABC):
    """Base class for handlers storing potentials produced by SMIRNOFF force fields."""

    is_plugin: bool = False

    def modify_openmm_forces(self, *args, **kwargs):
        """Optionally modify, create, or delete forces. Currently only available to plugins."""
        raise NotImplementedError()

    class Config:
        """Default configuration options for SMIRNOFF potential handlers."""

        json_dumps = dump_collection
        json_loads = collection_loader
        validate_assignment = True
        arbitrary_types_allowed = True

    @classmethod
    @abc.abstractmethod
    def allowed_parameter_handlers(cls):
        """Return a list of allowed types of ParameterHandler classes."""
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def supported_parameters(cls):
        """Return a list of parameter attributes supported by this handler."""
        raise NotImplementedError()

    @classmethod
    def potential_parameters(cls):
        """Return a subset of `supported_parameters` that are meant to be included in potentials."""
        raise NotImplementedError()

    @classmethod
    def check_supported_parameters(cls, parameter_handler: ParameterHandler):
        """Verify that a parameter handler is in an allowed list of handlers."""
        for parameter in parameter_handler.parameters:
            for parameter_attribute in parameter._get_defined_parameter_attributes():
                if parameter_attribute == "parent_id":
                    continue
                if parameter_attribute not in cls.supported_parameters():
                    raise SMIRNOFFParameterAttributeNotImplementedError(
                        parameter_attribute,
                    )

    @classmethod
    def check_openmm_requirements(cls, combine_nonbonded_forces: bool) -> None:
        """Run through a list of assertions about what is compatible when exporting this to OpenMM."""

    def store_matches(
        self,
        parameter_handler: ParameterHandler,
        topology: "Topology",
    ) -> None:
        """Populate self.key_map with key-val pairs of [TopologyKey, PotentialKey]."""
        if self.key_map:
            # TODO: Should the key_map always be reset, or should we be able to partially
            # update it? Also Note the duplicated code in the child classes
            self.key_map: dict[
                Union[TopologyKey, LibraryChargeTopologyKey],
                PotentialKey,
            ] = dict()
        matches = parameter_handler.find_matches(topology)
        for key, val in matches.items():
            topology_key = TopologyKey(atom_indices=key)
            potential_key = PotentialKey(
                id=val.parameter_type.smirks,
                associated_handler=parameter_handler.TAGNAME,
            )
            self.key_map[topology_key] = potential_key

        if self.__class__.__name__ in [
            "SMIRNOFFBondCollection",
            "SMIRNOFFAngleCollection",
        ]:
            valence_terms = self.valence_terms(topology)

            _check_all_valence_terms_assigned(
                handler=parameter_handler,
                assigned_terms=matches,
                topology=topology,
                valence_terms=valence_terms,
            )

    def store_potentials(self, parameter_handler: TP):
        """
        Populate self.potentials with key-val pairs of [PotentialKey, Potential].
        """
        raise NotImplementedError()

    @classmethod
    def create(
        cls: type[T],
        parameter_handler: TP,
        topology: "Topology",
    ) -> T:
        """
        Create a SMIRNOFFCOllection from toolkit data.

        """
        if type(parameter_handler) not in cls.allowed_parameter_handlers():
            raise InvalidParameterHandlerError(type(parameter_handler))

        handler = cls()
        if hasattr(handler, "fractional_bondorder_method"):
            if getattr(parameter_handler, "fractional_bondorder_method", None):
                handler.fractional_bond_order_method = (  # type: ignore[attr-defined]
                    parameter_handler.fractional_bondorder_method
                )
                handler.fractional_bond_order_interpolation = (  # type: ignore[attr-defined]
                    parameter_handler.fractional_bondorder_interpolation
                )
        handler.store_matches(parameter_handler=parameter_handler, topology=topology)
        handler.store_potentials(parameter_handler=parameter_handler)

        return handler

    def __repr__(self) -> str:
        return (
            f"Handler '{self.type}' with expression '{self.expression}', {len(self.key_map)} slots, "
            f"and {len(self.potentials)} potentials"
        )
