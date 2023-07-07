"""Models for storing applied force field parameters."""
import ast
import json
import warnings
from typing import Callable, Optional, Union

import numpy
from openff.models.models import DefaultModel
from openff.models.types import ArrayQuantity, FloatQuantity
from openff.units import unit
from openff.utilities.utilities import has_package, requires_package
from pydantic import Field, PrivateAttr, validator

from openff.interchange.exceptions import MissingParametersError
from openff.interchange.models import (
    LibraryChargeTopologyKey,
    PotentialKey,
    TopologyKey,
)
from openff.interchange.warnings import InterchangeDeprecationWarning

if has_package("jax"):
    from jax import numpy as jax_numpy

from numpy.typing import ArrayLike

if has_package("jax"):
    from jax import Array


def __getattr__(name: str):
    if name == "PotentialHandler":
        warnings.warn(
            "`PotentialHandler` has been renamed to `Collection`. "
            "Importing `Collection` instead.",
            InterchangeDeprecationWarning,
        )
        return Collection

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def potential_loader(data: str) -> dict:
    """Load a JSON blob dumped from a `Collection`."""
    tmp: dict[str, Union[int, bool, str, dict]] = {}

    for key, val in json.loads(data).items():
        if isinstance(val, (str, type(None))):
            tmp[key] = val  # type: ignore
        elif isinstance(val, dict):
            if key == "parameters":
                tmp["parameters"] = dict()

                for key_, val_ in val.items():
                    loaded = json.loads(val_)
                    tmp["parameters"][key_] = unit.Quantity(  # type: ignore[index]
                        loaded["val"],
                        loaded["unit"],
                    )

    return tmp


class Potential(DefaultModel):
    """Base class for storing applied parameters."""

    parameters: dict[str, FloatQuantity] = dict()
    map_key: Optional[int] = None

    class Config:
        """Pydantic configuration."""

        json_encoders: dict[type, Callable] = DefaultModel.Config.json_encoders
        json_loads: Callable = potential_loader
        validate_assignment: bool = True
        arbitrary_types_allowed: bool = True

    @validator("parameters")
    def validate_parameters(
        cls,
        v: dict[str, Union[ArrayQuantity, FloatQuantity]],
    ) -> dict[str, FloatQuantity]:
        for key, val in v.items():
            if isinstance(val, list):
                v[key] = ArrayQuantity.validate_type(val)
            else:
                v[key] = FloatQuantity.validate_type(val)
        return v

    def __hash__(self) -> int:
        return hash(tuple(self.parameters.values()))


class WrappedPotential(DefaultModel):
    """Model storing other Potential model(s) inside inner data."""

    class InnerData(DefaultModel):
        """The potentials being wrapped."""

        data: dict[Potential, float]

    _inner_data: InnerData = PrivateAttr()

    def __init__(self, data: Union[Potential, dict]) -> None:
        if isinstance(data, Potential):
            self._inner_data = self.InnerData(data={data: 1.0})
        elif isinstance(data, dict):
            self._inner_data = self.InnerData(data=data)

    @property
    def parameters(self) -> dict[str, FloatQuantity]:
        """Get the parameters as represented by the stored potentials and coefficients."""
        keys: set[str] = {
            param_key
            for pot in self._inner_data.data.keys()
            for param_key in pot.parameters.keys()
        }
        params = dict()
        for key in keys:
            params.update(
                {
                    key: sum(
                        coeff * pot.parameters[key]
                        for pot, coeff in self._inner_data.data.items()
                    ),
                },
            )
        return params

    def __repr__(self) -> str:
        return str(self._inner_data.data)


class Collection(DefaultModel):
    """Base class for storing parametrized force field data."""

    type: str = Field(..., description="The type of potentials this handler stores.")
    is_plugin: bool = Field(
        False,
        description="Whether this collection is defined as a plugin.",
    )
    expression: str = Field(
        ...,
        description="The analytical expression governing the potentials in this handler.",
    )
    key_map: dict[Union[TopologyKey, LibraryChargeTopologyKey], PotentialKey] = Field(
        dict(),
        description="A mapping between TopologyKey objects and PotentialKey objects.",
    )
    potentials: dict[PotentialKey, Union[Potential, WrappedPotential]] = Field(
        dict(),
        description="A mapping between PotentialKey objects and Potential objects.",
    )

    @property
    def independent_variables(self) -> set[str]:
        """
        Return a set of variables found in the expression but not in any potentials.
        """
        vars_in_potentials = set([*self.potentials.values()][0].parameters.keys())
        vars_in_expression = {
            node.id
            for node in ast.walk(ast.parse(self.expression))
            if isinstance(node, ast.Name)
        }
        return vars_in_expression - vars_in_potentials

    def _get_parameters(self, atom_indices: tuple[int]) -> dict:
        for topology_key in self.key_map:
            if topology_key.atom_indices == atom_indices:
                potential_key = self.key_map[topology_key]
                potential = self.potentials[potential_key]
                parameters = potential.parameters
                return parameters
        raise MissingParametersError(
            f"Could not find parameter in parameter in handler {self.type} "
            f"associated with atoms {atom_indices}",
        )

    def get_force_field_parameters(
        self,
        use_jax: bool = False,
    ) -> Union["ArrayLike", "Array"]:
        """Return a flattened representation of the force field parameters."""
        # TODO: Handle WrappedPotential
        if any(
            isinstance(potential, WrappedPotential)
            for potential in self.potentials.values()
        ):
            raise NotImplementedError

        if use_jax:
            return jax_numpy.array(
                [
                    [v.m for v in p.parameters.values()]
                    for p in self.potentials.values()
                ],
            )
        else:
            return numpy.array(
                [
                    [v.m for v in p.parameters.values()]
                    for p in self.potentials.values()
                ],
            )

    def set_force_field_parameters(self, new_p: "ArrayLike") -> None:
        """Set the force field parameters from a flattened representation."""
        mapping = self.get_mapping()
        if new_p.shape[0] != len(mapping):  # type: ignore
            raise RuntimeError

        for potential_key, potential_index in self.get_mapping().items():
            potential = self.potentials[potential_key]
            if len(new_p[potential_index, :]) != len(potential.parameters):  # type: ignore
                raise RuntimeError

            for parameter_index, parameter_key in enumerate(potential.parameters):
                parameter_units = potential.parameters[parameter_key].units
                modified_parameter = new_p[potential_index, parameter_index]  # type: ignore

                self.potentials[potential_key].parameters[parameter_key] = (
                    modified_parameter * parameter_units
                )

    def get_system_parameters(
        self,
        p=None,
        use_jax: bool = False,
    ) -> Union["ArrayLike", "Array"]:
        """
        Return a flattened representation of system parameters.

        These values are effectively force field parameters as applied to a chemical topology.
        """
        # TODO: Handle WrappedPotential
        if any(
            isinstance(potential, WrappedPotential)
            for potential in self.potentials.values()
        ):
            raise NotImplementedError

        if p is None:
            p = self.get_force_field_parameters(use_jax=use_jax)
        mapping = self.get_mapping()

        q: list = list()
        for potential_key in self.key_map.values():
            index = mapping[potential_key]
            q.append(p[index])

        if use_jax:
            return jax_numpy.array(q)
        else:
            return numpy.array(q)

    def get_mapping(self) -> dict[PotentialKey, int]:
        """Get a mapping between potentials and array indices."""
        mapping: dict = dict()
        index = 0
        for potential_key in self.key_map.values():
            if potential_key not in mapping:
                mapping[potential_key] = index
                index += 1

        return mapping

    def parametrize(
        self,
        p=None,
        use_jax: bool = True,
    ) -> Union["ArrayLike", "Array"]:
        """Return an array of system parameters, given an array of force field parameters."""
        if p is None:
            p = self.get_force_field_parameters(use_jax=use_jax)

        return self.get_system_parameters(p=p, use_jax=use_jax)

    def parametrize_partial(self):
        """Return a function that will call `self.parametrize()` with arguments specified by `self.mapping`."""
        from functools import partial

        return partial(
            self.parametrize,
            mapping=self.get_mapping(),
        )

    @requires_package("jax")
    def get_param_matrix(self) -> Union["Array", "ArrayLike"]:
        """Get a matrix representing the mapping between force field and system parameters."""
        from functools import partial

        import jax

        p = self.get_force_field_parameters(use_jax=True)

        parametrize_partial = partial(
            self.parametrize,
        )

        jac_parametrize = jax.jacfwd(parametrize_partial)
        jac_res = jac_parametrize(p)

        return jac_res.reshape(-1, p.flatten().shape[0])  # type: ignore[union-attr]

    def __getattr__(self, attr: str):
        if attr == "slot_map":
            warnings.warn(
                "The `slot_map` attribute is deprecated. Use `key_map` instead.",
                InterchangeDeprecationWarning,
            )
            return self.key_map
        else:
            return super().__getattribute__(attr)
