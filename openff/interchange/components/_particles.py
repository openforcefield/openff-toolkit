"""
Temporary module for second-class virtual site objects.
"""
import abc
import math
from typing import Literal

from openff.models.models import DefaultModel
from openff.models.types import FloatQuantity
from openff.units import unit


class _VirtualSite(DefaultModel, abc.ABC):
    type: str
    distance: FloatQuantity["nanometer"]
    orientations: tuple[int, ...]

    @abc.abstractproperty
    def local_frame_weights(self) -> tuple[list[float], ...]:
        raise NotImplementedError

    def local_frame_positions(self) -> unit.Quantity:
        raise NotImplementedError


class _BondChargeVirtualSite(_VirtualSite):
    type: Literal["BondCharge"]
    distance: FloatQuantity["nanometer"]
    orientations: tuple[int, ...]

    @property
    def local_frame_weights(self) -> tuple[list[float], ...]:
        originwt = [1.0, 0.0]  # first atom is origin
        xdir = [-1.0, 1.0]
        ydir = [-1.0, 1.0]

        return originwt, xdir, ydir

    @property
    def local_frame_positions(self) -> unit.Quantity:
        distance_unit = self.distance.units
        return unit.Quantity(
            [-self.distance.m_as(distance_unit), 0.0, 0.0],
            distance_unit,
        )


class _MonovalentLonePairVirtualSite(_VirtualSite):
    type: Literal["MonovalentLonePair"]
    distance: FloatQuantity["nanometer"]
    out_of_plane_angle: FloatQuantity["degree"]
    in_plane_angle: FloatQuantity["degree"]
    orientations: tuple[int, ...]

    @property
    def local_frame_weights(self) -> tuple[list[float], ...]:
        originwt = [1.0, 0.0, 0.0]
        xdir = [-1.0, 1.0, 0.0]
        ydir = [-1.0, 0.0, 1.0]

        return originwt, xdir, ydir

    @property
    def local_frame_positions(self) -> unit.Quantity:
        theta = self.in_plane_angle.m_as(unit.radian)
        phi = self.out_of_plane_angle.m_as(unit.radian)

        distance_unit = self.distance.units

        return unit.Quantity(
            [
                self.distance.m_as(distance_unit) * math.cos(theta) * math.cos(phi),
                self.distance.m_as(distance_unit) * math.sin(theta) * math.cos(phi),
                self.distance.m_as(distance_unit) * math.sin(phi),
            ],
            distance_unit,
        )


class _DivalentLonePairVirtualSite(_VirtualSite):
    type: Literal["DivalentLonePair"]
    distance: FloatQuantity["nanometer"]
    out_of_plane_angle: FloatQuantity["degree"]
    orientations: tuple[int, ...]

    @property
    def local_frame_weights(self) -> tuple[list[float], ...]:
        originwt = [1.0, 0.0, 0.0]
        xdir = [-1.0, 0.5, 0.5]
        ydir = [-1.0, 1.0, 0.0]

        return originwt, xdir, ydir

    @property
    def local_frame_positions(self) -> unit.Quantity:
        theta = self.out_of_plane_angle.m_as(unit.radian)

        distance_unit = self.distance.units

        return unit.Quantity(
            [
                -self.distance.m_as(distance_unit) * math.cos(theta),
                0.0,
                self.distance.m_as(distance_unit) * math.sin(theta),
            ],
            distance_unit,
        )


class _TrivalentLonePairVirtualSite(_VirtualSite):
    type: Literal["TrivalentLonePair"]
    distance: FloatQuantity["nanometer"]
    orientations: tuple[int, ...]

    @property
    def local_frame_weights(self) -> tuple[list[float], ...]:
        originwt = [1.0, 0.0, 0.0, 0.0]
        xdir = [-1.0, 1 / 3, 1 / 3, 1 / 3]
        ydir = [-1.0, 1.0, 0.0, 0.0]  # Not used

        return originwt, xdir, ydir

    @property
    def local_frame_positions(self) -> unit.Quantity:
        distance_unit = self.distance.units
        return unit.Quantity(
            [-self.distance.m_as(distance_unit), 0.0, 0.0],
            distance_unit,
        )
