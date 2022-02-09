"""Re-exports of concrete ParameterTypes

openff.toolkit.typing.engines.smirnoff.parameters defines a number of parameter
types within class definitions of the corresponding ParameterHandler. This module
re-exports them for discoverability and ease of use."""

import openff.toolkit.typing.engines.smirnoff.parameters as parameters

ParameterType = parameters.ParameterType

ConstraintType = parameters.ConstraintHandler.ConstraintType
BondType = parameters.BondHandler.BondType
AngleType = parameters.AngleHandler.AngleType
ProperTorsionType = parameters.ProperTorsionHandler.ProperTorsionType
ImproperTorsionType = parameters.ImproperTorsionHandler.ImproperTorsionType
vdWType = parameters.vdWHandler.vdWType
LibraryChargeType = parameters.LibraryChargeHandler.LibraryChargeType
GBSAType = parameters.GBSAHandler.GBSAType
ChargeIncrementType = parameters.ChargeIncrementModelHandler.ChargeIncrementType
VirtualSiteBondChargeType = parameters.VirtualSiteHandler.VirtualSiteBondChargeType
VirtualSiteMonovalentLonePairType = (
    parameters.VirtualSiteHandler.VirtualSiteMonovalentLonePairType
)
VirtualSiteDivalentLonePairType = (
    parameters.VirtualSiteHandler.VirtualSiteDivalentLonePairType
)
VirtualSiteTrivalentLonePairType = (
    parameters.VirtualSiteHandler.VirtualSiteTrivalentLonePairType
)

__all__ = [
    "ParameterType",
    "ConstraintType",
    "BondType",
    "AngleType",
    "ProperTorsionType",
    "ImproperTorsionType",
    "vdWType",
    "LibraryChargeType",
    "GBSAType",
    "ChargeIncrementType",
    "VirtualSiteBondChargeType",
    "VirtualSiteMonovalentLonePairType",
    "VirtualSiteDivalentLonePairType",
    "VirtualSiteTrivalentLonePairType",
]
