__all__ = (
    "DEFAULT_AROMATICITY_MODEL",
    "ALLOWED_AROMATICITY_MODELS",
    "DEFAULT_FRACTIONAL_BOND_ORDER_MODEL",
    "ALLOWED_FRACTIONAL_BOND_ORDER_MODELS",
    "DEFAULT_CHARGE_MODEL",
    "ALLOWED_CHARGE_MODELS",
)


# =============================================================================================
# SUPPORTED MODELS
#
# TODO: We may no longer need these since we now require SMIRNOFF to specify these models explicitly.
# =============================================================================================

DEFAULT_AROMATICITY_MODEL = "OEAroModel_MDL"  # TODO: Is there a more specific name and reference for the aromaticity model?
ALLOWED_AROMATICITY_MODELS = ["OEAroModel_MDL"]

DEFAULT_FRACTIONAL_BOND_ORDER_MODEL = "Wiberg"  # TODO: Is there a more specific name and reference for the fractional bond order models?
ALLOWED_FRACTIONAL_BOND_ORDER_MODELS = ["Wiberg"]

DEFAULT_CHARGE_MODEL = "AM1-BCC"  # TODO: Should this be `AM1-BCC`, or should we encode BCCs explicitly via AM1-CM2 preprocessing?
ALLOWED_CHARGE_MODELS = ["AM1-BCC"]  # TODO: Which models do we want to support?
