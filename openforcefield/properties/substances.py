# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Substances API.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>
* Levi N. Naden <levi.naden@choderalab.org>
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from pydantic import BaseModel, validator
from pydantic.validators import dict_validator
from typing import List


# =============================================================================================
# Component
# =============================================================================================

class Component(BaseModel):
    """Represents a chemical component.

    .. todo:: Refactor this out?

     Attributes
     ----------
     smiles : str
         SMILES descriptor of the component
     """

    smiles: str = None

    @classmethod
    def __get_validators__(cls):
        # yield dict_validator
        yield cls.validate

    @classmethod
    def validate(cls, value):
        # A dirty hack to ensure proper inheritance..

        if isinstance(value, cls):
            return value
        else:

            if 'mole_fraction' in value and 'impurity' in value:
                return Mixture.MixtureComponent(**value)

            cls(**dict_validator(value))


# =============================================================================================
# SUBSTANCE
# =============================================================================================

class Substance(BaseModel):
    """
    A substance, can be a pure chemical, or could be a Mixture.

    This class is not specific enough to be a chemical species all on its own
    """

    def __hash__(self):
        raise NotImplementedError('A Substance is a purely abstract base class.')

    def __eq__(self, other):
        raise NotImplementedError('A Substance is a purely abstract base class.')

    def __ne__(self, other):
        raise NotImplementedError('A Substance is a purely abstract base class.')

    @classmethod
    def __get_validators__(cls):
        # yield dict_validator
        yield cls.validate

    @classmethod
    def validate(cls, value):
        # A dirty hack to ensure proper inheritance..

        if isinstance(value, cls):
            return value
        else:

            if 'components' in value:
                return Mixture(**value)

            cls(**dict_validator(value))


# =============================================================================================
# MIXTURE
# =============================================================================================

# TODO: The name is perhaps misleading as a mixture can be pure...
# SystemComposition or just Composition perhaps?
class Mixture(Substance):
    """Defines the components and their amounts in a mixture.

    Examples
    --------

    A neat liquid has only one component:

    >>> liquid = Mixture()
    >>> liquid.add_component('water')

    A binary mixture has two components:

    >>> binary_mixture = Mixture()
    >>> binary_mixture.add_component('water', mole_fraction=0.2)
    >>> binary_mixture.add_component('methanol') # assumed to be rest of mixture if no mole_fraction specified

    A ternary mixture has three components:

    >>> ternary_mixture = Mixture()
    >>> binary_mixture.add_component('ethanol', mole_fraction=0.2)
    >>> binary_mixture.add_component('methanol', mole_fraction=0.2)
    >>> ternary_mixture.add_component('water')

    The infinite dilution of one solute within a solvent or mixture is also specified as a `Mixture`, where the solute
    has is treated as an impurity, and so only 1 atom is added:

    >>> infinite_dilution = Mixture()
    >>> infinite_dilution.add_component('phenol', impurity=True) # infinite dilution
    >>> infinite_dilution.add_component('water')

    """

    class MixtureComponent(Component):
        """Subclass of Component which has mole_fractions and impurity"""

        mole_fraction: float = 0.0
        impurity: bool = False

        def __str__(self):

            hash_value = self.smiles

            if self.mole_fraction is not None:
                hash_value += "{%s}" % str(self.mole_fraction)
            elif self.impurity is not None:
                hash_value += "{%s}" % str(self.impurity)

            return hash_value

        def __hash__(self):
            return hash(str(self))

        def __eq__(self, other):

            return (self.smiles == other.self and
                    self.mole_fraction == other.mole_fraction and
                    self.impurity == other.impurity)

        def __ne__(self, other):
            return not (self == other)

    components: List[MixtureComponent] = list()

    @property
    def total_mole_fraction(self):
        """Compute the total mole fraction.
        """
        return sum([component.mole_fraction for component in self.components])

    @property
    def number_of_components(self):
        return len(self.components)

    @property
    def number_of_impurities(self):
        return sum([1 for component in self.components if component.impurity is True])

    def add_component(self, smiles, mole_fraction, impurity=False):
        """Add a component to the mixture.

        Parameters
        ----------
        smiles : str
            SMILES pattern of the component
        mole_fraction : float or None, optional, default=None
            If specified, the mole fraction of this component as a float on the domain [0,1]
            If not specified, this will be the last or only component of the mixture.
        impurity : bool, optional, default=False
            If True, the component represents an impurity (single molecule).
            This is distinct from 0 mole fraction
        """

        mole_fraction, impurity = self._validate_mol_fraction(mole_fraction, impurity)

        component = self.MixtureComponent(smiles=smiles, mole_fraction=mole_fraction, impurity=impurity)
        self.components.append(component)

    def get_component(self, smiles: str):
        """Retrieve component by name.

        Parameters
        ----------
        smiles : str
            The smiles of the component to retrieve

        """
        for component in self._components:

            if component.smiles != smiles:
                continue

            return component

        raise Exception("No component with smiles '{0:s}' found.".format(smiles))

    def _validate_mol_fraction(self, mole_fraction, impurity):
        """
        Validates the mole_fraction and impurity, setting the defaults if need be.
        See :func:``add_component`` for parameters.
        """
        if not impurity and mole_fraction is None:
            raise ValueError("Either mole_fraction or impurity must be specified!")
        elif impurity and mole_fraction != 0:
            raise ValueError('Mole fraction must be 0.0 or None for impurities. '
                             'Specified mole fraction of {0:f}'.format(mole_fraction))
        elif mole_fraction is not None and not 0.0 <= mole_fraction <= 1.0:
            raise ValueError('Mole fraction must be positive; specified {0:f}.'.format(mole_fraction))
        if impurity:
            mole_fraction = 0.0
        if mole_fraction is None:
            mole_fraction = 1.0 - self.total_mole_fraction
        if (self.total_mole_fraction + mole_fraction) > 1.0:
            raise ValueError("Total mole fraction would exceed "
                             "unity ({0:f}); specified {1:f}".format(self.total_mole_fraction, mole_fraction))

        return mole_fraction, impurity

    def __str__(self):

        hash_tags = [str(component) for component in self.components]
        hash_tags.sort()

        return "|".join(hash_tags)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __ne__(self, other):
        return not (self == other)
