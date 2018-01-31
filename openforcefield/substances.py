#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Substances API.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>
* Levi N. Naden <levi.naden@choderalab.org>

TODO
----
* Add methods that construct real System and Topology objects for a specified system size, following the Mobley SolvationToolkit:
  https://github.com/MobleyLab/SolvationToolkit

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import copy
import numpy as np

from simtk import unit
from simtk.openmm import app
from openforcefield.packmol import pack_box
from openeye import oechem, oeiupac

from typing import Union, Tuple


# =============================================================================================
# Component
# =============================================================================================

class Component(object):

    _cached_molecules = dict()  # store cached molecules by IUPAC name

    def _create_molecule(self, iupac_name: str) -> oechem.OEMol:
        """Create molecule from IUPAC name.

        Best practices for creating initial coordinates can be applied here.

        Parameters
        ----------
        iupac_name : str
            IUPAC name

        Returns
        -------
        molecule : OEMol
            OEMol with 3D coordinates, but no charges

        """
        # Check cache
        if iupac_name in self._cached_molecules:
            return copy.deepcopy(self._cached_molecules[iupac_name])

        # Create molecule from IUPAC name.
        molecule = oechem.OEMol()
        if not oeiupac.OEParseIUPACName(molecule, iupac_name):
            # raise ValueError("The supplied IUPAC name '{}' could not be parsed.".format(iupac_name))
            # Pass for now to get something
            molecule.SetTitle(iupac_name)
            return molecule

        # Set molecule name
        molecule.SetTitle(iupac_name)

        # Normalize molecule
        oechem.OEAssignAromaticFlags(molecule, oechem.OEAroModelOpenEye)
        oechem.OEAddExplicitHydrogens(molecule)
        oechem.OETriposAtomNames(molecule)

        return molecule

    @staticmethod
    def _get_iupac_from_molecule(molecule: oechem.OEMol) -> str:
        """

        Parameters
        ----------
        molecule : oechem.OEMolecule
            Input OpenEye molecule to derive IUPAC name from

        Returns
        -------
        name : string
            IUPAC Name as a string

        """
        return oeiupac.OECreateIUPACName(molecule)

    def __init__(self, name):
        """Create a chemical component.

        Parameters
        ----------
        name : str
            IUPAC name of component
        """
        # TODO: Check that IUPAC name is not stereochemically ambiguous.

        self.molecule = self._create_molecule(name)
        self.name = name
        self.iupac_name = self._get_iupac_from_molecule(self.molecule)


# =============================================================================================
# SUBSTANCE
# =============================================================================================

class Substance(object):
    """
    A substance, can be a pure chemical, or could be a Mixture.

    This class is not specific enough to be a chemical species all on its own
    """
    pass


# =============================================================================================
# MIXTURE
# =============================================================================================

class Mixture(Substance):
    """A liquid or gas mixture.

    Properties
    ----------
    components : dict
        components[iupac_name] is the mole fraction of the specified component

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
        def __init__(self, name, mole_fraction=0.0, impurity=False):
            self.mole_fraction = mole_fraction
            self.impurity = impurity
            super().__init__(name)

    def __init__(self):
        """Create a Mixture.
        """
        self.components = list()

    @property
    def total_mole_fraction(self) -> float:
        """Compute the total mole fraction.
        """
        return sum([component.mole_fraction for component in self.components])

    @property
    def n_components(self) -> int:
        return len(self.components)

    @property
    def n_impurities(self) -> int:
        return sum([1 for component in self.components if component.impurity is True])

    def add_component(self, name: str, mole_fraction: Union[None, float]=None, impurity: bool=False):
        """Add a component to the mixture.

        Parameters
        ----------
        name : str
            Name of the component, either common ThermoML or IUPAC name
        mole_fraction : float or None, optional, default=None
            If specified, the mole fraction of this component as a float on the domain [0,1]
            If not specified, this will be the last or only component of the mixture.
        impurity : bool, optional, default=False
            If True, the component represents an impurity (single molecule).
            This is distinct from 0 mole fraction
        """

        mole_fraction, impurity = self._validate_mol_fraction(mole_fraction, impurity)

        component = self.MixtureComponent(name, mole_fraction=mole_fraction, impurity=impurity)
        self.components.append(component)

    def get_component(self, name: str) -> MixtureComponent:
        """Retrieve component by name.

        Parameters
        ----------
        name : str
            The name of the component to retrieve
            Accepts IUPAC or common name

        """
        for component in self.components:
            found = False
            if component.iupac_name == name:
                found = True
            elif component.name == name:
                found = True
            if found:
                return component
        raise Exception("No component with name '{0:s}' found.".format(name))

    def build(self, n_molecules: int=1000,
              mass_density: Union[None, float, unit.Quantity]=None
              ) -> Tuple[app.Topology, unit.Quantity, unit.Quantity]:
        """Build an instance of this mixture.

        Parameters
        ----------
        n_molecules : int, optional, default=True
            The number of molecules in the system to be created.
        mass_density : float, simtk.unit.Quantity, or None; optional, default=None
            If provided, will aid in the selecting an initial box size.

        Returns
        -------
        topology : simtk.openmm.Topology
            Topology object describing the system.
        molecules : list of oechem.OEMol
            The molecules describing the system (not guaranteed to have charges or 3D coordinates).
            These are copies of molecules.
        positions : simtk.unit.Quantity wrapping [n_atoms,3] numpy array with units compatible with angstroms
            Positions of all atoms in the system.

        Notes
        -----
        The number of molecules of each component need not be deterministic.
        Repeated calls may generate different numbers of molecules, orders, and positions.
        Impurities will have exactly one molecule per impurity.

        """

        # Create deep copies of molecules.
        molecules = [copy.deepcopy(component.molecule) for component in self.components]

        # Determine how many molecules of each type will be present in the system.
        mole_fractions = np.array([component.mole_fraction for component in self.components])
        n_copies = np.random.multinomial(n_molecules - self.n_impurities, pvals=mole_fractions)

        # Each impurity must have exactly one molecule
        for (index, component) in enumerate(self.components):
            if component.impurity:
                n_copies[index] = 1

        # Create packed box
        topology, positions = pack_box(molecules, n_copies, mass_density=mass_density)

        return topology, molecules, positions

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
