"""
A recreation of the various GB variants implemented via CustomGBForce

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2016 University of Virginia and the Authors.
Authors: Christoph Klein, Michael R. Shirts
Contributors: Jason M. Swails, Peter Eastman, Justin L. MacCallum

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


from __future__ import division, absolute_import

from simtk.openmm import CustomGBForce
from simtk import unit
from math import pi



#=============================================================================================
# CONSTANTS
#=============================================================================================

ONE_4PI_EPS0 = 138.935456 # OpenMM constant for Coulomb interactions (openmm/platforms/reference/include/SimTKOpenMMRealType.h) in OpenMM units
                          # TODO: Replace this with an import from simtk.openmm.constants once these constants are available there

OFFSET = 0.009 # Radius offset (in nm) for all GB models

#=============================================================================================
# SUBROUTINES
#=============================================================================================

def strip_unit(value, unit):
    """Strip off any units and return value in unit"""
    if not unit.is_quantity(value):
        return value
    return value.value_in_unit(unit)

def _get_option_stripped(kwargs, name, default, dtype=None, compatible_units=None):
    """Return the specified option, converted to float in md_unit_system units.

    Parameters
    ----------
    kwargs : dict
       Dictionary from which options are taken.
    name : str
       Name of the option to be retrieved.
    default : simtk.unit.Quantity or float
       Default value
    dtype : type
       If specified, will force to this type
    compatible_units : simtk.unit.Unit
       If not None, will ensure that quantity is compatible with these units.
    """
    if name in kwargs:
        x = kwargs[name]
        # Force to specified type
        if dtype is unit.Quantity:
            x = eval(x, unit.__dict__)
        elif dtype is not None:
            x = dtype(x)
    else:
        x = default

    if not unit.is_quantity(x):
        return x
    else:
        if (compatible_units is not None):
            # Check unit compatibility, raising exception if not compatible
            x.in_units_of(compatible_units)

        return default.value_in_unit_system(unit.md_unit_system)

#=============================================================================================
# GBSA MODELS
#=============================================================================================

class CustomAmberGBForceBase(CustomGBForce):
    """Base class for all of the Amber custom GBSA forces.

    Should not be instantiated directly, use one of its derived classes instead.
    """

    def __init__(self, **kwargs):
        CustomGBForce.__init__(self)

        self.addPerParticleParameter("charge")
        self.addPerParticleParameter("radius")
        self.addPerParticleParameter("scale")

        self.offset_terms_single = ("sr=scale*or;"
                                    "or=(radius-OFFSET);"
                                    "OFFSET=%.16f;" % OFFSET)

        self.offset_terms_pair = ("sr1=scale1*or1;"
                                  "or1=(radius1-OFFSET);"
                                  "sr2=scale2*or2;"
                                  "or2=(radius2-OFFSET);"
                                  "OFFSET=%.16f;" % OFFSET)

    def _createGBEnergyTerms(self, **kwargs):
        """Add energy terms for the GB model to the CustomGBForce.

        Parametes
        ---------
        cutoff : simtk.unit.Quantity with units compatible with distance
           If a cutoff is not None, the cutoff for GB interactions.
        kwargs : dict
           Optional arguments required for GB models (e.g. 'solventDielectric', 'soluteDielectric', 'kappa')

        """

        # Construct GB energy function
        energy_expression = "-0.5*ONE_4PI_EPS0*(1/soluteDielectric-kappa_coeff/solventDielectric)*charge^2/B"

        # Add dielectric constants
        solvent_dielectric = _get_option_stripped(kwargs, 'solvent_dielectric', 78.5, dtype=float)
        solute_dielectric = _get_option_stripped(kwargs, 'solute_dielectric', 1.0, dtype=float)
        energy_expression += "; solventDielectric=%.16g; soluteDielectric=%.16g" % (solvent_dielectric, solute_dielectric)

        # Salt screening term
        kappa = _get_option_stripped(kwargs, 'kappa', None, dtype=unit.Quantity, compatible_units=1.0/unit.nanometers)
        if kappa is not None:
            if (kappa < 0):
                raise ValueError('kappa/ionic strength must be >= 0')
            energy_expression += "; kappa_coeff = exp(-kappa*B); kappa=%.16f" % (kappa)
        else:
            energy_expression += "; kappa_coeff = 1"

        # Add constants
        energy_expression += "; ONE_4PI_EPS0=%.16g" % (ONE_4PI_EPS0)

        # Add force term
        self.addEnergyTerm(energy_expression, CustomGBForce.SingleParticle)

    def _createSAEnergyTerms(self, sa_model=None, **kwargs):
        """Add the energy terms for the SA model to the CustomGBForce.

        """

        if sa_model=='ACE':
            surface_area_penalty = _get_option_stripped(kwargs, 'surface_area_penalty', 5.4 * unit.calories / unit.mole / unit.angstrom**2,
                dtype=unit.Quantity, compatible_units=unit.calories / unit.mole / unit.angstrom**2)
            solvent_radius = _get_option_stripped(kwargs, 'solvent_radius', 0.14 * unit.nanometers,
                dtype=unit.Quantity, compatible_units=unit.nanometers)
            energy_expression  = 'surface_area_penalty*4*pi*(radius+solvent_radius)^2*(radius/B)^6'
            energy_expression += '; pi=%.16g;' % pi
            energy_expression += '; surface_area_penalty=%.16g;' % surface_area_penalty
            energy_expression += '; solvent_radius=%.16f' % solvent_radius
            self.addEnergyTerm(energy_expression, CustomGBForce.SingleParticle)

        elif sa_model is not None:
            raise ValueError("Unknown surface area method '%s'. Must be one of ['ACE', None]" % (sa_model))

class HCT(CustomAmberGBForceBase):
    """This class is equivalent to Amber ``igb=1``

    The list of parameters to ``addParticle`` is: ``[charge, radius, scale]``.
    """
    def __init__(self, **kwargs):
        CustomAmberGBForceBase.__init__(self, **kwargs)

        I_expression = ("step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
                        "U=r+sr2;"
                        "L=max(or1, D);"
                        "D=abs(r-sr2);") + self.offset_terms_pair
        self.addComputedValue("I", I_expression, CustomGBForce.ParticlePairNoExclusions)

        B_expression = "1/(1/or-I);" + self.offset_terms_single
        self.addComputedValue("B", B_expression, CustomGBForce.SingleParticle)

        self._createGBEnergyTerms(**kwargs)
        self._createSAEnergyTerms(**kwargs)

class OBC1(CustomAmberGBForceBase):
    """This class is equivalent to Amber ``igb=2``

    The list of parameters to ``addParticle`` is: ``[charge, radius, scale]``.
    """
    def __init__(self, **kwargs):
        CustomAmberGBForceBase.__init__(self, **kwargs)

        I_expression = ("step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
                        "U=r+sr2;"
                        "L=max(or1, D);"
                        "D=abs(r-sr2);") + self.offset_terms_pair
        self.addComputedValue("I",  I_expression, CustomGBForce.ParticlePairNoExclusions)

        B_expression = ("1/(1/or-tanh(0.8*psi+2.909125*psi^3)/radius);"
                        "psi=I*or;") + self.offset_terms_single
        self.addComputedValue("B", B_expression, CustomGBForce.SingleParticle)

        self._createGBEnergyTerms(**kwargs)
        self._createSAEnergyTerms(**kwargs)

class OBC2(CustomAmberGBForceBase):
    """This class is equivalent to Amber ``igb=5``

    The list of parameters to ``addParticle`` is: ``[charge, radius, scale]``.
    """
    def __init__(self, **kwargs):
        CustomAmberGBForceBase.__init__(self, **kwargs)

        I_expression = ("step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
                        "U=r+sr2;"
                        "L=max(or1, D);"
                        "D=abs(r-sr2);") + self.offset_terms_pair
        self.addComputedValue("I",  I_expression, CustomGBForce.ParticlePairNoExclusions)

        B_expression = ("1/(1/or-tanh(psi-0.8*psi^2+4.85*psi^3)/radius);"
                        "psi=I*or;") + self.offset_terms_single
        self.addComputedValue("B", B_expression, CustomGBForce.SingleParticle)

        self._createGBEnergyTerms(**kwargs)
        self._createSAEnergyTerms(**kwargs)
