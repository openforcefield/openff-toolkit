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

from collections import defaultdict
import copy
from simtk.openmm.app import element as E
from simtk.openmm import CustomGBForce, Discrete2DFunction
import simtk.unit as u
from math import floor


def strip_unit(value, unit):
    """Strip off any units and return value in unit"""
    if not u.is_quantity(value):
        return value
    return value.value_in_unit(unit)

def _createEnergyTerms(force, solventDielectric, soluteDielectric, SA_model, cutoff, kappa, offset):
    """Add the energy terms to the CustomGBForce.

    These are identical for all the GB models.

    """
    params = "; solventDielectric=%.16g; soluteDielectric=%.16g; kappa=%.16g; offset=%.16g" % (solventDielectric, soluteDielectric, kappa, offset)
    if cutoff is not None:
        params += "; cutoff=%.16g" % cutoff
    if kappa > 0:
        force.addEnergyTerm("-0.5*138.935485*(1/soluteDielectric-exp(-kappa*B)/solventDielectric)*charge^2/B"+params,
                CustomGBForce.SingleParticle)
    elif kappa < 0:
        # Do kappa check here to avoid repeating code everywhere
        raise ValueError('kappa/ionic strength must be >= 0')
    else:
        force.addEnergyTerm("-0.5*138.935485*(1/soluteDielectric-1/solventDielectric)*charge^2/B"+params,
                CustomGBForce.SingleParticle)
    if SA_model=='ACE':
        force.addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/B)^6; radius=or+offset"+params, CustomGBForce.SingleParticle)
    elif SA_model is not None:
        raise ValueError('Unknown surface area method: '+SA_model)
    if cutoff is None:
        if kappa > 0:
            force.addEnergyTerm("-138.935485*(1/soluteDielectric-exp(-kappa*f)/solventDielectric)*charge1*charge2/f;"
                                "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
        else:
            force.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*charge1*charge2/f;"
                                "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
    else:
        if kappa > 0:
            force.addEnergyTerm("-138.935485*(1/soluteDielectric-exp(-kappa*f)/solventDielectric)*charge1*charge2*(1/f-"+str(1/cutoff)+");"
                                "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
        else:
            force.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*charge1*charge2*(1/f-"+str(1/cutoff)+");"
                                "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)

class CustomAmberGBForceBase(CustomGBForce):
    """Base class for all of the Amber custom GB forces.

    Should not be instantiated directly, use one of its
    derived classes instead.

    """
    OFFSET = 0.009
    RADIUS_ARG_POSITION = 1
    SCREEN_POSITION = 2

    def __init__(self):
        CustomGBForce.__init__(self)
        self.parameters = []

    # def addParticle(self, params):
    #     """Add a particle to the force
    #
    #     Particles are added in order. The number of particles
    #     added must match the number of particles in the system.
    #
    #     Parameters
    #     ----------
    #     params : list
    #         A list of parameters to add to the force. The meaning
    #         parameters depends on the model.
    #
    #     Returns
    #     -------
    #     list
    #         The list of parameters after stripping off units and
    #         modifying SCREEN.
    #
    #     """
    #     params = copy.deepcopy(params)
    #     params[self.RADIUS_ARG_POSITION] = strip_unit(params[self.RADIUS_ARG_POSITION], u.nanometer) - self.OFFSET
    #     params[self.SCREEN_POSITION] *= params[self.RADIUS_ARG_POSITION]
    #     self.parameters.append(params)
    #     return params
    #
    # def finalize(self):
    #     """Finalize this force so it can be added to a system.
    #
    #     This method must be called before the force is added
    #     to the system.
    #
    #     """
    #     self._addParticles()
    #
    # def _addParticles(self):
    #     for params in self.parameters:
    #         CustomGBForce.addParticle(self, params)


class HCT(CustomAmberGBForceBase):
    """This class is equivalent to Amber ``igb=1``

    The list of parameters to ``addParticle`` is: ``[charge, radius, scale]``.

    Parameters
    ----------
    solventDielectric: float
        Dielectric constant for the solvent
    soluteDielectric: float
        Dielectric constant for the solute
    SA_model: string or None
        Surface area model to use ['ACE', None]
    cutoff: float or Quantity or None
        Cutoff distance to use. If float, value is in nm. If ``None``,
        then no cutoffs are used.
    kappa: float or Quantity
        Debye kappa parameter related to modelling salt in GB. It has
        units of 1 / length with 1 / nanometer assumed if a float
        is given. A value of zero corresponds to zero salt concentration.

    """
    def __init__(self, solventDielectric=78.5, soluteDielectric=1, SA_model=None,
                 cutoff=None, kappa=0.0):
        CustomAmberGBForceBase.__init__(self)

        self.addPerParticleParameter("charge")
        self.addPerParticleParameter("radius") # Offset radius
        self.addPerParticleParameter("scale") # Scaled offset radius
        self.addComputedValue("I", "step(r+scale2-radius1)*0.5*(1/L-1/U+0.25*(r-scale2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
                                   "U=r+scale2;"
                                   "L=max(radius1, D);"
                                   "D=abs(r-scale2)",
                              CustomGBForce.ParticlePairNoExclusions)

        self.addComputedValue("B", "1/(1/radius-I)", CustomGBForce.SingleParticle)
        _createEnergyTerms(self, solventDielectric, soluteDielectric, SA_model, cutoff, kappa, 0.009)

class OBC1(CustomAmberGBForceBase):
    """This class is equivalent to Amber ``igb=2``

    The list of parameters to ``addParticle`` is: ``[charge, radius, scale]``.

    Parameters
    ----------
    solventDielectric: float
        Dielectric constant for the solvent
    soluteDielectric: float
        Dielectric constant for the solute
    SA_model: string or None
        Surface area model to use ['ACE', None]
    cutoff: float or Quantity or None
        Cutoff distance to use. If float, value is in nm. If ``None``,
        then no cutoffs are used.
    kappa: float or Quantity
        Debye kappa parameter related to modelling salt in GB. It has
        units of 1 / length with 1 / nanometer assumed if a float
        is given. A value of zero corresponds to zero salt concentration.

    """
    def __init__(self, solventDielectric=78.5, soluteDielectric=1, SA_model=None,
                 cutoff=None, kappa=0.0):

        CustomAmberGBForceBase.__init__(self)

        self.addPerParticleParameter("charge")
        self.addPerParticleParameter("radius") # Offset radius
        self.addPerParticleParameter("scale") # Scaled offset radius
        self.addComputedValue("I",  "step(r+scale2-radius1)*0.5*(1/L-1/U+0.25*(r-scale2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
                                    "U=r+scale2;"
                                    "L=max(radius1, D);"
                                    "D=abs(r-scale2)", CustomGBForce.ParticlePairNoExclusions)

        self.addComputedValue("B", "1/(1/radius-tanh(0.8*psi+2.909125*psi^3)/radius);"
                                   "psi=I*radius; radius=radius+offset; offset=0.009", CustomGBForce.SingleParticle)
        _createEnergyTerms(self, solventDielectric, soluteDielectric, SA_model, cutoff, kappa, 0.009)

class OBC2(OBC1):
    """This class is equivalent to Amber ``igb=5``

    The list of parameters to ``addParticle`` is: ``[charge, radius, scale]``.

    Parameters
    ----------
    solventDielectric: float
        Dielectric constant for the solvent
    soluteDielectric: float
        Dielectric constant for the solute
    SA_model: string or None
        Surface area model to use ['ACE', None]
    cutoff: float or Quantity or None
        Cutoff distance to use. If float, value is in nm. If ``None``,
        then no cutoffs are used.
    kappa: float or Quantity
        Debye kappa parameter related to modelling salt in GB. It has
        units of 1 / length with 1 / nanometer assumed if a float
        is given. A value of zero corresponds to zero salt concentration.

    """
    def __init__(self, solventDielectric=78.5, soluteDielectric=1, SA_model=None,
                 cutoff=None, kappa=0.0):

        CustomAmberGBForceBase.__init__(self)

        self.addPerParticleParameter("charge")
        self.addPerParticleParameter("radius") # Offset radius
        self.addPerParticleParameter("scale") # Scaled offset radius
        self.addComputedValue("I",  "step(r+scale2-radius1)*0.5*(1/L-1/U+0.25*(r-scale2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
                                    "U=r+scale2;"
                                    "L=max(radius1, D);"
                                    "D=abs(r-scale2)", CustomGBForce.ParticlePairNoExclusions)

        self.addComputedValue("B", "1/(1/radius-tanh(psi-0.8*psi^2+4.85*psi^3)/radius);"
                                     "psi=I*radius; radius=radius+offset; offset=0.009", CustomGBForce.SingleParticle)
        _createEnergyTerms(self, solventDielectric, soluteDielectric, SA_model, cutoff, kappa, 0.009)
