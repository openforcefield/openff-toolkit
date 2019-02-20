#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Test classes and function in module openforcefield.typing.engines.smirnoff.parameters.

"""


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================


from openforcefield.typing.engines.smirnoff.parameters import ParameterList, ParameterType, BondHandler, SMIRNOFFSpecError

import pytest

#=============================================================================================
# PARAMETER LIST
#=============================================================================================


class TestParameterList:
    """Test capabilities of ParameterList for accessing and manipulating SMIRNOFF parameter definitions.
    """

    def test_create(self):
        """Test creation of a parameter list.
        """
        p1 = ParameterType(smirks='[*:1]')
        p2 = ParameterType(smirks='[#1:1]')
        parameters = ParameterList([p1, p2])

    def test_getitem(self):
        """Test ParameterList __getitem__ overloading.
        """
        p1 = ParameterType(smirks='[*:1]')
        p2 = ParameterType(smirks='[#1:1]')
        parameters = ParameterList([p1, p2])
        assert parameters[0] == p1
        assert parameters[1] == p2
        assert parameters[p1.smirks] == p1
        assert parameters[p2.smirks] == p2

        # Note that this call access __getitem__, not __setitem__.
        parameters['[*:1]'].smirks = '[*X4:1]'
        assert parameters[0].smirks == '[*X4:1]'
        assert p1.smirks == '[*X4:1]'

    def test_contains(self):
        """Test ParameterList __contains__ overloading.
        """
        p1 = ParameterType(smirks='[*:1]')
        p2 = ParameterType(smirks='[#1:1]')
        p3 = ParameterType(smirks='[#7:1]')
        parameters = ParameterList([p1, p2])
        assert p1 in parameters
        assert p2 in parameters
        assert p3 not in parameters
        assert p1.smirks in parameters
        assert p2.smirks in parameters
        assert p3.smirks not in parameters

    def test_del(self):
        """
        Test ParameterList __del__ overloading.
        """
        p1 = ParameterType(smirks='[*:1]')
        p2 = ParameterType(smirks='[#1:1]')
        p3 = ParameterType(smirks='[#7:1]')
        parameters = ParameterList([p1, p2, p3])

        # Test that original list deletion behavior is preserved.
        del parameters[2]
        assert len(parameters) == 2
        assert p1 in parameters
        assert p2 in parameters
        assert p3 not in parameters

        # Test that we can delete elements by their smirks.
        del parameters['[#1:1]']
        assert len(parameters) == 1
        assert p1 in parameters
        assert p2 not in parameters

class TestParameterType:

    def test_base_parametertype_to_dict(self):
        """
        Test ParameterType to_dict.
        """
        p1 = ParameterType(smirks='[*:1]')
        param_dict, attached_units = p1.to_dict()
        assert param_dict['smirks'] == '[*:1]'
        assert len(param_dict.keys()) == 1

    def test_bondtype_to_dict(self):
        """
        Test BondType to_dict.
        """
        from simtk import unit

        p1 = BondHandler.BondType(smirks='[*:1]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        param_dict, attached_units = p1.to_dict()
        assert param_dict == {'smirks': '[*:1]',
                              'length': 1.02,
                              'k': 5,}
        assert attached_units == {'length_unit': [('angstrom', 1)],
                                  'k_unit': [('angstrom', -2), ('mole', -1), ('kilocalorie', 1)]
                                  }


    def test_bondtype_to_dict_custom_output_units(self):
        """
        Test BondType to_dict with custom output units.
        """
        from simtk import unit
        p1 = BondHandler.BondType(smirks='[*:1]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        param_dict, attached_units = p1.to_dict(output_units={'length': unit.nanometer})
        assert attached_units['length_unit'] == [('nanometer', 1)]
        assert abs(param_dict['length'] - 0.102) < 1e-10


    def test_bondtype_to_dict_invalid_output_units(self):
        """
        Test ParameterType to_dict with invalid output units.
        """
        from simtk import unit
        p1 = BondHandler.BondType(smirks='[*:1]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        with pytest.raises(ValueError,
                           match='Requested output unit calorie is not compatible with quantity unit angstrom .'
                           ) as context:
            param_dict, attached_units = p1.to_dict(output_units={'length': unit.calorie})

    def test_read_write_cosmetic_parameter_attribute(self):
        """
        Test ParameterTypes' ability to store and write out cosmetic attributes passed to __init__()
        """
        from simtk import unit

        p1 = BondHandler.BondType(smirks='[*:1]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2,
                                  pilot='alice',
                                  permit_cosmetic_attributes=True
                                  )
        param_dict, attached_units = p1.to_dict(return_cosmetic_attributes=True)
        assert ('pilot', 'alice') in param_dict.items()

    def test_read_but_dont_write_cosmetic_parameter_attribute(self):
        """
        Test ParameterTypes' ability to ignore cosmetic attributes passed to __init__() if instructed
        """
        from simtk import unit

        p1 = BondHandler.BondType(smirks='[*:1]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2,
                                  pilot='alice',
                                  permit_cosmetic_attributes=True
                                  )
        param_dict, attached_units = p1.to_dict(return_cosmetic_attributes=False)
        assert ('pilot', 'alice') not in param_dict

    def test_error_cosmetic_parameter_attribute(self):
        """
        Test that ParameterTypes raise an error on receiving unexpected attributes passed to __init__()
        """
        from simtk import unit

        with pytest.raises(SMIRNOFFSpecError, match="Incompatible kwarg {'pilot': 'alice'}") as context:
            p1 = BondHandler.BondType(smirks='[*:1]',
                                      length=1.02*unit.angstrom,
                                      k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2,
                                      pilot='alice',
                                      permit_cosmetic_attributes=False
                                      )
