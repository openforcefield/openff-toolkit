#!/usr/bin/env python

#======================================================================
# MODULE DOCSTRING
#======================================================================

"""
Test classes and function in module openforcefield.typing.engines.smirnoff.parameters.

"""


#======================================================================
# GLOBAL IMPORTS
#======================================================================

from simtk import unit
import pytest

from openforcefield.typing.engines.smirnoff.parameters import (
    _ParameterAttribute, ParameterList, ParameterType, BondHandler,
    ParameterHandler, ProperTorsionHandler, ImproperTorsionHandler,
    ToolkitAM1BCCHandler, SMIRNOFFSpecError
)
from openforcefield.typing.engines.smirnoff import SMIRNOFFVersionError
from openforcefield.utils import detach_units, IncompatibleUnitError


#======================================================================
# Test ParameterAttribute descriptor
#======================================================================

class TestParameterAttribute:
    """Test cases for the descriptor _ParameterAttribute."""

    def test_default_value(self):
        """Default values are assigned correctly on initialization."""
        class MyParameter:
            attr_optional = _ParameterAttribute(default=2)
        my_par = MyParameter()
        assert my_par.attr_optional == 2

    def test_none_default_value(self):
        """None is a valid default value for ParameterAttribute."""
        class MyParameter:
            attr_optional = _ParameterAttribute(default=None)
        my_par = MyParameter()
        assert my_par.attr_optional is None

    def test_mandatory_value(self):
        """AttributeError is raised if a mandatory attribute is accessed before being initialized."""
        class MyParameter:
            attr_mandatory = _ParameterAttribute()
        my_par = MyParameter()
        with pytest.raises(AttributeError):
            my_par.attr_mandatory

    def test_unit_validation(self):
        """ParameterAttributes attached to a unit are validated correctly."""
        class MyParameter:
            attr_unit = _ParameterAttribute(unit=unit.kilocalories_per_mole/unit.angstrom**2)
        my_par = MyParameter()

        # TypeError is raised when setting a unit-less value.
        with pytest.raises(TypeError, match='should have units of'):
            my_par.attr_unit = 3.0
        # TypeError is raised when setting a value with incorrect units.
        with pytest.raises(TypeError, match='should have units of'):
            my_par.attr_unit = 3.0 * unit.kilocalories_per_mole

        # Otherwise the attribute is assigned correctly.
        value = 3.0 * unit.kilocalories_per_mole/unit.angstrom**2
        my_par.attr_unit = value
        assert my_par.attr_unit == value

    def test_custom_validator(self):
        """Custom validators of ParameterAttributes are executed correctly."""
        class MyParameter:
            attr_all_to_float = _ParameterAttribute(validator=float)
            attr_int_to_float = _ParameterAttribute()
            @attr_int_to_float.validator
            def attr_int_to_float(self, value):
                """Convert only integers to float"""
                if isinstance(value, int):
                    return float(value)
                elif not isinstance(value, float):
                    raise TypeError()
                return value

        my_par = MyParameter()

        # Both strings and integers are converted to floats when casted with float().
        my_par.attr_all_to_float = '1.0'
        assert isinstance(my_par.attr_all_to_float, float) and my_par.attr_all_to_float == 1.0
        my_par.attr_all_to_float = 2
        assert isinstance(my_par.attr_all_to_float, float) and my_par.attr_all_to_float == 2.0

        # Only integers are converted with the custom validator
        with pytest.raises(TypeError):
            my_par.attr_int_to_float = '1.0'
        my_par.attr_int_to_float = 2
        assert isinstance(my_par.attr_int_to_float, float) and my_par.attr_int_to_float == 2.0

    def test_incosistent_default(self):
        """An exception is raised when a default that doesn't pass validation is used."""
        with pytest.raises(TypeError, match='default value None does not pass validation'):
            class MyParameter:
                attr_inconsistent = _ParameterAttribute(default=None, validator=float)


#======================================================================
# Test ParameterHandler
#======================================================================

class TestParameterHandler:

    def test_different_units_to_dict(self):
        """Test ParameterHandler.to_dict() function when some parameters are in
        different units (proper behavior is to convert all quantities to the last-
        read unit)
        """
        from simtk import unit
        bh = BondHandler(skip_version_check=True)
        bh.add_parameter({'smirks': '[*:1]-[*:2]',
                          'length': 1*unit.angstrom,
                          'k': 10*unit.kilocalorie_per_mole/unit.angstrom**2})
        bh.add_parameter({'smirks': '[*:1]=[*:2]',
                          'length': 0.2*unit.nanometer,
                          'k': 0.4*unit.kilojoule_per_mole/unit.nanometer**2})
        bh_dict = bh.to_dict()
        assert bh_dict['Bond'][0]['length'] == unit.Quantity(value=1, unit=unit.angstrom)
        assert bh_dict['Bond'][1]['length'] == unit.Quantity(value=2, unit=unit.angstrom)

    def test_to_dict_maintain_units(self):
        """Test ParameterHandler.to_dict() function when parameters were provided in different units
        """
        from simtk import unit
        bh = BondHandler(skip_version_check=True)
        bh.add_parameter({'smirks': '[*:1]-[*:2]',
                          'length': 1*unit.angstrom,
                          'k': 10*unit.kilocalorie_per_mole/unit.angstrom**2})
        bh.add_parameter({'smirks': '[*:1]=[*:2]',
                          'length': 0.2*unit.nanometer,
                          'k': 0.4*unit.kilojoule_per_mole/unit.nanometer**2})
        bh_dict = bh.to_dict()
        assert bh_dict['Bond'][0]['length'] == unit.Quantity(1., unit.angstrom)
        assert bh_dict['Bond'][0]['length'].unit == unit.angstrom
        assert bh_dict['Bond'][1]['length'] == unit.Quantity(0.2, unit.nanometer)
        assert bh_dict['Bond'][1]['length'].unit == unit.nanometer


    def test_missing_section_version(self):
        """Test that exceptions are raised if invalid or improper section versions are provided during intialization"""
        # Generate a SMIRNOFFSpecError by not providing a section version
        with pytest.raises(SMIRNOFFSpecError, match='Missing version while trying to construct') as excinfo:
            ph = ParameterHandler()
        # Successfully create ParameterHandler by skipping version check
        ph = ParameterHandler(skip_version_check=True)

        # Successfully create ParameterHandler by providing max supported version
        ph = ParameterHandler(version=ParameterHandler._MAX_SUPPORTED_SECTION_VERSION)

        # Successfully create ParameterHandler by providing min supported version
        ph = ParameterHandler(version=ParameterHandler._MIN_SUPPORTED_SECTION_VERSION)

        # Generate a SMIRNOFFSpecError ParameterHandler by providing a value higher than the max supported
        with pytest.raises(SMIRNOFFVersionError, match='SMIRNOFF offxml file was written with version 1000.0, '
                                                    'but this version of ForceField only supports version') as excinfo:
            ph = ParameterHandler(version='1000.0')


        # Generate a SMIRNOFFSpecError ParameterHandler by providing a value lower than the min supported
        with pytest.raises(SMIRNOFFVersionError, match='SMIRNOFF offxml file was written with version 0.1, '
                                                    'but this version of ForceField only supports version') as excinfo:
            ph = ParameterHandler(version='0.1')

    def test_add_delete_cosmetic_attributes(self):
        """Test ParameterHandler.to_dict() function when some parameters are in
        different units (proper behavior is to convert all quantities to the last-
        read unit)
        """
        from simtk import unit
        bh = BondHandler(skip_version_check=True)
        bh.add_parameter({'smirks': '[*:1]-[*:2]',
                          'length': 1*unit.angstrom,
                          'k': 10*unit.kilocalorie_per_mole/unit.angstrom**2})
        bh.add_parameter({'smirks': '[*:1]=[*:2]',
                          'length': 0.2*unit.nanometer,
                          'k': 0.4*unit.kilojoule_per_mole/unit.nanometer**2})

        # Ensure the cosmetic attribute is present by default during output
        bh.add_cosmetic_attribute('pilot', 'alice')
        param_dict = bh.to_dict()
        assert ('pilot', 'alice') in param_dict.items()

        # Ensure the cosmetic attribute isn't present if we request that it be discarded
        param_dict = bh.to_dict(discard_cosmetic_attributes=True)
        assert ('pilot', 'alice') not in param_dict.items()

        # Manually delete the cosmetic attribute and ensure it doesn't get written out
        bh.delete_cosmetic_attribute('pilot')
        param_dict = bh.to_dict()
        assert ('pilot', 'alice') not in param_dict.items()



class TestParameterList:
    """Test capabilities of ParameterList for accessing and manipulating SMIRNOFF parameter definitions.
    """

    def test_create(self):
        """Test creation of a parameter list.
        """
        p1 = ParameterType(smirks='[*:1]')
        p2 = ParameterType(smirks='[#1:1]')
        parameters = ParameterList([p1, p2])

    @pytest.mark.wip(reason="Until ChemicalEnvironment won't be refactored to use the ToolkitRegistry "
                            "API, the smirks assignment will fail with RDKit.")
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

    def test_index(self):
        """
        Tests the ParameterList.index() function by attempting lookups by SMIRKS and by ParameterType equivalence.

        """
        p1 = ParameterType(smirks='[*:1]')
        p2 = ParameterType(smirks='[#1:1]')
        p3 = ParameterType(smirks='[#7:1]')
        parameters = ParameterList([p1, p2, p3])
        assert parameters.index(p1) == 0
        assert parameters.index(p2) == 1
        assert parameters.index(p3) == 2
        assert parameters.index('[*:1]') == 0
        assert parameters.index('[#1:1]') == 1
        assert parameters.index('[#7:1]') == 2
        with pytest.raises(IndexError, match=r'SMIRKS \[#2:1\] not found in ParameterList') as excinfo:
            parameters.index('[#2:1]')

        p4 = ParameterType(smirks='[#2:1]')
        with pytest.raises(ValueError, match='is not in list') as excinfo:
            parameters.index(p4)

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

        with pytest.raises(IndexError, match='list assignment index out of range'):
            del parameters[4]
        with pytest.raises(IndexError, match=r'SMIRKS \[#6:1\] not found in ParameterList'):
            del parameters['[#6:1]']

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


    def test_append(self):
        """
        Test ParameterList.append, ensuring that the new parameter was added to the bottom of the list
        and that it is properly recorded as the most recently-added.
        """
        p1 = ParameterType(smirks='[*:1]-[*:2]')
        p2 = ParameterType(smirks='[*:1]=[*:2]')
        param_list = ParameterList()
        param_list.append(p1)
        assert len(param_list) == 1
        assert '[*:1]-[*:2]' in param_list
        param_list.append(p2)
        assert len(param_list) == 2
        assert '[*:1]=[*:2]' in param_list
        assert param_list[-1] == p2

    def test_insert(self):
        """
        Test ParameterList.insert, ensuring that the new parameter was added to the proper spot in
        the list and that it is propertly recorded as the most recently added.
        """
        p1 = ParameterType(smirks='[*:1]-[*:2]')
        p2 = ParameterType(smirks='[*:1]=[*:2]')
        p3 = ParameterType(smirks='[*:1]#[*:2]')
        param_list = ParameterList([p1, p2])
        param_list.insert(1, p3)
        assert param_list[1] == p3

    def test_extend(self):
        """
        Test ParameterList.extend, ensuring that the new parameter was added to the proper spot in
        the list and that it is propertly recorded as the most recently added.
        """
        p1 = ParameterType(smirks='[*:1]-[*:2]')
        p2 = ParameterType(smirks='[*:1]=[*:2]')
        param_list1 = ParameterList()
        param_list2 = ParameterList([p1, p2])

        param_list1.extend(param_list2)
        assert len(param_list1) == 2
        assert '[*:1]-[*:2]' in param_list1
        assert '[*:1]=[*:2]' in param_list1
        assert param_list1[-1] == p2

    def test_to_list(self):
        """Test basic ParameterList.to_list() function, ensuring units are preserved"""
        from simtk import unit
        p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                  length=1.01 * unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        p2 = BondHandler.BondType(smirks='[*:1]=[*:2]',
                                  length=1.02 * unit.angstrom,
                                  k=6 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        p3 = BondHandler.BondType(smirks='[*:1]#[*:3]',
                                  length=1.03 * unit.angstrom,
                                  k=7 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        parameter_list = ParameterList([p1, p2, p3])
        ser_param_list = parameter_list.to_list()
        assert len(ser_param_list) == 3
        assert ser_param_list[0]['length'] == 1.01 * unit.angstrom

    def test_round_trip(self):
        """Test basic ParameterList.to_list() function and constructor"""
        from simtk import unit
        p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                  length=1.01 * unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        p2 = BondHandler.BondType(smirks='[*:1]=[*:2]',
                                  length=1.02 * unit.angstrom,
                                  k=6 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        p3 = BondHandler.BondType(smirks='[*:1]#[*:3]',
                                  length=1.03 * unit.angstrom,
                                  k=7 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        parameter_list = ParameterList([p1, p2, p3])
        param_dict_list = parameter_list.to_list()
        parameter_list_2 = ParameterList()
        for param_dict in param_dict_list:
            new_parameter = BondHandler.BondType(**param_dict)
            parameter_list_2.append(new_parameter)
        assert parameter_list.to_list() == parameter_list_2.to_list()


class TestParameterType:

    def test_base_parametertype_to_dict(self):
        """
        Test ParameterType to_dict.
        """
        p1 = ParameterType(smirks='[*:1]')
        param_dict = p1.to_dict()
        assert param_dict['smirks'] == '[*:1]'
        assert len(param_dict.keys()) == 1

    def test_bondtype_to_dict(self):
        """
        Test BondType to_dict.
        """
        from simtk import unit

        p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                  length=1.02 * unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        param_dict = p1.to_dict()
        param_dict_unitless, attached_units = detach_units(param_dict)
        assert param_dict_unitless == {'smirks': '[*:1]-[*:2]',
                                       'length': 1.02,
                                       'k': 5,}
        assert attached_units == {'length_unit': unit.angstrom,
                                  'k_unit': (unit.angstrom ** -2) * (unit.mole ** -1) * (unit.kilocalorie ** 1)
                                  }


    def test_bondtype_to_dict_custom_output_units(self):
        """
        Test BondType to_dict with custom output units.
        """
        from simtk import unit
        p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        param_dict = p1.to_dict()
        param_dict_unitless, attached_units=detach_units(param_dict, output_units={'length_unit':
                                                                                       unit.nanometer})
        assert attached_units['length_unit'] == unit.nanometer
        assert abs(param_dict_unitless['length'] - 0.102) < 1e-10


    def test_bondtype_to_dict_invalid_output_units(self):
        """
        Test ParameterType to_dict with invalid output units.
        """
        from simtk import unit
        p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2
                                  )
        param_dict = p1.to_dict()
        with pytest.raises(ValueError,
                           match='Requested output unit calorie is not compatible with quantity unit angstrom .'
                           ) as context:
            param_dict_unitless, attached_units = detach_units(param_dict, output_units = {'length_unit':
                                                                                               unit.calorie})


    def test_read_write_optional_parameter_attribute(self):
        """
        Test ParameterTypes' ability to store and write out optional attributes passed to __init__()
        """
        from simtk import unit

        p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2,
                                  id='b1'
                                  )
        param_dict= p1.to_dict()
        assert ('id', 'b1') in param_dict.items()

    def test_read_write_cosmetic_parameter_attribute(self):
        """
        Test ParameterTypes' ability to store and write out cosmetic attributes passed to __init__()
        """
        from simtk import unit

        p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2,
                                  pilot='alice',
                                  allow_cosmetic_attributes=True
                                  )
        param_dict= p1.to_dict(discard_cosmetic_attributes=False)
        assert ('pilot', 'alice') in param_dict.items()

    def test_read_but_dont_write_cosmetic_parameter_attribute(self):
        """
        Test ParameterTypes' ability to ignore cosmetic attributes passed to __init__() if instructed
        """
        from simtk import unit

        p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2,
                                  pilot='alice',
                                  allow_cosmetic_attributes=True
                                  )
        param_dict = p1.to_dict(discard_cosmetic_attributes=True)
        assert ('pilot', 'alice') not in param_dict.items()

    def test_error_cosmetic_parameter_attribute(self):
        """
        Test that ParameterTypes raise an error on receiving unexpected attributes passed to __init__()
        """
        from simtk import unit

        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg (pilot: alice)*") as context:
            p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                      length=1.02*unit.angstrom,
                                      k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2,
                                      pilot='alice',
                                      allow_cosmetic_attributes=False
                                      )


    def test_add_delete_cosmetic_attrib(self):
        """
        Test adding and deleting cosmetic attributes for already-initialized ParameterType objects
        """
        from simtk import unit

        p1 = BondHandler.BondType(smirks='[*:1]-[*:2]',
                                  length=1.02*unit.angstrom,
                                  k=5 * unit.kilocalorie_per_mole / unit.angstrom ** 2,
                                  )
        # Ensure the cosmetic attribute is present by default during output
        p1.add_cosmetic_attribute('pilot', 'alice')
        param_dict = p1.to_dict()
        assert ('pilot', 'alice') in param_dict.items()

        # Ensure the cosmetic attribute isn't present if we request that it be discarded
        param_dict = p1.to_dict(discard_cosmetic_attributes=True)
        assert ('pilot', 'alice') not in param_dict.items()

        # Manually delete the cosmetic attribute and ensure it doesn't get written out
        p1.delete_cosmetic_attribute('pilot')
        param_dict = p1.to_dict()
        assert ('pilot', 'alice') not in param_dict.items()




    def test_single_term_proper_torsion(self):
        """
        Test creation and serialization of a single-term proper torsion
        """
        from simtk import unit

        p1 = ProperTorsionHandler.ProperTorsionType(smirks='[*:1]-[*:2]-[*:3]-[*:4]',
                                                    phase1=30 * unit.degree,
                                                    periodicity1=2,
                                                    k1=5 * unit.kilocalorie_per_mole
                                                    )
        param_dict = p1.to_dict()
        assert ('k1', 5 * unit.kilocalorie_per_mole) in param_dict.items()
        assert ('phase1', 30 * unit.degree) in param_dict.items()
        assert ('periodicity1', 2) in param_dict.items()
        assert 'idivf' not in param_dict

    def test_single_term_proper_torsion_w_idivf(self):
        """
        Test creation and serialization of a single-term proper torsion
        """
        from simtk import unit

        p1 = ProperTorsionHandler.ProperTorsionType(smirks='[*:1]-[*:2]-[*:3]-[*:4]',
                                                    phase1=30 * unit.degree,
                                                    periodicity1=2,
                                                    k1=5 * unit.kilocalorie_per_mole,
                                                    idivf1=4
                                                    )

        param_dict = p1.to_dict()
        assert ('k1', 5 * unit.kilocalorie_per_mole) in param_dict.items()
        assert ('phase1', 30 * unit.degree) in param_dict.items()
        assert ('periodicity1', 2) in param_dict.items()
        assert ('idivf1', 4) in param_dict.items()


    def test_multi_term_proper_torsion(self):
        """
        Test creation and serialization of a multi-term proper torsion
        """
        from simtk import unit

        p1 = ProperTorsionHandler.ProperTorsionType(smirks='[*:1]-[*:2]-[*:3]-[*:4]',
                                                    phase1=30 * unit.degree,
                                                    periodicity1=2,
                                                    k1=5 * unit.kilocalorie_per_mole,
                                                    phase2=31 * unit.degree,
                                                    periodicity2=3,
                                                    k2=6 * unit.kilocalorie_per_mole,
                                                    )
        param_dict = p1.to_dict()
        assert ('k1', 5 * unit.kilocalorie_per_mole) in param_dict.items()
        assert ('phase1', 30 * unit.degree) in param_dict.items()
        assert ('periodicity1', 2) in param_dict.items()
        assert ('k2', 6 * unit.kilocalorie_per_mole) in param_dict.items()
        assert ('phase2', 31 * unit.degree) in param_dict.items()
        assert ('periodicity2', 3) in param_dict.items()

    def test_multi_term_proper_torsion_skip_index(self):
        """
        Test creation and serialization of a multi-term proper torsion where
        the indices are not consecutive and a SMIRNOFFSpecError is raised
        """
        from simtk import unit

        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg (phase3: 31 deg)*") as context:
            p1 = ProperTorsionHandler.ProperTorsionType(smirks='[*:1]-[*:2]-[*:3]-[*:4]',
                                                        phase1=30 * unit.degree,
                                                        periodicity1=2,
                                                        k1=5 * unit.kilocalorie_per_mole,
                                                        phase3=31 * unit.degree,
                                                        periodicity3=3,
                                                        k3=6 * unit.kilocalorie_per_mole,
                                                        )

    def test_multi_term_proper_torsion_bad_units(self):
        """
        Test creation and serialization of a multi-term proper torsion where
        one of the terms has incorrect units
        """
        from simtk import unit

        with pytest.raises(IncompatibleUnitError, match="__init__ function.  phase with value 31 A, is incompatible")\
                as context:
            p1 = ProperTorsionHandler.ProperTorsionType(smirks='[*:1]-[*:2]-[*:3]-[*:4]',
                                                        phase1=30 * unit.degree,
                                                        periodicity1=2,
                                                        k1=5 * unit.kilocalorie_per_mole,
                                                        phase2=31 * unit.angstrom, # This should be caught
                                                        periodicity2=3,
                                                        k2=6 * unit.kilocalorie_per_mole,
                                                        )


    def test_torsion_handler_potential_setting(self):
        """
        Test creation of TorsionHandlers with the deprecated 0.2 potential value "charmm" instead of the current
        supported potential value "fourier".
        """
        # Test creating ProperTorsionHandlers
        with pytest.raises(SMIRNOFFSpecError, match="ProperTorsionHandler given 'potential' value of 'charmm'. "
                                                    "Supported options are "
                                                    "[[]'k[*][(]1[+]cos[(]periodicity[*]theta[-]phase[)][)]'[]].")\
                as context:
            ph1 = ProperTorsionHandler(potential='charmm', skip_version_check=True)
        ph1 = ProperTorsionHandler(potential='k*(1+cos(periodicity*theta-phase))', skip_version_check=True)

        # Same test, but with ImproperTorsionHandler
        with pytest.raises(SMIRNOFFSpecError, match="ImproperTorsionHandler given 'potential' value of 'charmm'. "
                                                    "Supported options are "
                                                    "[[]'k[*][(]1[+]cos[(]periodicity[*]theta[-]phase[)][)]'[]].")\
                as context:
            ph1 = ImproperTorsionHandler(potential='charmm', skip_version_check=True)
        ph1 = ImproperTorsionHandler(potential='k*(1+cos(periodicity*theta-phase))', skip_version_check=True)


        #     p1 = ProperTorsionHandler.ProperTorsionType(smirks='[*:1]-[*:2]-[*:3]-[*:4]',
        #                                                 phase1=30 * unit.degree,
        #                                                 periodicity1=2,
        #                                                 k1=5 * unit.kilocalorie_per_mole,
        #                                                 phase2=31 * unit.angstrom, # This should be caught
        #                                                 periodicity2=3,
        #                                                 k2=6 * unit.kilocalorie_per_mole,
        #                                                 )

# TODO: test_nonbonded_settings (ensure that choices in Electrostatics and vdW tags resolve
#                                to correct openmm.NonbondedForce subtypes, that setting different cutoffs raises
#                                exceptions, etc)
# TODO: test_(all attributes of all ParameterTypes)
# TODO: test_add_parameter_fractional_bondorder
# TODO: test_get_indexed_attrib
# TODO: test_set_unitbearing_attrib (requires implementing __getattr__ and __setattr__)
# TODO: test_parametertype_unit_getattr
# TODO: test_parametertype_unit_setattr
# TODO: test_optional_attribs
# TODO: test_optional_indexed_attribs
# TODO: test_attach_units
# TODO: test_detach_units
# TODO: test_(X)handler_compatibility, where X is all handlers