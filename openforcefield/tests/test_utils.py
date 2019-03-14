#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for utility methods

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

from openforcefield import utils

import os

#=============================================================================================
# TESTS
#=============================================================================================

def test_subclasses(self):
    """Test that all subclasses (and descendents) are correctly identified by all_subclasses()"""
    class Foo(object):
        pass
    class FooSubclass1(Foo):
        pass
    class FooSubclass2(Foo):
        pass
    class FooSubSubclass(FooSubclass1):
        pass

    subclass_names = [ cls.__name__ for cls in utils.all_subclasses(Foo) ]
    assert set(subclass_names) == set(['FooSubclass1', 'FooSubclass2', 'FooSubSubclass'])

def test_temporary_cd():
    """Test temporary_cd() context manager"""
    initial_dir = os.getcwd()
    temporary_dir = '/'
    from openforcefield.utils import temporary_cd
    with temporary_cd(temporary_dir):
        assert os.getcwd() == temporary_dir
    assert os.getcwd() == initial_dir

def test_temporary_directory():
    """Test temporary_directory() context manager"""
    from openforcefield.utils import temporary_directory
    initial_dir = os.getcwd()
    with temporary_directory() as tmp_dir:
        # Make sure the temporary directory is not the current directory
        assert tmp_dir != initial_dir

        # Make sure the temporary directory is writeable
        output_filename = os.path.join(tmp_dir, 'test.out')
        outfile = open(output_filename,'w')
        outfile.write('test')
        outfile.close()

        # Make sure the file we created exists
        assert os.path.exists(output_filename)

    # Make sure the directory has been cleaned up
    assert not os.path.exists(tmp_dir), "Temporary directory was not automatically cleaned up."
    assert not os.path.exists(output_filename), "Temporary directory was not automatically cleaned up."

def test_get_data_filename():
    """Test get_data_filename()"""
    from openforcefield.utils import get_data_filename
    filename = get_data_filename('forcefield/tip3p.offxml')
    assert os.path.exists(filename)

def test_one_dimension_serialize_quantity():
    """Test utils.serialize_quantity's behavior for a single-dimension quantity"""
    from simtk import unit
    quantity =unit.Quantity(5, unit.angstrom)
    serialized = utils.serialize_quantity(quantity)
    assert serialized == {'unitless_value': 5, 'unit': [('angstrom', 1)]}

def test_two_dimension_serialize_quantity():
    """Test utils.serialize_quantity's behavior for a single-dimension quantity"""
    from simtk import unit
    quantity = unit.Quantity(5, unit.angstrom * unit.kilocalorie_per_mole)
    serialized = utils.serialize_quantity(quantity)
    assert serialized['unitless_value'] ==  5
    assert ('angstrom', 1) in serialized['unit']
    assert ('mole', -1) in serialized['unit']
    assert ('kilocalorie', 1) in serialized['unit']

def test_two_dimension_serialize_quantity_custom_output():
    """Test utils.serialize_quantity's behavior for a single-dimension quantity"""
    from simtk import unit
    quantity = unit.Quantity(5, unit.angstrom * unit.kilocalorie_per_mole)
    serialized = utils.serialize_quantity(quantity,
                                          # output as eV*nm/mol
                                          output_unit=unit.nanometer*unit.elementary_charge*unit.volt/unit.mole)
    assert serialized['unitless_value'] == 1.3057238144187346e+22
    assert ('nanometer', 1) in serialized['unit']
    assert ('mole', -1) in serialized['unit']
    assert ('elementary charge', 1) in serialized['unit']
    assert ('volt', 1) in serialized['unit']


@pytest.mark.parametrize('unit_string,expected_unit',[
    ('kilocalories_per_mole', unit.kilocalories_per_mole),
    ('kilocalories_per_mole/angstrom**2', unit.kilocalories_per_mole/unit.angstrom**2),
    ('joule/(mole * nanometer**2)', unit.joule/(unit.mole * unit.nanometer**2)),
    ('picosecond**(-1)', unit.picosecond**(-1)),
    ('300.0 * kelvin', 300*unit.kelvin),
    ('1 * kilojoule + 500 * joule', 1.5*unit.kilojoule),
    ('1 / meter', 1.0 / unit.meter)
])
def test_ast_unit_eval(unit_string, expected_unit):
    """Test that _ast_unit_eval() correctly parses string quantities."""
    ast_root_node = ast.parse(unit_string, mode='eval').body
    parsed_units = _ast_unit_eval(ast_root_node)
    assert parsed_units == expected_unit


@pytest.mark.parametrize('attributes,expected',[
    ({'not_parsed': 'blabla'},
        {}),
    ({'not_parsed': 'blabla', 'attr_unit': 'angstrom/femtosecond'},
        {'attr': unit.angstrom/unit.femtosecond}),
    ({'not_parsed': 'blabla', 'attr1_unit': 'meter', 'attr2_unit': 'kilojoule_per_mole'},
        {'attr1': unit.meter, 'attr2': unit.kilojoule_per_mole}),
])
def test_extract_attached_units(attributes, expected):
    """Test that _extract_attached_units() correctly parses the correct."""
    assert _extract_attached_units(attributes) == expected


@pytest.mark.parametrize('attributes',[
    {'attr_unit': '300.0 * kelvin'},
    {'attr_unit': '1 / picosecond'}
])
def test_extract_attached_units_raises(attributes):
    """Test that _extract_attached_units() raises an error when a quantity is specified instead of a unit."""
    with pytest.raises(ValueError, match='associated to a quantity rather than only units'):
        _extract_attached_units(attributes)
