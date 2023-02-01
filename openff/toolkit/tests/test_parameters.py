"""
Test classes and function in module openff.toolkit.typing.engines.smirnoff.parameters.

"""
from collections import defaultdict
from inspect import isabstract, isclass
from typing import Dict, List, Tuple

import numpy
import pytest
from numpy.testing import assert_almost_equal
from openff.units import unit
from packaging.version import Version

import openff.toolkit.typing.engines.smirnoff.parameters
from openff.toolkit.tests.mocking import VirtualSiteMocking
from openff.toolkit.tests.utils import does_not_raise
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.typing.engines.smirnoff.parameters import (
    BondHandler,
    ChargeIncrementModelHandler,
    ElectrostaticsHandler,
    GBSAHandler,
    ImproperTorsionHandler,
    IndexedParameterAttribute,
    LibraryChargeHandler,
    ParameterAttribute,
    ParameterHandler,
    ParameterList,
    ParameterType,
    ProperTorsionHandler,
    VirtualSiteHandler,
    _cal_mol_a2,
    _linear_inter_or_extrapolate,
    _ParameterAttributeHandler,
    vdWHandler,
)
from openff.toolkit.utils.collections import ValidatedList
from openff.toolkit.utils.exceptions import (
    DuplicateParameterError,
    IncompatibleParameterError,
    IncompatibleUnitError,
    MissingIndexedAttributeError,
    NotEnoughPointsForInterpolationError,
    ParameterLookupError,
    SMIRNOFFSpecError,
    SMIRNOFFSpecUnimplementedError,
    SMIRNOFFVersionError,
)


class TestParameterAttribute:
    """Test cases for the descriptor ParameterAttribute."""

    def test_default_value(self):
        """Default values are assigned correctly on initialization."""

        class MyParameter:
            attr_optional = ParameterAttribute(default=2)

        my_par = MyParameter()
        assert my_par.attr_optional == 2

    def test_none_default_value(self):
        """None is a valid default value for ParameterAttribute."""

        class MyParameter:
            attr_optional = ParameterAttribute(default=None)

        my_par = MyParameter()
        assert my_par.attr_optional is None

    def test_required_value(self):
        """AttributeError is raised if a required attribute is accessed before being initialized."""

        class MyParameter:
            attr_required = ParameterAttribute()

        my_par = MyParameter()
        with pytest.raises(AttributeError):
            my_par.attr_required

    def test_unit_validation(self):
        """ParameterAttributes attached to a unit are validated correctly."""

        class MyParameter:
            attr_unit = ParameterAttribute(
                unit=unit.kilocalorie / unit.mole / unit.angstrom**2
            )

        my_par = MyParameter()

        # TypeError is raised when setting a unit-less value.
        with pytest.raises(IncompatibleUnitError, match="should have units of"):
            my_par.attr_unit = 3.0
        # TypeError is raised when setting a value with incorrect units.
        with pytest.raises(IncompatibleUnitError, match="should have units of"):
            my_par.attr_unit = 3.0 * unit.kilocalorie / unit.mole
        # Otherwise the attribute is assigned correctly.
        value = 3.0 * unit.kilocalorie / unit.mole / unit.angstrom**2
        my_par.attr_unit = value
        assert my_par.attr_unit == value
        assert my_par.attr_unit.units == value.units

    def test_quantity_string_parsing(self):
        """ParameterAttributes attached to units convert strings into Quantity objects."""

        class MyParameter:
            attr_unit = ParameterAttribute(unit=unit.meter / unit.second**2)

        my_par = MyParameter()

        my_par.attr_unit = "3.0*meter/second**2"
        assert my_par.attr_unit == 3.0 * unit.meter / unit.second**2
        assert my_par.attr_unit.units == unit.meter / unit.second**2

        # Assigning incorrect units still raises an error.
        with pytest.raises(IncompatibleUnitError, match="should have units of"):
            my_par.attr_unit = "3.0"
        with pytest.raises(IncompatibleUnitError, match="should have units of"):
            my_par.attr_unit = "3.0*meter/second"

    def test_custom_converter(self):
        """Custom converters of ParameterAttributes are executed correctly."""

        class MyParameter:
            attr_all_to_float = ParameterAttribute(converter=float)
            attr_int_to_float = ParameterAttribute()

            @attr_int_to_float.converter
            def attr_int_to_float(self, attr, value):
                """Convert only integers to float"""
                if isinstance(value, int):
                    return float(value)
                elif not isinstance(value, float):
                    raise TypeError()
                return value

        my_par = MyParameter()

        # Both strings and integers are converted to floats when casted with float().
        my_par.attr_all_to_float = "1.0"
        assert (
            isinstance(my_par.attr_all_to_float, float)
            and my_par.attr_all_to_float == 1.0
        )
        my_par.attr_all_to_float = 2
        assert (
            isinstance(my_par.attr_all_to_float, float)
            and my_par.attr_all_to_float == 2.0
        )

        # Only integers are converted with the custom converter function.
        with pytest.raises(TypeError):
            my_par.attr_int_to_float = "1.0"
        my_par.attr_int_to_float = 2
        assert (
            isinstance(my_par.attr_int_to_float, float)
            and my_par.attr_int_to_float == 2.0
        )

    def test_default_pass_validation(self):
        """The default value of ParameterAttribute is always allowed regardless of the validator/converter."""

        class MyParameter:
            attr = ParameterAttribute(
                default=None, unit=unit.angstrom, converter=unit.Quantity
            )

        my_par = MyParameter()
        my_par.attr = 3.0 * unit.nanometer
        my_par.attr = None
        assert my_par.attr is None

    def test_get_descriptor_object(self):
        """When the descriptor is called from the class, the ParameterAttribute descriptor is returned."""

        class MyParameter:
            attr = ParameterAttribute()

        assert isinstance(MyParameter.attr, ParameterAttribute)


class TestIndexedParameterAttribute:
    """Tests for the IndexedParameterAttribute descriptor."""

    def test_tuple_conversion(self):
        """IndexedParameterAttribute converts internally sequences to ValidatedList."""

        class MyParameter:
            attr_indexed = IndexedParameterAttribute()

        my_par = MyParameter()
        my_par.attr_indexed = [1, 2, 3]
        assert isinstance(my_par.attr_indexed, ValidatedList)

    def test_indexed_default(self):
        """IndexedParameterAttribute handles default values correctly."""

        class MyParameter:
            attr_indexed_optional = IndexedParameterAttribute(default=None)

        my_par = MyParameter()
        assert my_par.attr_indexed_optional is None

        # Assigning the default is allowed.
        my_par.attr_indexed_optional = None
        assert my_par.attr_indexed_optional is None

    def test_units_on_all_elements(self):
        """IndexedParameterAttribute validates every single element of the sequence."""

        class MyParameter:
            attr_indexed_unit = IndexedParameterAttribute(unit=unit.gram)

        my_par = MyParameter()

        # Strings are correctly converted.
        my_par.attr_indexed_unit = ["1.0*gram", 2 * unit.gram]
        assert my_par.attr_indexed_unit == [1.0 * unit.gram, 2 * unit.gram]

        # Incompatible units on a single elements are correctly caught.
        with pytest.raises(IncompatibleUnitError, match="should have units of"):
            my_par.attr_indexed_unit = [3.0, 2 * unit.gram]
        with pytest.raises(IncompatibleUnitError, match="should have units of"):
            my_par.attr_indexed_unit = [2 * unit.gram, 4.0 * unit.meter]

    def test_converter_on_all_elements(self):
        """IndexedParameterAttribute calls custom converters on every single element of the sequence."""

        class MyParameter:
            attr_indexed_converter = IndexedParameterAttribute(converter=float)

        my_par = MyParameter()

        my_par.attr_indexed_converter = [1, "2.0", "1e-3", 4.0]
        assert my_par.attr_indexed_converter == [1.0, 2.0, 1e-3, 4.0]

    def test_validate_new_elements(self):
        """New elements set in the list are correctly validated."""

        class MyParameter:
            attr_indexed = IndexedParameterAttribute(converter=int)

        my_par = MyParameter()
        my_par.attr_indexed = (1, 2, 3)

        # Modifying one or more elements of the list should validate them.
        my_par.attr_indexed[2] = "4"
        assert my_par.attr_indexed[2] == 4
        my_par.attr_indexed[0:3] = ["2", "3", 4]
        assert my_par.attr_indexed == [2, 3, 4]

        # Same for append().
        my_par.attr_indexed.append("5")
        assert my_par.attr_indexed[3] == 5

        # And extend.
        my_par.attr_indexed.extend([6, "7"])
        assert my_par.attr_indexed[5] == 7
        my_par.attr_indexed += ["8", 9]
        assert my_par.attr_indexed[6] == 8

        # And insert.
        my_par.attr_indexed.insert(5, "10")
        assert my_par.attr_indexed[5] == 10


class TestInterpolation:
    """Test method(s) that are used for functionality like fractional bond order-dependent parameter interpolation"""

    @pytest.mark.parametrize(
        ("fractional_bond_order", "k_interpolated"),
        [(1.6, 1.48), (0.7, 0.76), (2.3, 2.04)],
    )
    def test_linear_inter_or_extrapolate(self, fractional_bond_order, k_interpolated):
        """Test that linear interpolation works as expected"""
        k_bondorder = {
            1: 1 * unit.kilocalorie / unit.mole,
            2: 1.8 * unit.kilocalorie / unit.mole,
        }

        k = _linear_inter_or_extrapolate(k_bondorder, fractional_bond_order)
        assert_almost_equal(k.m, k_interpolated)

    def test_linear_inter_or_extrapolate_one_point(self):
        """Test that linear interpolation raises an error if attempted with just one point"""
        k_bondorder = {
            2: 1.8 * unit.kilocalorie / unit.mole,
        }
        with pytest.raises(NotEnoughPointsForInterpolationError):
            _linear_inter_or_extrapolate(k_bondorder, 1)

    @pytest.mark.parametrize(
        ("fractional_bond_order", "k_interpolated"),
        [(1.6, 1.48), (0.7, 0.76), (2.3, 2.01), (3.1, 2.57)],
    )
    def test_linear_inter_or_extrapolate_3_terms(
        self, fractional_bond_order, k_interpolated
    ):
        """Test that linear interpolation works as expected for three terms"""
        k_bondorder = {
            1: 1 * unit.kilocalorie / unit.mole,
            2: 1.8 * unit.kilocalorie / unit.mole,
            3: 2.5 * unit.kilocalorie / unit.mole,
        }

        k = _linear_inter_or_extrapolate(k_bondorder, fractional_bond_order)
        assert_almost_equal(k.m, k_interpolated)

    def test_linear_inter_or_extrapolate_below_zero(self):
        """Test that linear interpolation does not error if resulting k less than 0"""
        k_bondorder = {
            1: 1 * unit.kilocalorie / unit.mole,
            2: 2.3 * unit.kilocalorie / unit.mole,
        }

        fractional_bond_order = 0.2
        k = _linear_inter_or_extrapolate(k_bondorder, fractional_bond_order)

        assert k.m < 0


class TestParameterAttributeHandler:
    """Test suite for the base class _ParameterAttributeHandler."""

    def test_access_get_set_single_indexed_attribute_legacy(self):
        """Single indexed attributes such as k1 can be accessed through normal attribute syntax."""

        class MyParameterType(_ParameterAttributeHandler):
            k = IndexedParameterAttribute()

        my_parameter = MyParameterType(k=[1, 2, 3])

        # Getting the attribute works.
        assert my_parameter.k1 == 1
        assert my_parameter.k2 == 2
        assert my_parameter.k3 == 3

        # So does setting it.
        my_parameter.k2 = 5
        assert my_parameter.k2 == 5
        assert my_parameter.k == [1, 5, 3]

        # Accessing k4 raises an index error.
        with pytest.raises(
            IndexError, match="'k4' is out of bounds for indexed attribute 'k'"
        ):
            my_parameter.k4
        with pytest.raises(
            IndexError, match="'k4' is out of bounds for indexed attribute 'k'"
        ):
            my_parameter.k4 = 2

        # For other attributes, the behavior is normal.
        with pytest.raises(AttributeError, match="has no attribute 'x'"):
            my_parameter.x
        # Monkey-patching.
        my_parameter.x = 3

    def test_access_get_set_single_indexed_attribute(self):
        """Single indexed attributes such as k1 can be accessed through normal attribute syntax."""

        class MyParameterType(_ParameterAttributeHandler):
            k = IndexedParameterAttribute()

        my_parameter = MyParameterType(k=[1, 2, 3])

        # Getting the attribute works.
        assert my_parameter.k1 == 1
        assert my_parameter.k2 == 2
        assert my_parameter.k3 == 3

        # So does setting it.
        my_parameter.k2 = 5
        assert my_parameter.k2 == 5
        assert my_parameter.k == [1, 5, 3]

        # Accessing k4 raises an index error.
        with pytest.raises(
            MissingIndexedAttributeError,
            match="'k4' is out of bounds for indexed attribute 'k'",
        ):
            my_parameter.k4
        with pytest.raises(
            MissingIndexedAttributeError,
            match="'k4' is out of bounds for indexed attribute 'k'",
        ):
            my_parameter.k4 = 2

        # For other attributes, the behavior is normal.
        with pytest.raises(AttributeError, match="has no attribute 'x'"):
            my_parameter.x
        # Monkey-patching.
        my_parameter.x = 3

    def test_hasattr(self):
        """Single indexed attributes such as k1 can be accessed through normal attribute syntax."""

        class MyParameterType(_ParameterAttributeHandler):
            k = IndexedParameterAttribute()

        my_parameter = MyParameterType(k=[1, 2, 3])
        assert hasattr(my_parameter, "k3")
        assert not hasattr(my_parameter, "k4")

    def test_mro_access_get_set_single_indexed_attribute(self):
        """Attribute access is forwarded correctly to the next MRO classes."""

        class MixIn:
            """Utility class to keep track of whether __get/setattr__ are called."""

            data = {}

            def __getattr__(self, item):
                self.getattr_flag = True
                try:
                    return self.data[item]
                except KeyError:
                    raise AttributeError()

            def __setattr__(self, key, value):
                self.data[key] = value
                super().__setattr__("setattr_flag", True)

            def assert_getattr(self):
                assert self.getattr_flag is True
                self.getattr_flag = False

            def assert_setattr(self):
                assert self.setattr_flag is True
                super().__setattr__("setattr_flag", False)

        class MyParameterType(_ParameterAttributeHandler, MixIn):
            k = IndexedParameterAttribute()

        my_parameter = MyParameterType(k=[1, 2, 3])

        # Non-existing parameters.
        my_parameter.a = 2
        my_parameter.assert_setattr()
        my_parameter.a1 = 4
        my_parameter.assert_setattr()

        my_parameter.a
        my_parameter.assert_getattr()
        my_parameter.a1
        my_parameter.assert_getattr()


class TestParameterHandler:
    length = 1 * unit.angstrom
    k = 10 * unit.kilocalorie / unit.mole / unit.angstrom**2

    def test_tagname(self):
        """Test the TAGNAME getter and default behavior"""
        ph = ParameterHandler(skip_version_check=True)
        assert ph.TAGNAME is None

        bh = BondHandler(skip_version_check=True)
        assert bh.TAGNAME == "Bonds"

    def test_add_parameter(self):
        """Test the behavior of add_parameter"""
        bh = BondHandler(skip_version_check=True)
        param1 = {
            "smirks": "[*:1]-[*:2]",
            "length": self.length,
            "k": self.k,
            "id": "b1",
        }
        param2 = {
            "smirks": "[*:1]=[*:2]",
            "length": self.length,
            "k": self.k,
            "id": "b2",
        }
        param3 = {
            "smirks": "[*:1]#[*:2]",
            "length": self.length,
            "k": self.k,
            "id": "b3",
        }

        bh.add_parameter(param1)
        bh.add_parameter(param2)
        bh.add_parameter(param3)

        assert [p.id for p in bh._parameters] == ["b1", "b2", "b3"]

        param_duplicate_smirks = {
            "smirks": param2["smirks"],
            "length": 2 * self.length,
            "k": 2 * self.k,
        }

        # Ensure a duplicate parameter cannot be added
        with pytest.raises(DuplicateParameterError):
            bh.add_parameter(param_duplicate_smirks)

        dict_to_add_by_smirks = {
            "smirks": "[#1:1]-[#6:2]",
            "length": self.length,
            "k": self.k,
            "id": "d1",
        }
        dict_to_add_by_index = {
            "smirks": "[#1:1]-[#8:2]",
            "length": self.length,
            "k": self.k,
            "id": "d2",
        }

        param_to_add_by_smirks = BondHandler.BondType(
            **{
                "smirks": "[#6:1]-[#6:2]",
                "length": self.length,
                "k": self.k,
                "id": "p1",
            }
        )
        param_to_add_by_index = BondHandler.BondType(
            **{
                "smirks": "[#6:1]=[#8:2]",
                "length": self.length,
                "k": self.k,
                "id": "p2",
            }
        )

        param_several_apart = {
            "smirks": "[#1:1]-[#7:2]",
            "length": self.length,
            "k": self.k,
            "id": "s0",
        }

        # The `before` parameter should come after the `after` parameter
        # in the parameter list; i.e. in this list of ['-', '=', '#'], it is
        # impossible to add a new parameter after '=' *and* before '-'
        with pytest.raises(ValueError):
            # Test invalid parameter order by SMIRKS
            bh.add_parameter(
                dict_to_add_by_smirks, after="[*:1]=[*:2]", before="[*:1]-[*:2]"
            )

        with pytest.raises(ValueError):
            # Test invalid parameter order by index
            bh.add_parameter(dict_to_add_by_index, after=1, before=0)

        # Add d1 before param b2
        bh.add_parameter(dict_to_add_by_smirks, before="[*:1]=[*:2]")

        assert [p.id for p in bh._parameters] == ["b1", "d1", "b2", "b3"]

        # Add d2 after index 2 (which is also param b2)
        bh.add_parameter(dict_to_add_by_index, after=2)

        assert [p.id for p in bh._parameters] == ["b1", "d1", "b2", "d2", "b3"]

        # Add p1 before param b3
        bh.add_parameter(parameter=param_to_add_by_smirks, before="[*:1]=[*:2]")

        assert [p.id for p in bh._parameters] == ["b1", "d1", "p1", "b2", "d2", "b3"]

        # Add p2 after index 2 (which is param p1)
        bh.add_parameter(parameter=param_to_add_by_index, after=2)

        assert [p.id for p in bh._parameters] == [
            "b1",
            "d1",
            "p1",
            "p2",
            "b2",
            "d2",
            "b3",
        ]

        # Add s0 between params that are several positions apart
        bh.add_parameter(param_several_apart, after=1, before=6)

        assert [p.id for p in bh._parameters] == [
            "b1",
            "d1",
            "s0",
            "p1",
            "p2",
            "b2",
            "d2",
            "b3",
        ]

    def test_different_units_to_dict(self):
        """Test ParameterHandler.to_dict() function when some parameters are in
        different units (proper behavior is to convert all quantities to the last-
        read unit)
        """
        bh = BondHandler(skip_version_check=True)
        bh.add_parameter(
            {
                "smirks": "[*:1]-[*:2]",
                "length": 1 * unit.angstrom,
                "k": 10 * unit.kilocalorie / unit.mole / unit.angstrom**2,
            }
        )
        bh.add_parameter(
            {
                "smirks": "[*:1]=[*:2]",
                "length": 0.2 * unit.nanometer,
                "k": 0.4 * unit.kilojoule / unit.mole / unit.nanometer**2,
            }
        )
        bh_dict = bh.to_dict()
        assert bh_dict["Bond"][0]["length"] == unit.Quantity(
            value=1, units=unit.angstrom
        )
        assert bh_dict["Bond"][1]["length"] == unit.Quantity(
            value=2, units=unit.angstrom
        )

    def test_to_dict_maintain_units(self):
        """Test ParameterHandler.to_dict() function when parameters were provided in different units"""
        bh = BondHandler(skip_version_check=True)
        bh.add_parameter(
            {
                "smirks": "[*:1]-[*:2]",
                "length": 1 * unit.angstrom,
                "k": 10 * unit.kilocalorie / unit.mole / unit.angstrom**2,
            }
        )
        bh.add_parameter(
            {
                "smirks": "[*:1]=[*:2]",
                "length": 0.2 * unit.nanometer,
                "k": 0.4 * unit.kilojoule / unit.mole / unit.nanometer**2,
            }
        )
        bh_dict = bh.to_dict()
        assert bh_dict["Bond"][0]["length"] == unit.Quantity(1.0, unit.angstrom)
        assert bh_dict["Bond"][0]["length"].units == unit.angstrom
        assert bh_dict["Bond"][1]["length"] == unit.Quantity(0.2, unit.nanometer)
        assert bh_dict["Bond"][1]["length"].units == unit.nanometer

    def test_missing_section_version(self):
        """Test that exceptions are raised if invalid or improper section versions are provided during intialization"""
        # Generate a SMIRNOFFSpecError by not providing a section version
        with pytest.raises(
            SMIRNOFFSpecError, match="Missing version while trying to construct"
        ):
            ParameterHandler()

        # Successfully create ParameterHandler by skipping version check
        ParameterHandler(skip_version_check=True)

        # Successfully create ParameterHandler by providing max supported version
        ParameterHandler(version=ParameterHandler._MAX_SUPPORTED_SECTION_VERSION)

        # Successfully create ParameterHandler by providing min supported version
        ParameterHandler(version=ParameterHandler._MIN_SUPPORTED_SECTION_VERSION)

        # Generate a SMIRNOFFSpecError by providing a value higher than the max supported
        with pytest.raises(
            SMIRNOFFVersionError,
            match="SMIRNOFF offxml file was written with version 1000.0, "
            "but this version of ForceField only supports version",
        ):
            ParameterHandler(version="1000.0")

        # Generate a SMIRNOFFSpecError by providing a value lower than the min supported
        with pytest.raises(
            SMIRNOFFVersionError,
            match="SMIRNOFF offxml file was written with version 0.1, "
            "but this version of ForceField only supports version",
        ):
            ParameterHandler(version="0.1")

    def test_supported_version_range(self):
        """
        Ensure that version values in various formats can be correctly parsed and validated
        """

        class MyPHSubclass(ParameterHandler):
            _MIN_SUPPORTED_SECTION_VERSION = Version("0.3")
            _MAX_SUPPORTED_SECTION_VERSION = Version("2")

        with pytest.raises(SMIRNOFFVersionError):
            MyPHSubclass(version=0.1)
        with pytest.raises(Exception, match="Could not convert .*list"):
            MyPHSubclass(version=[0])
        MyPHSubclass(version=0.3)
        MyPHSubclass(version=1)
        MyPHSubclass(version="1.9")
        MyPHSubclass(version=2.0)
        with pytest.raises(SMIRNOFFVersionError):
            MyPHSubclass(version=2.1)

    def test_write_same_version_as_was_set(self):
        """Ensure that a ParameterHandler remembers the version that was set when it was initialized."""

        class MyPHSubclass(ParameterHandler):
            _MIN_SUPPORTED_SECTION_VERSION = Version("0.3")
            _MAX_SUPPORTED_SECTION_VERSION = Version("2")

        my_ph = MyPHSubclass(version=1.234)
        assert my_ph.to_dict()["version"] == str(Version("1.234"))

    def test_add_delete_cosmetic_attributes(self):
        """Test ParameterHandler.to_dict() function when some parameters are in
        different units (proper behavior is to convert all quantities to the last-
        read unit)
        """
        bh = BondHandler(skip_version_check=True)
        bh.add_parameter(
            {
                "smirks": "[*:1]-[*:2]",
                "length": 1 * unit.angstrom,
                "k": 10 * unit.kilocalorie / unit.mole / unit.angstrom**2,
            }
        )
        bh.add_parameter(
            {
                "smirks": "[*:1]=[*:2]",
                "length": 0.2 * unit.nanometer,
                "k": 0.4 * unit.kilojoule / unit.mole / unit.nanometer**2,
            }
        )

        assert not (bh.attribute_is_cosmetic("pilot"))

        # Ensure the cosmetic attribute is present by default during output
        bh.add_cosmetic_attribute("pilot", "alice")
        param_dict = bh.to_dict()
        assert ("pilot", "alice") in param_dict.items()
        assert bh.attribute_is_cosmetic("pilot")

        # Ensure the cosmetic attribute isn't present if we request that it be discarded
        param_dict = bh.to_dict(discard_cosmetic_attributes=True)
        assert "pilot" not in param_dict

        # Manually delete the cosmetic attribute and ensure it doesn't get written out
        bh.delete_cosmetic_attribute("pilot")
        param_dict = bh.to_dict()
        assert "pilot" not in param_dict
        assert not (bh.attribute_is_cosmetic("pilot"))

    def test_get_parameter(self):
        """Test that ParameterHandler.get_parameter can lookup function"""
        bh = BondHandler(skip_version_check=True, allow_cosmetic_attributes=True)

        bh.add_parameter(
            {
                "smirks": "[*:1]-[*:2]",
                "length": 1 * unit.angstrom,
                "k": 10 * unit.kilocalorie / unit.mole / unit.angstrom**2,
                "id": "b0",
            }
        )
        bh.parameters[0].add_cosmetic_attribute("foo", "bar")

        # Check base behavior
        params = bh.get_parameter({"smirks": "[*:1]-[*:2]"})

        assert params[0].length == unit.Quantity(1.0, unit.angstrom)
        assert params[0].k == unit.Quantity(
            10.0, unit.kilocalorie / unit.mole / unit.angstrom**2
        )

        # Ensure a query with no matches returns an empty list
        assert not bh.get_parameter({"smirks": "xyz"})

        # Ensure searching for a nonexistent attr does not raise an exception
        assert not bh.get_parameter({"bAdAttR": "0"})

        # Check for optional and cosmetic attrs
        optional_params = bh.get_parameter({"id": "b0"})
        cosmetic_params = bh.get_parameter({"foo": "bar"})

        assert optional_params[0].id == "b0"
        assert cosmetic_params[0]._foo == "bar"

        # Ensure selection behaves a "OR" not "AND"
        bh.add_parameter(
            {
                "smirks": "[#1:1]-[#6:2]",
                "length": 1 * unit.angstrom,
                "k": 10 * unit.kilocalorie / unit.mole / unit.angstrom**2,
                "id": "b1",
            }
        )

        params = bh.get_parameter({"id": "b0", "smirks": "[#1:1]-[#6:2]"})

        assert "b0" in [param.id for param in params]
        assert "[*:1]-[*:2]" in [param.smirks for param in params]

        # Ensure selection does not return duplicates if multiple matches
        params = bh.get_parameter({"id": "b1", "smirks": "[#1:1]-[#6:2]"})

        assert len(params) == 1

    def test_create_force(self):
        class MyParameterHandler(ParameterHandler):
            pass

        handler = MyParameterHandler(version=0.3)

        with pytest.raises(
            NotImplementedError, match="no longer create OpenMM forces."
        ):
            handler.create_force()


class TestParameterList:
    """Test capabilities of ParameterList for accessing and manipulating SMIRNOFF parameter definitions."""

    def test_create(self):
        """Test creation of a parameter list."""
        p1 = ParameterType(smirks="[*:1]")
        p2 = ParameterType(smirks="[#1:1]")
        ParameterList([p1, p2])

    @pytest.mark.wip(
        reason="Until ChemicalEnvironment won't be refactored to use the ToolkitRegistry "
        "API, the smirks assignment will fail with RDKit."
    )
    def test_getitem(self):
        """Test ParameterList __getitem__ overloading."""
        p1 = ParameterType(smirks="[*:1]")
        p2 = ParameterType(smirks="[#1:1]")
        parameters = ParameterList([p1, p2])
        assert parameters[0] == p1
        assert parameters[1] == p2
        assert parameters[p1.smirks] == p1
        assert parameters[p2.smirks] == p2

        # Note that this call access __getitem__, not __setitem__.
        parameters["[*:1]"].smirks = "[*X4:1]"
        assert parameters[0].smirks == "[*X4:1]"
        assert p1.smirks == "[*X4:1]"

    def test_index(self):
        """
        Tests the ParameterList.index() function by attempting lookups by SMIRKS and by ParameterType equivalence.

        """
        p1 = ParameterType(smirks="[*:1]")
        p2 = ParameterType(smirks="[#1:1]")
        p3 = ParameterType(smirks="[#7:1]")
        parameters = ParameterList([p1, p2, p3])
        assert parameters.index(p1) == 0
        assert parameters.index(p2) == 1
        assert parameters.index(p3) == 2
        assert parameters.index("[*:1]") == 0
        assert parameters.index("[#1:1]") == 1
        assert parameters.index("[#7:1]") == 2
        with pytest.raises(
            ParameterLookupError, match=r"SMIRKS \[#2:1\] not found in ParameterList"
        ):
            parameters.index("[#2:1]")

        p4 = ParameterType(smirks="[#2:1]")
        with pytest.raises(ValueError, match="is not in list"):
            parameters.index(p4)

    def test_contains(self):
        """Test ParameterList __contains__ overloading."""
        p1 = ParameterType(smirks="[*:1]")
        p2 = ParameterType(smirks="[#1:1]")
        p3 = ParameterType(smirks="[#7:1]")
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
        p1 = ParameterType(smirks="[*:1]")
        p2 = ParameterType(smirks="[#1:1]")
        p3 = ParameterType(smirks="[#7:1]")
        parameters = ParameterList([p1, p2, p3])

        with pytest.raises(IndexError, match="list assignment index out of range"):
            del parameters[4]
        with pytest.raises(
            ParameterLookupError,
            match=r"SMIRKS \[#6:1\] not found in ParameterList",
        ):
            del parameters["[#6:1]"]

        # Test that original list deletion behavior is preserved.
        del parameters[2]
        assert len(parameters) == 2
        assert p1 in parameters
        assert p2 in parameters
        assert p3 not in parameters

        # Test that we can delete elements by their smirks.
        del parameters["[#1:1]"]
        assert len(parameters) == 1
        assert p1 in parameters
        assert p2 not in parameters

    def test_append(self):
        """
        Test ParameterList.append, ensuring that the new parameter was added to the bottom of the list
        and that it is properly recorded as the most recently-added.
        """
        p1 = ParameterType(smirks="[*:1]-[*:2]")
        p2 = ParameterType(smirks="[*:1]=[*:2]")
        param_list = ParameterList()
        param_list.append(p1)
        assert len(param_list) == 1
        assert "[*:1]-[*:2]" in param_list
        param_list.append(p2)
        assert len(param_list) == 2
        assert "[*:1]=[*:2]" in param_list
        assert param_list[-1] == p2

    def test_insert(self):
        """
        Test ParameterList.insert, ensuring that the new parameter was added to the proper spot in
        the list and that it is propertly recorded as the most recently added.
        """
        p1 = ParameterType(smirks="[*:1]-[*:2]")
        p2 = ParameterType(smirks="[*:1]=[*:2]")
        p3 = ParameterType(smirks="[*:1]#[*:2]")
        param_list = ParameterList([p1, p2])
        param_list.insert(1, p3)
        assert param_list[1] == p3

    def test_extend(self):
        """
        Test ParameterList.extend, ensuring that the new parameter was added to the proper spot in
        the list and that it is propertly recorded as the most recently added.
        """
        p1 = ParameterType(smirks="[*:1]-[*:2]")
        p2 = ParameterType(smirks="[*:1]=[*:2]")
        param_list1 = ParameterList()
        param_list2 = ParameterList([p1, p2])

        param_list1.extend(param_list2)
        assert len(param_list1) == 2
        assert "[*:1]-[*:2]" in param_list1
        assert "[*:1]=[*:2]" in param_list1
        assert param_list1[-1] == p2

    def test_to_list(self):
        """Test basic ParameterList.to_list() function, ensuring units are preserved"""
        p1 = BondHandler.BondType(
            smirks="[*:1]-[*:2]",
            length=1.01 * unit.angstrom,
            k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        )
        p2 = BondHandler.BondType(
            smirks="[*:1]=[*:2]",
            length=1.02 * unit.angstrom,
            k=6 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        )
        p3 = BondHandler.BondType(
            smirks="[*:1]#[*:3]",
            length=1.03 * unit.angstrom,
            k=7 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        )
        parameter_list = ParameterList([p1, p2, p3])
        ser_param_list = parameter_list.to_list()
        assert len(ser_param_list) == 3
        assert ser_param_list[0]["length"] == 1.01 * unit.angstrom

    def test_round_trip(self):
        """Test basic ParameterList.to_list() function and constructor"""
        p1 = BondHandler.BondType(
            smirks="[*:1]-[*:2]",
            length=1.01 * unit.angstrom,
            k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        )
        p2 = BondHandler.BondType(
            smirks="[*:1]=[*:2]",
            length=1.02 * unit.angstrom,
            k=6 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        )
        p3 = BondHandler.BondType(
            smirks="[*:1]#[*:3]",
            length=1.03 * unit.angstrom,
            k=7 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        )
        parameter_list = ParameterList([p1, p2, p3])
        param_dict_list = parameter_list.to_list()
        parameter_list_2 = ParameterList()
        for param_dict in param_dict_list:
            new_parameter = BondHandler.BondType(**param_dict)
            parameter_list_2.append(new_parameter)
        assert parameter_list.to_list() == parameter_list_2.to_list()


class TestParameterType:
    def test_find_all_parameter_attrs(self):
        """ParameterType find all ParameterAttributes in the declared order."""

        class MyParameter(ParameterType):
            attr = ParameterAttribute()
            indexed = IndexedParameterAttribute()

        parameter_attributes = MyParameter._get_parameter_attributes()

        # The function should find also the parent's attributes and in the correct order.
        expected_attributes = ["smirks", "id", "parent_id", "attr", "indexed"]
        assert list(parameter_attributes.keys()) == expected_attributes

        # The keys map to the descriptor instances.
        assert type(parameter_attributes["attr"]) is ParameterAttribute
        assert type(parameter_attributes["indexed"]) is IndexedParameterAttribute

    def test_find_all_indexed_parameter_attrs(self):
        """ParameterType find all IndexedParameterAttributes."""

        class MyParameter(ParameterType):
            attr = ParameterAttribute()
            indexed = IndexedParameterAttribute()
            attr2 = ParameterAttribute()
            indexed2 = IndexedParameterAttribute(default=None)

        expected_names = ["indexed", "indexed2"]
        parameter_attributes = MyParameter._get_indexed_parameter_attributes()
        assert list(parameter_attributes.keys()) == expected_names
        assert all(
            isinstance(parameter_attributes[name], IndexedParameterAttribute)
            for name in expected_names
        )

    def test_find_all_required_and_optional_parameter_attrs(self):
        """ParameterType distinguish between required and optional ParameterAttributes."""

        class MyParameter(ParameterType):
            required = ParameterAttribute()
            optional = ParameterAttribute(default=1)
            required_indexed = IndexedParameterAttribute()
            optional_indexed2 = IndexedParameterAttribute(default=None)

        expected_names = ["smirks", "required", "required_indexed"]
        parameter_attributes = MyParameter._get_required_parameter_attributes()
        assert list(parameter_attributes.keys()) == expected_names

        expected_names = ["id", "parent_id", "optional", "optional_indexed2"]
        parameter_attributes = MyParameter._get_optional_parameter_attributes()
        assert list(parameter_attributes.keys()) == expected_names

    def test_required_attribute_on_init(self):
        """ParameterType raises TypeError if a required attribute is not specified on construction."""

        class MyParameter(ParameterType):
            required = ParameterAttribute()
            optional = ParameterAttribute(default=None)

        with pytest.raises(
            SMIRNOFFSpecError, match="require the following missing parameters"
        ):
            MyParameter(smirks="[*:1]", optional=1)

    def test_add_delete_cosmetic_attributes(self):
        """
        Test ParameterType.add_cosmetic_attribute, delete_cosmetic_attribute,
        attribute_is_cosmetic, and to_dict() functions for proper behavior
        """

        class MyParameter(ParameterType):
            required = ParameterAttribute()

        my_par = MyParameter(smirks="[*:1]", required="aaa")
        assert not (my_par.attribute_is_cosmetic("pilot"))

        # Ensure the cosmetic attribute is present by default during output
        my_par.add_cosmetic_attribute("pilot", "alice")
        param_dict = my_par.to_dict()
        assert ("pilot", "alice") in param_dict.items()
        assert my_par.attribute_is_cosmetic("pilot")

        # Ensure the cosmetic attribute isn't present if we request that it be discarded
        param_dict = my_par.to_dict(discard_cosmetic_attributes=True)
        assert "pilot" not in param_dict

        # Manually delete the cosmetic attribute and ensure it doesn't get written out
        my_par.delete_cosmetic_attribute("pilot")
        param_dict = my_par.to_dict()
        assert "pilot" not in param_dict
        assert not (my_par.attribute_is_cosmetic("pilot"))

    def test_indexed_attrs(self):
        """ParameterType handles indexed attributes correctly."""

        class MyParameter(ParameterType):
            a = IndexedParameterAttribute()
            b = IndexedParameterAttribute()

        my_par = MyParameter(smirks="[*:1]", a1=1, a3=3, a2=2, b1=4, b2=5, b3=6)
        assert my_par.a == [1, 2, 3]
        assert my_par.b == [4, 5, 6]

    def test_sequence_init_indexed_attr(self):
        """ParameterType handle indexed attributes initialized with sequences correctly."""

        class MyParameter(ParameterType):
            a = IndexedParameterAttribute()

        my_par = MyParameter(smirks="[*:1]", a=(1, 2))
        assert my_par.a == [1, 2]

    def test_same_length_indexed_attrs(self):
        """ParameterType raises TypeError if indexed attributes of different lengths are given."""

        class MyParameter(ParameterType):
            a = IndexedParameterAttribute()
            b = IndexedParameterAttribute()

        with pytest.raises(
            TypeError, match="indexed attributes have different lengths"
        ):
            MyParameter(smirks="[*:1]", a1=1, a2=2, a3=3, b1=1, b2=2)

    def test_error_single_value_plus_index(self):
        """ParameterType raises an error if an indexed attribute is specified with and without index."""

        class MyParameter(ParameterType):
            a = IndexedParameterAttribute()

        with pytest.raises(
            TypeError, match="'a' has been specified with and without index"
        ):
            MyParameter(smirks="[*:1]", a=[1], a1=2)

    def test_find_all_defined_parameter_attrs(self):
        """ParameterType._get_defined_attributes() discards None default-value attributes."""

        class MyParameter(ParameterType):
            required1 = ParameterAttribute()
            optional1 = ParameterAttribute(default=None)
            optional2 = IndexedParameterAttribute(default=None)
            optional3 = ParameterAttribute(default=5)
            required2 = IndexedParameterAttribute()
            optional4 = ParameterAttribute(default=2)

        my_par = MyParameter(smirks="[*:1]", required1=0, optional1=10, required2=[0])

        # _get_defined_parameter_attributes discards only the attribute
        # that are set to None as a default value.
        expected_names = [
            "smirks",
            "required1",
            "required2",
            "optional1",
            "optional3",
            "optional4",
        ]
        parameter_attributes = my_par._get_defined_parameter_attributes()
        assert list(parameter_attributes.keys()) == expected_names

    def test_base_parametertype_to_dict(self):
        """
        Test ParameterType to_dict.
        """
        p1 = ParameterType(smirks="[*:1]")
        param_dict = p1.to_dict()
        assert param_dict["smirks"] == "[*:1]"
        assert len(param_dict.keys()) == 1

    def test_repr(self):
        class NamedType(ParameterType):
            pass

        assert NamedType(smirks="[*:1]").__repr__().startswith("<NamedType")


class TestBondType:
    """Tests for the BondType class."""

    def test_bondtype_to_dict(self):
        """
        Test BondType to_dict.
        """
        p1 = BondHandler.BondType(
            smirks="[*:1]-[*:2]",
            length=1.02 * unit.angstrom,
            k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        )
        param_dict = p1.to_dict()

        assert param_dict == {
            "smirks": "[*:1]-[*:2]",
            "length": unit.Quantity(1.02, unit.angstrom),
            "k": unit.Quantity(5, unit.kilocalorie_per_mole / unit.angstrom**2),
        }

    def test_bondtype_partial_bondorders(self):
        """
        Test the parsing of a BondType with k_bondorder1/2/3 definitions
        """
        length = 1.4 * unit.angstrom
        k1 = 101 * unit.kilocalorie / unit.mole / unit.angstrom**2
        k2 = 202 * unit.kilocalorie / unit.mole / unit.angstrom**2
        k3 = 303 * unit.kilocalorie / unit.mole / unit.angstrom**2

        param = BondHandler.BondType(
            smirks="[*:1]-[*:2]",
            length=length,
            k_bondorder1=k1,
            k_bondorder2=k2,
            k_bondorder3=k3,
        )

        assert param.k_bondorder == {1: k1, 2: k2, 3: k3}

    def test_bondtype_bad_params(self):
        """
        Test the over/underspecification of k/k_bondorderN are caught
        """
        length = 1.4 * unit.angstrom
        length1 = 1.5 * unit.angstrom
        length2 = 1.3 * unit.angstrom
        k = 50 * unit.kilocalorie / unit.mole / unit.angstrom**2
        k1 = 101 * unit.kilocalorie / unit.mole / unit.angstrom**2
        k2 = 202 * unit.kilocalorie / unit.mole / unit.angstrom**2

        with pytest.raises(SMIRNOFFSpecError, match="Either k or k_bondorder"):
            BondHandler.BondType(
                smirks="[*:1]-[*:2]",
                length=length,
            )

        with pytest.raises(SMIRNOFFSpecError, match="BOTH k and k_bondorder"):
            BondHandler.BondType(
                smirks="[*:1]-[*:2]",
                length=length,
                k=k,
                k_bondorder1=k1,
                k_bondorder2=k2,
            )

        with pytest.raises(
            SMIRNOFFSpecError, match="Either length or length_bondorder"
        ):
            BondHandler.BondType(
                smirks="[*:1]-[*:2]",
                k=k,
            )

        with pytest.raises(SMIRNOFFSpecError, match="BOTH length and length_bondorder"):
            BondHandler.BondType(
                smirks="[*:1]-[*:2]",
                length=length,
                k=k,
                length_bondorder1=length1,
                length_bondorder2=length2,
            )

    def test_bondtype_to_dict_custom_output_units(self):
        """
        Test BondType to_dict with custom output units.
        """
        p1 = BondHandler.BondType(
            smirks="[*:1]-[*:2]",
            length=1.02 * unit.angstrom,
            k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        )

        param_dict = p1.to_dict()

        assert abs(param_dict["length"].m_as(unit.nanometer) - 0.102) < 1e-10

    def test_read_write_optional_parameter_attribute(self):
        """
        Test ParameterTypes' ability to store and write out optional attributes passed to __init__()
        """
        p1 = BondHandler.BondType(
            smirks="[*:1]-[*:2]",
            length=1.02 * unit.angstrom,
            k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
            id="b1",
        )
        param_dict = p1.to_dict()
        assert ("id", "b1") in param_dict.items()

    def test_read_write_cosmetic_parameter_attribute(self):
        """
        Test ParameterTypes' ability to store and write out cosmetic attributes passed to __init__()
        """
        p1 = BondHandler.BondType(
            smirks="[*:1]-[*:2]",
            length=1.02 * unit.angstrom,
            k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
            pilot="alice",
            allow_cosmetic_attributes=True,
        )
        param_dict = p1.to_dict(discard_cosmetic_attributes=False)
        assert ("pilot", "alice") in param_dict.items()

    def test_read_but_dont_write_cosmetic_parameter_attribute(self):
        """
        Test ParameterTypes' ability to ignore cosmetic attributes passed to __init__() if instructed
        """
        p1 = BondHandler.BondType(
            smirks="[*:1]-[*:2]",
            length=1.02 * unit.angstrom,
            k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
            pilot="alice",
            allow_cosmetic_attributes=True,
        )
        param_dict = p1.to_dict(discard_cosmetic_attributes=True)
        assert ("pilot", "alice") not in param_dict.items()

    def test_error_cosmetic_parameter_attribute(self):
        """
        Test that ParameterTypes raise an error on receiving unexpected attributes passed to __init__()
        """
        with pytest.raises(SMIRNOFFSpecError, match="Unexpected kwarg (pilot: alice)*"):
            BondHandler.BondType(
                smirks="[*:1]-[*:2]",
                length=1.02 * unit.angstrom,
                k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
                pilot="alice",
                allow_cosmetic_attributes=False,
            )

    def test_add_delete_cosmetic_attrib(self):
        """
        Test adding and deleting cosmetic attributes for already-initialized ParameterType objects
        """
        p1 = BondHandler.BondType(
            smirks="[*:1]-[*:2]",
            length=1.02 * unit.angstrom,
            k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        )
        # Ensure the cosmetic attribute is present by default during output
        p1.add_cosmetic_attribute("pilot", "alice")
        param_dict = p1.to_dict()
        assert ("pilot", "alice") in param_dict.items()

        # Ensure the cosmetic attribute isn't present if we request that it be discarded
        param_dict = p1.to_dict(discard_cosmetic_attributes=True)
        assert ("pilot", "alice") not in param_dict.items()

        # Manually delete the cosmetic attribute and ensure it doesn't get written out
        p1.delete_cosmetic_attribute("pilot")
        param_dict = p1.to_dict()
        assert ("pilot", "alice") not in param_dict.items()


class TestBondHandler:
    @pytest.mark.parametrize(
        ("fractional_bond_order", "k_interpolated", "length_interpolated"),
        [
            (1.0, 101, 1.4),
            (1.5, 112.0, 1.35),
            (1.99, 122.78, 1.301),
            (2.1, 125.2, 1.29),
        ],
    )
    def test_linear_interpolate(
        self, fractional_bond_order, k_interpolated, length_interpolated
    ):
        """Test that linear interpolation works as expected"""
        k_bondorder = {
            1: 101 * unit.kilocalorie / unit.mole / unit.angstrom**2,
            2: 123 * unit.kilocalorie / unit.mole / unit.angstrom**2,
        }

        length_bondorder = {
            1: 1.4 * unit.angstrom,
            2: 1.3 * unit.angstrom,
        }

        k = _linear_inter_or_extrapolate(k_bondorder, fractional_bond_order)
        length = _linear_inter_or_extrapolate(length_bondorder, fractional_bond_order)
        assert_almost_equal(k.m, k_interpolated, 1)
        assert_almost_equal(length.m, length_interpolated, 2)

    def test_different_defaults_03_04(self):
        """Ensure that the 0.3 and 0.4 versions' defaults are correctly set"""
        bh = BondHandler(version=0.3)
        assert bh.fractional_bondorder_method == "none"
        assert bh.potential == "harmonic"
        bh2 = BondHandler(version=0.4)
        assert bh2.fractional_bondorder_method == "AM1-Wiberg"
        assert bh2.potential == "(k/2)*(r-length)^2"
        bh3 = BondHandler(version=0.3, fractional_bondorder_method="AM1-Wiberg")
        assert bh3.fractional_bondorder_method == "AM1-Wiberg"
        assert bh3.fractional_bondorder_interpolation == "linear"
        assert bh3.potential == "harmonic"

    def test_harmonic_potentials_are_compatible(self):
        """
        Ensure that handlers with `potential="harmonic"` evaluate as compatible with handlers with
        potential="(k/2)*(r-length)^2"
        """
        bh1 = BondHandler(skip_version_check=True)
        bh2 = BondHandler(skip_version_check=True)
        bh1.potential = "harmonic"
        bh2.potential = "(k/2)*(r-length)^2"
        # This comparison should pass, since the potentials defined above are compatible
        bh1.check_handler_compatibility(bh2)

    def test_am1_wiberg_combine_error(self):
        """Reproduce issue #719 and a suggested fix."""
        v3 = BondHandler(version=0.3)
        v4 = BondHandler(version=0.4)

        v3.add_parameter(
            parameter=BondHandler.BondType(
                smirks="[#6:1]-[*:2]",
                length=1.0 * unit.angstrom,
                k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
            )
        )

        v4.add_parameter(
            parameter=BondHandler.BondType(
                smirks="[#8:1]-[*:2]",
                length=2.0 * unit.angstrom,
                k=5 * unit.kilocalorie / unit.mole / unit.angstrom**2,
            )
        )

        ff_v3 = ForceField()
        ff_v4 = ForceField()

        ff_v3.register_parameter_handler(v3)
        ff_v4.register_parameter_handler(v4)

        with pytest.raises(
            IncompatibleParameterError,
            match="This likely results from mixing bond handlers",
        ):
            ff_v3._load_smirnoff_data(ff_v4._to_smirnoff_data())

        # Simulate a user following the recommendation to switch a
        # version 0.3 BondHandler to use 'AM1-Wiberg'
        ff_v3["Bonds"].fractional_bondorder_method = "AM1-Wiberg"

        ff_v3._load_smirnoff_data(ff_v4._to_smirnoff_data())


class TestProperTorsionType:
    """Tests for the ProperTorsionType class."""

    def test_single_term_proper_torsion(self):
        """
        Test creation and serialization of a single-term proper torsion
        """
        p1 = ProperTorsionHandler.ProperTorsionType(
            smirks="[*:1]-[*:2]-[*:3]-[*:4]",
            phase1=30 * unit.degree,
            periodicity1=2,
            k1=5 * unit.kilocalorie / unit.mole,
        )
        param_dict = p1.to_dict()
        assert ("k1", 5 * unit.kilocalorie / unit.mole) in param_dict.items()
        assert ("phase1", 30 * unit.degree) in param_dict.items()
        assert ("periodicity1", 2) in param_dict.items()
        assert "idivf" not in param_dict

    def test_single_term_proper_torsion_w_idivf(self):
        """
        Test creation and serialization of a single-term proper torsion
        """
        p1 = ProperTorsionHandler.ProperTorsionType(
            smirks="[*:1]-[*:2]-[*:3]-[*:4]",
            phase1=30 * unit.degree,
            periodicity1=2,
            k1=5 * unit.kilocalorie / unit.mole,
            idivf1=4,
        )

        param_dict = p1.to_dict()
        assert ("k1", 5 * unit.kilocalorie / unit.mole) in param_dict.items()
        assert ("phase1", 30 * unit.degree) in param_dict.items()
        assert ("periodicity1", 2) in param_dict.items()
        assert ("idivf1", 4) in param_dict.items()

    def test_multi_term_proper_torsion(self):
        """
        Test creation and serialization of a multi-term proper torsion
        """
        p1 = ProperTorsionHandler.ProperTorsionType(
            smirks="[*:1]-[*:2]-[*:3]-[*:4]",
            phase1=30 * unit.degree,
            periodicity1=2,
            k1=5 * unit.kilocalorie / unit.mole,
            phase2=31 * unit.degree,
            periodicity2=3,
            k2=6 * unit.kilocalorie / unit.mole,
        )
        param_dict = p1.to_dict()
        assert param_dict["k1"] == 5 * unit.kilocalorie / unit.mole
        assert param_dict["phase1"] == 30 * unit.degree
        assert param_dict["periodicity1"] == 2
        assert param_dict["k2"] == 6 * unit.kilocalorie / unit.mole
        assert param_dict["phase2"] == 31 * unit.degree
        assert param_dict["periodicity2"] == 3

    def test_multi_term_proper_torsion_skip_index(self):
        """
        Test creation and serialization of a multi-term proper torsion where
        the indices are not consecutive and a SMIRNOFFSpecError is raised
        """
        with pytest.raises(
            SMIRNOFFSpecError, match=r"Unexpected kwarg \(phase3: 31 deg\)*."
        ):
            ProperTorsionHandler.ProperTorsionType(
                smirks="[*:1]-[*:2]-[*:3]-[*:4]",
                phase1=30 * unit.degree,
                periodicity1=2,
                k1=5 * unit.kilocalorie / unit.mole,
                phase3=31 * unit.degree,
                periodicity3=3,
                k3=6 * unit.kilocalorie / unit.mole,
            )

    def test_multi_term_proper_torsion_bad_units(self):
        """
        Test creation and serialization of a multi-term proper torsion where
        one of the terms has incorrect units
        """
        with pytest.raises(IncompatibleUnitError, match="should have units of"):
            ProperTorsionHandler.ProperTorsionType(
                smirks="[*:1]-[*:2]-[*:3]-[*:4]",
                phase1=30 * unit.degree,
                periodicity1=2,
                k1=5 * unit.kilocalorie / unit.mole,
                phase2=31 * unit.angstrom,  # This should be caught
                periodicity2=3,
                k2=6 * unit.kilocalorie / unit.mole,
            )

    def test_single_term_proper_torsion_bo(self):
        """
        Test creation and serialization of a single-term proper torsion with bond order interpolation.
        """
        p1 = ProperTorsionHandler.ProperTorsionType(
            smirks="[*:1]-[*:2]~[*:3]-[*:4]",
            phase1=30 * unit.degree,
            periodicity1=2,
            k1_bondorder1=1 * unit.kilocalorie / unit.mole,
            k1_bondorder2=1.8 * unit.kilocalorie / unit.mole,
        )
        param_dict = p1.to_dict()
        assert ("k1_bondorder1", 1 * unit.kilocalorie / unit.mole) in param_dict.items()
        assert (
            "k1_bondorder2",
            1.8 * unit.kilocalorie / unit.mole,
        ) in param_dict.items()
        assert ("phase1", 30 * unit.degree) in param_dict.items()
        assert ("periodicity1", 2) in param_dict.items()
        assert "idivf" not in param_dict

        assert len(p1.k_bondorder) == 1
        assert len(p1.k_bondorder[0]) == 2
        assert {1, 2} == set(p1.k_bondorder[0].keys())

    def test_single_term_proper_torsion_bo_w_idivf(self):
        """
        Test creation and serialization of a single-term proper torsion with bond order interpolation.

        With `idivf1` specified.

        """
        p1 = ProperTorsionHandler.ProperTorsionType(
            smirks="[*:1]-[*:2]-[*:3]-[*:4]",
            phase1=30 * unit.degree,
            periodicity1=2,
            k1_bondorder1=1 * unit.kilocalorie / unit.mole,
            k1_bondorder2=1.8 * unit.kilocalorie / unit.mole,
            idivf1=4,
        )

        param_dict = p1.to_dict()
        assert ("k1_bondorder1", 1 * unit.kilocalorie / unit.mole) in param_dict.items()
        assert (
            "k1_bondorder2",
            1.8 * unit.kilocalorie / unit.mole,
        ) in param_dict.items()
        assert ("phase1", 30 * unit.degree) in param_dict.items()
        assert ("periodicity1", 2) in param_dict.items()
        assert ("idivf1", 4) in param_dict.items()

    def test_multi_term_proper_torsion_bo(self):
        """
        Test creation and serialization of a multi-term proper torsion with bond order interpolation.
        """
        p1 = ProperTorsionHandler.ProperTorsionType(
            smirks="[*:1]-[*:2]-[*:3]-[*:4]",
            phase1=30 * unit.degree,
            periodicity1=2,
            k1_bondorder1=1 * unit.kilocalorie / unit.mole,
            k1_bondorder2=1.8 * unit.kilocalorie / unit.mole,
            phase2=31 * unit.degree,
            periodicity2=3,
            k2_bondorder1=1.2 * unit.kilocalorie / unit.mole,
            k2_bondorder2=1.9 * unit.kilocalorie / unit.mole,
        )
        param_dict = p1.to_dict()
        assert param_dict["k1_bondorder1"] == 1 * unit.kilocalorie / unit.mole
        assert param_dict["k1_bondorder2"] == 1.8 * unit.kilocalorie / unit.mole
        assert param_dict["phase1"] == 30 * unit.degree
        assert param_dict["periodicity1"] == 2
        assert param_dict["k2_bondorder1"] == 1.2 * unit.kilocalorie / unit.mole
        assert param_dict["k2_bondorder2"] == 1.9 * unit.kilocalorie / unit.mole
        assert param_dict["phase2"] == 31 * unit.degree
        assert param_dict["periodicity2"] == 3

    def test_multi_term_proper_torsion_bo_getters_setters(self):
        """
        Test getters and setters of a multi-term proper torsion with bond order interpolation.
        """
        p1 = ProperTorsionHandler.ProperTorsionType(
            smirks="[*:1]-[*:2]-[*:3]-[*:4]",
            phase1=30 * unit.degree,
            periodicity1=2,
            k1_bondorder1=1 * unit.kilocalorie / unit.mole,
            k1_bondorder2=1.8 * unit.kilocalorie / unit.mole,
            phase2=31 * unit.degree,
            periodicity2=3,
            k2_bondorder1=1.2 * unit.kilocalorie / unit.mole,
            k2_bondorder2=1.9 * unit.kilocalorie / unit.mole,
        )

        assert p1.k1_bondorder1 == 1.0 * unit.kilocalorie / unit.mole
        p1.k1_bondorder1 = 2.0 * unit.kilocalorie / unit.mole
        assert p1.k1_bondorder1 == 2.0 * unit.kilocalorie / unit.mole

        assert p1.k2_bondorder2 == 1.9 * unit.kilocalorie / unit.mole
        p1.k2_bondorder2 = 2.9 * unit.kilocalorie / unit.mole
        assert p1.k2_bondorder2 == 2.9 * unit.kilocalorie / unit.mole

    def test_multi_term_proper_torsion_bo_skip_index(self):
        """
        Test creation and serialization of a multi-term proper torsion where
        the indices are not consecutive and a SMIRNOFFSpecError is raised
        AND we are doing bond order interpolation
        """
        with pytest.raises(
            SMIRNOFFSpecError, match=r"Unexpected kwarg \(k3_bondorder1*."
        ):
            ProperTorsionHandler.ProperTorsionType(
                smirks="[*:1]-[*:2]-[*:3]-[*:4]",
                phase1=30 * unit.degree,
                periodicity1=2,
                k1_bondorder1=1 * unit.kilocalorie / unit.mole,
                k1_bondorder2=1.8 * unit.kilocalorie / unit.mole,
                phase3=31 * unit.degree,
                periodicity3=3,
                k3_bondorder1=1.2 * unit.kilocalorie / unit.mole,
                k3_bondorder2=1.9 * unit.kilocalorie / unit.mole,
            )

    def test_single_term_single_bo_exception(self):
        """Test behavior where a single bond order term is specified for a single k"""
        # raises no error, as checks are handled at parameterization
        # we may add a `validate` method later that is called manually by user when they want it
        ProperTorsionHandler.ProperTorsionType(
            smirks="[*:1]-[*:2]~[*:3]-[*:4]",
            phase1=30 * unit.degree,
            periodicity1=2,
            k1_bondorder1=1 * unit.kilocalorie / unit.mole,
        )

    def test_multi_term_single_bo_exception(self):
        """Test behavior where a single bond order term is specified for each of multiple k"""
        # TODO : currently raises no error, as checks are handled at parameterization
        # is this a spec thing that we should be checking?
        # if so, it will be painful to implement
        ProperTorsionHandler.ProperTorsionType(
            smirks="[*:1]-[*:2]-[*:3]-[*:4]",
            phase1=30 * unit.degree,
            periodicity1=2,
            k1_bondorder1=1 * unit.kilocalorie / unit.mole,
            phase2=31 * unit.degree,
            periodicity2=3,
            k2_bondorder1=1.2 * unit.kilocalorie / unit.mole,
        )


class TestProperTorsionHandler:
    def test_torsion_handler_charmm_potential(self):
        """
        Test creation of TorsionHandlers with the deprecated 0.2 potential value "charmm" instead of the current
        supported potential value "fourier".
        """
        import re

        # Test creating ProperTorsionHandlers
        err_msg = re.escape(
            "Attempted to set ProperTorsionHandler.potential to charmm. Currently, "
            "only the following values are supported: ['k*(1+cos(periodicity*theta-phase))']."
        )
        with pytest.raises(SMIRNOFFSpecError, match=err_msg):
            ProperTorsionHandler(potential="charmm", skip_version_check=True)
        ProperTorsionHandler(
            potential="k*(1+cos(periodicity*theta-phase))", skip_version_check=True
        )

        # Same test, but with ImproperTorsionHandler
        err_msg = re.escape(
            "Attempted to set ImproperTorsionHandler.potential to charmm. Currently, "
            "only the following values are supported: ['k*(1+cos(periodicity*theta-phase))']."
        )
        with pytest.raises(SMIRNOFFSpecError, match=err_msg):
            ImproperTorsionHandler(potential="charmm", skip_version_check=True)
        ImproperTorsionHandler(
            potential="k*(1+cos(periodicity*theta-phase))", skip_version_check=True
        )


class TestvdWHandler:
    def test_add_param_str(self):
        """
        Ensure that string input is supported, given the added complication that the
        sigma/rmin_half setters silently set each other's value.
        See https://github.com/openforcefield/openff-toolkit/issues/788
        """
        vdw_handler = vdWHandler(version=0.3)
        param1 = {
            "epsilon": "0.5 * kilocalorie/mole",
            "rmin_half": "1.2 * angstrom",
            "smirks": "[*:1]",
            "id": "n99",
        }
        param2 = {
            "epsilon": "0.1 * kilocalorie/mole",
            "sigma": "0.8 * angstrom",
            "smirks": "[#1:1]",
            "id": "n00",
        }
        vdw_handler.add_parameter(param1)
        vdw_handler.add_parameter(param2)

        assert vdw_handler.get_parameter({"smirks": "[*:1]"})[0].id == "n99"
        assert vdw_handler.get_parameter({"smirks": "[#1:1]"})[0].id == "n00"

    def test_set_invalid_scale_factor(self):
        handler = vdWHandler(version=0.3)

        with pytest.raises(SMIRNOFFSpecError, match="unable to handle scale12"):
            handler.scale12 = 0.1

        with pytest.raises(SMIRNOFFSpecError, match="unable to handle scale13"):
            handler.scale13 = 0.1

        with pytest.raises(SMIRNOFFSpecError, match="unable to handle scale15"):
            handler.scale15 = 0.1


class TestvdWType:
    """
    Test the behavior of vdWType
    """

    def test_sigma_rmin_half(self):
        """Test the setter/getter behavior or sigma and rmin_half"""
        from openff.toolkit.typing.engines.smirnoff.parameters import vdWHandler

        data = {
            "smirks": "[*:1]",
            "rmin_half": 0.5 * unit.angstrom,
            "epsilon": 0.5 * unit.kilocalorie / unit.mole,
        }
        param = vdWHandler.vdWType(**data)

        assert param.sigma is not None
        assert param.rmin_half is not None
        assert numpy.isclose(
            param.sigma.m_as(unit.angstrom),
            (2.0 * param.rmin_half / 2 ** (1 / 6)).m_as(unit.angstrom),
        )
        assert "sigma" not in param.to_dict()
        assert "rmin_half" in param.to_dict()

        param.sigma = param.sigma

        assert numpy.isclose(param.rmin_half.m_as(unit.angstrom), 0.5)
        assert "sigma" in param.to_dict()
        assert "rmin_half" not in param.to_dict()

        param.rmin_half = param.rmin_half

        assert numpy.isclose(
            param.sigma.m_as(unit.angstrom),
            (2.0 * param.rmin_half / 2 ** (1 / 6)).m_as(unit.angstrom),
        )
        assert "sigma" not in param.to_dict()
        assert "rmin_half" in param.to_dict()


class TestElectrostaticsHandler:
    def test_solvent_dielectric(self):
        with pytest.raises(
            SMIRNOFFSpecUnimplementedError,
            match="make use of `solvent_d",
        ):
            ElectrostaticsHandler(version=0.3, method="PME", solvent_dielectric=4)

        handler = ElectrostaticsHandler(version=0.4)
        assert handler.solvent_dielectric is None

        handler.solvent_dielectric = None

        with pytest.raises(
            SMIRNOFFSpecUnimplementedError,
            match="make use of `solvent_d",
        ):
            handler.solvent_dielectric = 78.3

    def test_unknown_periodic_potential(self):
        handler = ElectrostaticsHandler(version=0.4)

        with pytest.raises(
            NotImplementedError,
            match="unexpected periodic potential",
        ):
            handler.periodic_potential = "PPPM"

    def test_set_invalid_scale_factor(self):
        handler = ElectrostaticsHandler(version=0.4)

        with pytest.raises(SMIRNOFFSpecError, match="unable to handle scale12"):
            handler.scale12 = 0.1

        with pytest.raises(SMIRNOFFSpecError, match="unable to handle scale13"):
            handler.scale13 = 0.1

        with pytest.raises(SMIRNOFFSpecError, match="unable to handle scale15"):
            handler.scale15 = 0.1


class TestElectrostaticsHandlerUpconversion:
    """
    Test the implementation of OFF-EP-0005:

    https://openforcefield.github.io/standards/enhancement-proposals/off-ep-0005/
    """

    default_reaction_field_expression = (
        "charge1*charge2/(4*pi*epsilon0)*(1/r + k_rf*r^2 - c_rf);"
        "k_rf=(cutoff^(-3))*(solvent_dielectric-1)/(2*solvent_dielectric+1);"
        "c_rf=cutoff^(-1)*(3*solvent_dielectric)/(2*solvent_dielectric+1)"
    )

    @pytest.mark.parametrize(
        ("old_method", "new_method"),
        [
            ("PME", "Ewald3D-ConductingBoundary"),
            ("Coulomb", "Coulomb"),
            ("reaction-field", default_reaction_field_expression),
        ],
    )
    def test_upconversion(self, old_method, new_method):
        handler = ElectrostaticsHandler(version=0.3, method=old_method)
        assert handler.version == Version("0.4")

        # Only `periodic_potential` is a function of the values in a version 0.3 handler ...
        assert handler.periodic_potential == new_method

        # ... for everything else, it's the same
        assert handler.nonperiodic_potential == "Coulomb"
        assert handler.exception_potential == "Coulomb"

    def test_invalid_0_4_kwargs(self):
        with pytest.raises(SMIRNOFFSpecError, match="removed in version 0.4 of the E"):
            ElectrostaticsHandler(version=0.4, method="PME")


class TestVirtualSiteHandler:
    """
    Test the creation of a VirtualSiteHandler and the implemented VirtualSiteTypes
    """

    @pytest.mark.parametrize(
        "parameter, expected_index",
        [
            (VirtualSiteMocking.bond_charge_parameter("[*:1][*:2]"), 0),
            (VirtualSiteMocking.monovalent_parameter("[*:1][*:2][*:3]"), 0),
            (VirtualSiteMocking.divalent_parameter("[*:2][*:1][*:3]", "once"), 0),
            (VirtualSiteMocking.trivalent_parameter("[*:1][*:2][*:3][*:4]"), 0),
        ],
    )
    def test_parent_index(self, parameter, expected_index):
        assert parameter.parent_index == expected_index
        assert (
            VirtualSiteHandler.VirtualSiteType.type_to_parent_index(parameter.type)
            == expected_index
        )

    @pytest.mark.parametrize(
        "kwargs, expected_raises",
        [
            (
                VirtualSiteMocking.bond_charge_parameter("[*:1][*:2]").to_dict(),
                does_not_raise(),
            ),
            (
                VirtualSiteMocking.monovalent_parameter("[*:1][*:2][*:3]").to_dict(),
                does_not_raise(),
            ),
            (
                VirtualSiteMocking.divalent_parameter(
                    "[*:1][*:2][*:3]", match="once", angle=0.0 * unit.degrees
                ).to_dict(),
                does_not_raise(),
            ),
            (
                VirtualSiteMocking.divalent_parameter(
                    "[*:1][*:2][*:3]",
                    match="all_permutations",
                    angle=2.0 * unit.degrees,
                ).to_dict(),
                does_not_raise(),
            ),
            (
                VirtualSiteMocking.trivalent_parameter(
                    "[*:1][*:2][*:3][*:4]"
                ).to_dict(),
                does_not_raise(),
            ),
            # Validate `type`
            (
                {},
                pytest.raises(SMIRNOFFSpecError, match="the `type` keyword is missing"),
            ),
            (
                {"type": "InvalidType"},
                pytest.raises(
                    SMIRNOFFSpecError,
                    match="'InvalidType' is not a supported virtual site type",
                ),
            ),
            # Validate `match`
            (
                {"type": "BondCharge"},
                pytest.raises(
                    SMIRNOFFSpecError, match="the `match` keyword is missing"
                ),
            ),
            (
                {"type": "BondCharge", "match": "once"},
                pytest.raises(
                    SMIRNOFFSpecError,
                    match="match='once' not supported with type='BondCharge'",
                ),
            ),
            (
                {"type": "BondCharge", "match": "once"},
                pytest.raises(
                    SMIRNOFFSpecError,
                    match="match='once' not supported with type='BondCharge'",
                ),
            ),
            (
                {"type": "MonovalentLonePair", "match": "once"},
                pytest.raises(
                    SMIRNOFFSpecError,
                    match="match='once' not supported with type='MonovalentLonePair'",
                ),
            ),
            (
                {
                    "type": "DivalentLonePair",
                    "match": "once",
                    "outOfPlaneAngle": 2.0 * unit.degrees,
                },
                pytest.raises(
                    SMIRNOFFSpecError,
                    match="match='once' not supported with "
                    "type='DivalentLonePair' and is_in_plane=False",
                ),
            ),
            (
                {"type": "TrivalentLonePair", "match": "all_permutations"},
                pytest.raises(
                    SMIRNOFFSpecError,
                    match="match='all_permutations' not supported with "
                    "type='TrivalentLonePair'",
                ),
            ),
        ],
    )
    def test_add_default_init_kwargs_validation(self, kwargs, expected_raises):
        with expected_raises:
            VirtualSiteHandler.VirtualSiteType._add_default_init_kwargs(kwargs)

    @pytest.mark.parametrize(
        "kwargs, expected_kwargs",
        [
            (
                {"type": "BondCharge", "match": "all_permutations"},
                {
                    "type": "BondCharge",
                    "match": "all_permutations",
                    "outOfPlaneAngle": None,
                    "inPlaneAngle": None,
                    "sigma": 0.0 * unit.angstrom,
                    "epsilon": 0.0 * unit.kilocalorie_per_mole,
                },
            ),
            (
                {
                    "type": "BondCharge",
                    "match": "all_permutations",
                    "rmin_half": 1.0 * unit.angstrom,
                },
                {
                    "type": "BondCharge",
                    "match": "all_permutations",
                    "outOfPlaneAngle": None,
                    "inPlaneAngle": None,
                    "rmin_half": 1.0 * unit.angstrom,
                    "epsilon": 0.0 * unit.kilocalorie_per_mole,
                },
            ),
            (
                {"type": "MonovalentLonePair", "match": "all_permutations"},
                {
                    "type": "MonovalentLonePair",
                    "match": "all_permutations",
                    "sigma": 0.0 * unit.angstrom,
                    "epsilon": 0.0 * unit.kilocalorie_per_mole,
                },
            ),
        ],
    )
    def test_add_default_init_kwargs_values(self, kwargs, expected_kwargs):
        assert kwargs != expected_kwargs
        VirtualSiteHandler.VirtualSiteType._add_default_init_kwargs(kwargs)
        assert kwargs == expected_kwargs

    @pytest.mark.parametrize(
        "parameter, in_plane_angle, expected_raises",
        [
            (
                VirtualSiteMocking.bond_charge_parameter("[*:1][*:2]"),
                0.0 * unit.degrees,
                pytest.raises(
                    SMIRNOFFSpecError,
                    match="'BondCharge' sites do not support `inPlaneAngle`",
                ),
            ),
            (
                VirtualSiteMocking.bond_charge_parameter("[*:1][*:2]"),
                None,
                does_not_raise(),
            ),
            (
                VirtualSiteMocking.monovalent_parameter("[*:1][*:2][*:3]"),
                1.0 * unit.angstrom,
                pytest.raises(IncompatibleUnitError),
            ),
            (
                VirtualSiteMocking.monovalent_parameter("[*:1][*:2][*:3]"),
                130.0 * unit.degrees,
                does_not_raise(),
            ),
        ],
    )
    def test_in_plane_angle_converter(self, parameter, in_plane_angle, expected_raises):
        parameter_dict = parameter.to_dict()
        parameter_dict["inPlaneAngle"] = in_plane_angle

        with expected_raises:
            new_parameter = VirtualSiteHandler.VirtualSiteType(**parameter_dict)
            assert new_parameter.inPlaneAngle == in_plane_angle

    @pytest.mark.parametrize(
        "parameter, out_of_plane_angle, expected_raises",
        [
            (
                VirtualSiteMocking.bond_charge_parameter("[*:1][*:2]"),
                0.0 * unit.degrees,
                pytest.raises(
                    SMIRNOFFSpecError,
                    match="'BondCharge' sites do not support `outOfPlaneAngle`",
                ),
            ),
            (
                VirtualSiteMocking.bond_charge_parameter("[*:1][*:2]"),
                None,
                does_not_raise(),
            ),
            (
                VirtualSiteMocking.monovalent_parameter("[*:1][*:2][*:3]"),
                130.0 * unit.degrees,
                does_not_raise(),
            ),
        ],
    )
    def test_out_of_plane_angle_converter(
        self, parameter, out_of_plane_angle, expected_raises
    ):
        parameter_dict = parameter.to_dict()
        parameter_dict["outOfPlaneAngle"] = out_of_plane_angle

        with expected_raises:
            new_parameter = VirtualSiteHandler.VirtualSiteType(**parameter_dict)
            assert new_parameter.outOfPlaneAngle == out_of_plane_angle

    def test_serialize_roundtrip(self):
        force_field = ForceField()

        handler = force_field.get_parameter_handler("VirtualSites")
        handler.add_parameter(
            parameter=VirtualSiteMocking.bond_charge_parameter("[*:1][*:2]")
        )
        handler.add_parameter(
            parameter=VirtualSiteMocking.monovalent_parameter("[*:1][*:2][*:3]")
        )
        handler.add_parameter(
            parameter=VirtualSiteMocking.divalent_parameter("[*:2][*:1][*:3]", "once")
        )
        handler.add_parameter(
            parameter=VirtualSiteMocking.trivalent_parameter("[*:1][*:2][*:3][*:4]")
        )
        offxml_string = force_field.to_string()
        roundtripped_force_field = ForceField(offxml_string)
        assert offxml_string == roundtripped_force_field.to_string()

    @pytest.mark.parametrize(
        "smiles, matched_indices, parameter, expected_raises",
        [
            (
                "[Cl:1][H:2]",
                (1, 2),
                VirtualSiteMocking.bond_charge_parameter("[Cl:1][H:2]"),
                does_not_raise(),
            ),
            (
                "[N:1]([H:2])([H:3])[H:4]",
                (1, 2, 3),
                VirtualSiteMocking.monovalent_parameter("[*:2][N:1][*:3]"),
                pytest.raises(NotImplementedError, match="please describe what it is"),
            ),
        ],
    )
    def test_validate_found_match(
        self, smiles, matched_indices, parameter, expected_raises
    ):
        molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        topology: Topology = molecule.to_topology()

        topology_atoms = {
            i: topology_atom for i, topology_atom in enumerate(topology.topology_atoms)
        }

        with expected_raises:
            VirtualSiteHandler._validate_found_match(
                topology_atoms, matched_indices, parameter
            )

    @pytest.mark.parametrize(
        "handler_a, handler_b, expected_raises",
        [
            # Currently, no other test cases as only one `exclusion_policy` is supported
            (
                VirtualSiteHandler(version="0.3"),
                VirtualSiteHandler(version="0.3"),
                does_not_raise(),
            )
        ],
    )
    def test_check_handler_compatability(self, handler_a, handler_b, expected_raises):
        with expected_raises:
            handler_a.check_handler_compatibility(handler_b)

    @pytest.mark.parametrize(
        "parameters, smiles, expected_matches",
        [
            # Check that a basic BondCharge virtual site can be applied
            (
                [VirtualSiteMocking.bond_charge_parameter("[Cl:1]-[C:2]")],
                "[Cl:2][C:1]([H:3])([H:4])[H:5]",
                {(1, 0): {("[Cl:1]-[C:2]", "EP")}},
            ),
            # Check that two bond charge vsites with different names can be applied
            # to the same atoms
            (
                [
                    VirtualSiteMocking.bond_charge_parameter("[Cl:1]-[C:2]", name="EP"),
                    VirtualSiteMocking.bond_charge_parameter("[Cl:1]-[C:2]", name="LP"),
                ],
                "[Cl:1][C:2]([H:3])([H:4])[H:5]",
                {(0, 1): {("[Cl:1]-[C:2]", "EP"), ("[Cl:1]-[C:2]", "LP")}},
            ),
            # Check that a bond charge vsite can be applied to a symmetric moiety
            (
                [VirtualSiteMocking.bond_charge_parameter("[C:1]#[C:2]")],
                "[H:1][C:2]#[C:3][Cl:4]",
                {(1, 2): {("[C:1]#[C:2]", "EP")}, (2, 1): {("[C:1]#[C:2]", "EP")}},
            ),
            # Check that a monovalent lone pair vsite can be applied to a relatively symmetric moiety
            (
                [VirtualSiteMocking.monovalent_parameter("[O:1]=[C:2]-[*:3]")],
                "[O:2]=[C:1]([H:3])[Cl:4]",
                {
                    (1, 0, 2): {("[O:1]=[C:2]-[*:3]", "EP")},
                    (1, 0, 3): {("[O:1]=[C:2]-[*:3]", "EP")},
                },
            ),
            # Check that monovalent lone pair vsites override correctly
            (
                [
                    VirtualSiteMocking.monovalent_parameter("[O:1]=[C:2]-[*:3]"),
                    VirtualSiteMocking.monovalent_parameter("[O:1]=[C:2]-[Cl:3]"),
                ],
                "[O:2]=[C:1]([H:3])[Cl:4]",
                {(1, 0, 3): {("[O:1]=[C:2]-[Cl:3]", "EP")}},
            ),
            # Check that a single-particle in-plane divalent lone pair vsite
            # can be applied
            (
                [
                    VirtualSiteMocking.divalent_parameter(
                        "[H:2][O:1][H:3]", match="once", angle=0.0 * unit.degree
                    )
                ],
                "[H:1][O:2][H:3]",
                {(1, 0, 2): {("[H:2][O:1][H:3]", "EP")}},
            ),
            # Check that a two-particle out-of-plane divalent lone pair vsite can be applied
            (
                [
                    VirtualSiteMocking.divalent_parameter(
                        "[H:2][O:1][H:3]",
                        match="all_permutations",
                        angle=30.0 * unit.degree,
                    )
                ],
                "[H:1][O:2][H:3]",
                {
                    (1, 0, 2): {("[H:2][O:1][H:3]", "EP")},
                    (1, 2, 0): {("[H:2][O:1][H:3]", "EP")},
                },
            ),
            # Check that a single-particle in-plane divalent lone pair site can
            # be applied to a symmetric molecule
            (
                [
                    VirtualSiteMocking.divalent_parameter(
                        "[*:2]-[N:1]~[*:3]", match="once", angle=0.0 * unit.degree
                    )
                ],
                "[H:1][N:2]=[N:3][H:4]",
                {
                    (1, 0, 2): {("[*:2]-[N:1]~[*:3]", "EP")},
                    (2, 3, 1): {("[*:2]-[N:1]~[*:3]", "EP")},
                },
            ),
            # Check that a two-particle out-of-plane divalent lone pair vsite
            # can be applied to a symmetric molecule
            (
                [
                    VirtualSiteMocking.divalent_parameter(
                        "[*:2]~[N:1]~[*:3]",
                        match="all_permutations",
                        angle=30.0 * unit.degree,
                    )
                ],
                "[H:1][N:2]=[N:3][H:4]",
                {
                    (1, 0, 2): {("[*:2]~[N:1]~[*:3]", "EP")},
                    (1, 2, 0): {("[*:2]~[N:1]~[*:3]", "EP")},
                    (2, 1, 3): {("[*:2]~[N:1]~[*:3]", "EP")},
                    (2, 3, 1): {("[*:2]~[N:1]~[*:3]", "EP")},
                },
            ),
            # Check that a trivalent lone pair vsite can be applied to a simple molecule
            (
                [VirtualSiteMocking.trivalent_parameter("[N:1]([H:2])([H:3])[H:4]")],
                "[N:2]([H:1])([H:3])[H:4]",
                {(1, 0, 2, 3): {("[N:1]([H:2])([H:3])[H:4]", "EP")}},
            ),
        ],
    )
    def test_find_matches(
        self,
        parameters: List[VirtualSiteHandler.VirtualSiteType],
        smiles: str,
        expected_matches: Dict[Tuple[int, ...], List[Tuple[str, str]]],
    ):
        molecule = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)

        handler = VirtualSiteHandler(version="0.3")

        for parameter in parameters:
            handler.add_parameter(parameter=parameter)

        matches = handler.find_matches(molecule.to_topology(), unique=False)

        matched_smirks = defaultdict(set)

        for match_list in matches.values():
            for match in match_list:
                matched_smirks[match.environment_match.topology_atom_indices].add(
                    (match.parameter_type.smirks, match.parameter_type.name)
                )

        assert {**matched_smirks} == expected_matches

    def test_find_matches_multiple_molecules(self):
        topology = Topology.from_molecules(
            [
                Molecule.from_mapped_smiles("[Cl:2][C:1]([H:3])([H:4])[H:5]"),
                Molecule.from_mapped_smiles("[O:2]=[C:1]([H:3])[F:4]"),
            ]
        )

        handler = VirtualSiteHandler(version="0.3")

        handler.add_parameter(
            parameter=VirtualSiteMocking.bond_charge_parameter("[Cl:1]-[C:2]")
        )
        handler.add_parameter(
            parameter=VirtualSiteMocking.monovalent_parameter("[O:1]=[C:2]-[*:3]")
        )

        matches = handler.find_matches(topology, unique=False)

        matched_smirks = defaultdict(set)

        for match_key in matches:
            match_list: List = matches[match_key]

            for match in match_list:
                matched_smirks[match.environment_match.topology_atom_indices].add(
                    (match.parameter_type.smirks, match.parameter_type.name)
                )

        expected_matches = {
            (1, 0): {("[Cl:1]-[C:2]", "EP")},
            (6, 5, 7): {("[O:1]=[C:2]-[*:3]", "EP")},
            (6, 5, 8): {("[O:1]=[C:2]-[*:3]", "EP")},
        }
        assert {**matched_smirks} == expected_matches

    @pytest.mark.parametrize(
        "query_parameter, query_key, expected_index",
        [
            (
                VirtualSiteMocking.monovalent_parameter("[*:1][*:2][*:3]", name="LP"),
                None,
                2,
            ),
            (None, ("BondCharge", "[*:1][*:2]", "EP"), 3),
            (None, ("BondCharge", "[*:1][*:2]", "LP"), 1),
            (None, ("MonovalentLonePair", "[*:1][*:2][*:3]", "EP"), 0),
            (None, ("MonovalentLonePair", "[*:1][*:2][*:3]", "LP"), 2),
            (None, ("DivalentLonePair", "[*:1][*:2][*:3]", "LP"), None),
        ],
    )
    def test_index_of_parameter(self, query_parameter, query_key, expected_index):
        handler = VirtualSiteHandler(version="0.3")

        for parameter in [
            VirtualSiteMocking.monovalent_parameter("[*:1][*:2][*:3]", name="EP"),
            VirtualSiteMocking.bond_charge_parameter("[*:1][*:2]", name="LP"),
            VirtualSiteMocking.monovalent_parameter("[*:1][*:2][*:3]", name="LP"),
            VirtualSiteMocking.bond_charge_parameter("[*:1][*:2]", name="EP"),
        ]:
            handler.add_parameter(parameter=parameter)

        assert handler._index_of_parameter(query_parameter, query_key) == expected_index

    def test_invalid_num_charge_increments(self):
        with pytest.raises(
            SMIRNOFFSpecError,
            match="'BondCharge' virtual sites expect exactly 2 charge increments,",
        ):
            VirtualSiteHandler.VirtualSiteType(
                type="BondCharge",
                smirks="[*]-[*]",
                name="EP",
                charge_increment=[0.1, 0.2, 0.3] * unit.elementary_charge,
                match="all_permutations",
                distance=2.0 * unit.angstrom,
            )


class TestLibraryChargeHandler:
    def test_create_library_charge_handler(self):
        """Test creation of an empty LibraryChargeHandler"""
        LibraryChargeHandler(skip_version_check=True)

    def test_library_charge_type_wrong_num_charges(self):
        """Ensure that an error is raised if a LibraryChargeType is initialized with a different number of
        tagged atoms and charges"""
        LibraryChargeHandler.LibraryChargeType(
            smirks="[#6:1]-[#7:2]",
            charge1=0.1 * unit.elementary_charge,
            charge2=-0.1 * unit.elementary_charge,
        )

        LibraryChargeHandler.LibraryChargeType(
            smirks="[#6:1]-[#7:2]-[#6]",
            charge1=0.1 * unit.elementary_charge,
            charge2=-0.1 * unit.elementary_charge,
        )

        with pytest.raises(
            SMIRNOFFSpecError,
            match="initialized with unequal number of tagged atoms and charges",
        ):
            LibraryChargeHandler.LibraryChargeType(
                smirks="[#6:1]-[#7:2]",
                charge1=0.05 * unit.elementary_charge,
                charge2=0.05 * unit.elementary_charge,
                charge3=-0.1 * unit.elementary_charge,
            )

        with pytest.raises(
            SMIRNOFFSpecError,
            match="initialized with unequal number of tagged atoms and charges",
        ):
            LibraryChargeHandler.LibraryChargeType(
                smirks="[#6:1]-[#7:2]-[#6]",
                charge1=0.05 * unit.elementary_charge,
                charge2=0.05 * unit.elementary_charge,
                charge3=-0.1 * unit.elementary_charge,
            )

        with pytest.raises(
            SMIRNOFFSpecError,
            match="initialized with unequal number of tagged atoms and charges",
        ):
            LibraryChargeHandler.LibraryChargeType(
                smirks="[#6:1]-[#7:2]-[#6]", charge1=0.05 * unit.elementary_charge
            )

    def test_library_charge_type_from_molecule(self):
        mol = Molecule.from_smiles("CCO")

        with pytest.raises(ValueError, match="missing partial"):
            LibraryChargeHandler.LibraryChargeType.from_molecule(mol)

        mol.partial_charges = numpy.linspace(-0.4, 0.4, 9) * unit.elementary_charge

        library_charges = LibraryChargeHandler.LibraryChargeType.from_molecule(mol)

        assert isinstance(library_charges, LibraryChargeHandler.LibraryChargeType)
        assert library_charges.smirks == mol.to_smiles(mapped=True)
        assert library_charges.charge == [*mol.partial_charges]


class TestChargeIncrementModelHandler:
    def test_create_charge_increment_model_handler(self):
        """Test creation of ChargeIncrementModelHandlers"""
        handler = ChargeIncrementModelHandler(skip_version_check=True)
        assert handler.number_of_conformers == 1
        assert handler.partial_charge_method == "AM1-Mulliken"
        handler = ChargeIncrementModelHandler(
            skip_version_check=True, number_of_conformers=10
        )
        handler = ChargeIncrementModelHandler(
            skip_version_check=True, number_of_conformers=1
        )
        handler = ChargeIncrementModelHandler(
            skip_version_check=True, number_of_conformers="10"
        )
        handler = ChargeIncrementModelHandler(
            skip_version_check=True, number_of_conformers=0
        )
        handler = ChargeIncrementModelHandler(
            skip_version_check=True, number_of_conformers="0"
        )
        with pytest.raises(TypeError):
            handler = ChargeIncrementModelHandler(
                skip_version_check=True, number_of_conformers=None
            )
        with pytest.raises(SMIRNOFFSpecError):
            handler = ChargeIncrementModelHandler(
                skip_version_check=True, n_conformers=[10]
            )
        handler = ChargeIncrementModelHandler(
            skip_version_check=True, partial_charge_method="AM1-Mulliken"
        )
        handler = ChargeIncrementModelHandler(
            skip_version_check=True, partial_charge_method="Gasteiger"
        )
        handler = ChargeIncrementModelHandler(
            skip_version_check=True, partial_charge_method=None
        )

    def test_charge_increment_model_handler_getters_setters(self):
        """Test ChargeIncrementModelHandler getters and setters"""
        handler = ChargeIncrementModelHandler(skip_version_check=True)
        assert handler.number_of_conformers == 1
        assert handler.partial_charge_method == "AM1-Mulliken"
        handler.number_of_conformers = 2
        assert handler.number_of_conformers == 2
        handler.number_of_conformers = "3"
        assert handler.number_of_conformers == 3
        with pytest.raises(ValueError):
            handler.number_of_conformers = "string that can't be cast to int"

    def test_charge_increment_model_handlers_are_compatible(self):
        """Test creation of ChargeIncrementModelHandlers"""
        handler1 = ChargeIncrementModelHandler(skip_version_check=True)
        handler2 = ChargeIncrementModelHandler(skip_version_check=True)
        handler1.check_handler_compatibility(handler2)

        handler3 = ChargeIncrementModelHandler(
            skip_version_check=True, number_of_conformers="9"
        )
        with pytest.raises(IncompatibleParameterError):
            handler1.check_handler_compatibility(handler3)

    def test_charge_increment_type_wrong_num_increments(self):
        """Ensure that an error is raised if a ChargeIncrementType is initialized with a different number of
        tagged atoms and chargeincrements"""
        ChargeIncrementModelHandler.ChargeIncrementType(
            smirks="[#6:1]-[#7:2]",
            charge_increment1=0.1 * unit.elementary_charge,
            charge_increment2=-0.1 * unit.elementary_charge,
        )

        ChargeIncrementModelHandler.ChargeIncrementType(
            smirks="[#6:1]-[#7:2]-[#6]",
            charge_increment1=0.1 * unit.elementary_charge,
            charge_increment2=-0.1 * unit.elementary_charge,
        )

        with pytest.raises(
            SMIRNOFFSpecError,
            match="an invalid combination of tagged atoms and charge increments",
        ):
            ChargeIncrementModelHandler.ChargeIncrementType(
                smirks="[#6:1]-[#7:2]",
                charge_increment1=0.05 * unit.elementary_charge,
                charge_increment2=0.05 * unit.elementary_charge,
                charge_increment3=-0.1 * unit.elementary_charge,
            )

        with pytest.raises(
            SMIRNOFFSpecError,
            match="an invalid combination of tagged atoms and charge increments",
        ):
            ChargeIncrementModelHandler.ChargeIncrementType(
                smirks="[#6:1]-[#7:2]-[#6]",
                charge_increment1=0.05 * unit.elementary_charge,
                charge_increment2=0.05 * unit.elementary_charge,
                charge_increment3=-0.1 * unit.elementary_charge,
            )

        ChargeIncrementModelHandler.ChargeIncrementType(
            smirks="[#6:1]-[#7:2]-[#6]",
            charge_increment1=0.05 * unit.elementary_charge,
        )

    def test_charge_increment_one_ci_missing(self):
        """Test creating a chargeincrement parameter with a missing value"""
        ChargeIncrementModelHandler.ChargeIncrementType(
            smirks="[*:1]-[*:2]",
            charge_increment=[0.1 * unit.elementary_charge],
        )

        ChargeIncrementModelHandler.ChargeIncrementType(
            smirks="[*:1]-[*:2]",
            charge_increment=[
                0.1 * unit.elementary_charge,
                -0.1 * unit.elementary_charge,
            ],
        )


class TestGBSAHandler:
    def test_create_default_gbsahandler(self):
        """Test creation of an empty GBSAHandler, with all default attributes"""
        gbsa_handler = GBSAHandler(skip_version_check=True)
        assert gbsa_handler.gb_model == "OBC1"
        assert gbsa_handler.solvent_dielectric == 78.5
        assert gbsa_handler.solute_dielectric == 1
        assert gbsa_handler.sa_model == "ACE"
        assert gbsa_handler.surface_area_penalty == 5.4 * _cal_mol_a2
        assert gbsa_handler.solvent_radius == 1.4 * unit.angstrom

    def test_gbsahandler_setters(self):
        """Test creation of an empty GBSAHandler, with all default attributes"""
        gbsa_handler = GBSAHandler(skip_version_check=True)

        gbsa_handler.gb_model = "OBC2"
        gbsa_handler.gb_model = "HCT"
        gbsa_handler.gb_model = "OBC1"
        with pytest.raises(SMIRNOFFSpecError):
            gbsa_handler.gb_model = "Something invalid"

        gbsa_handler.solvent_dielectric = 50.0
        gbsa_handler.solvent_dielectric = "50.0"
        with pytest.raises(ValueError):
            gbsa_handler.solvent_dielectric = "string that can not be cast to float"

        gbsa_handler.solute_dielectric = 2.5
        gbsa_handler.solute_dielectric = "3.5"
        with pytest.raises(ValueError):
            gbsa_handler.solute_dielectric = "string that can not be cast to float"

        gbsa_handler.sa_model = "ACE"

        # NOTE -- Right now, the SMIRNOFF spec will implicitly assume these are the same.
        gbsa_handler.sa_model = None
        gbsa_handler.sa_model = "None"

        with pytest.raises(TypeError):
            gbsa_handler.sa_model = "Invalid SA option"

        gbsa_handler.surface_area_penalty = (
            1.23 * unit.kilocalorie / unit.mole / unit.nanometer**2
        )
        with pytest.raises(IncompatibleUnitError):
            gbsa_handler.surface_area_penalty = (
                1.23 * unit.degree / unit.mole / unit.nanometer**2
            )

        gbsa_handler.solvent_radius = 300 * unit.femtometer
        with pytest.raises(IncompatibleUnitError):
            gbsa_handler.solvent_radius = 3000 * unit.radian

    def test_gbsahandlers_are_compatible(self):
        """
        Test the check_handler_compatibility function of GBSAHandler
        """
        gbsa_handler_1 = GBSAHandler(skip_version_check=True)
        gbsa_handler_2 = GBSAHandler(skip_version_check=True)

        # Perform a check which should pass
        gbsa_handler_1.check_handler_compatibility(gbsa_handler_2)

        # Perform a check which should fail
        gbsa_handler_3 = GBSAHandler(
            skip_version_check=True, solvent_radius=1.3 * unit.angstrom
        )
        with pytest.raises(
            IncompatibleParameterError, match="Difference between 'solvent_radius' "
        ):
            gbsa_handler_1.check_handler_compatibility(gbsa_handler_3)


class TestParameterTypeReExports:
    def test_parametertype_reexports(self):
        params_module = openff.toolkit.typing.engines.smirnoff.parameters

        def subclass_attrs(obj, classinfo):
            """Iterate over members of ``obj`` that are concrete, public subclasses of ``classinfo``"""
            return filter(
                lambda nv: (  # (name, value)
                    isclass(nv[1])
                    and issubclass(nv[1], classinfo)
                    and not isabstract(nv[1])
                    and not nv[0].startswith("_")
                ),
                vars(obj).items(),
            )

        for _, paramhandler in subclass_attrs(params_module, ParameterHandler):
            for paramtype_name, paramtype in subclass_attrs(
                paramhandler, ParameterType
            ):
                assert paramtype_name in vars(params_module), (
                    f"ParameterType {paramtype_name!r} is "
                    f"not re-exported from parameters module"
                )
                assert vars(params_module)[paramtype_name] is paramtype, (
                    f"Exported attribute parameters.{paramtype_name} "
                    f"does not match ParameterType {paramtype_name!r}"
                )
                assert paramtype_name in params_module.__all__, (
                    f"ParameterType {paramtype_name!r} "
                    f"missing from parameters.__all__"
                )


# TODO: test_(all attributes of all ParameterTypes)
# TODO: test_add_parameter_fractional_bondorder
# TODO: test_get_indexed_attrib
# TODO: test_set_unitbearing_attrib (requires implementing __getattr__ and __setattr__)
# TODO: test_parametertype_unit_getattr
# TODO: test_parametertype_unit_setattr
# TODO: test_optional_attribs
# TODO: test_optional_indexed_attribs
# TODO: test_(X)handler_compatibility, where X is all handlers
