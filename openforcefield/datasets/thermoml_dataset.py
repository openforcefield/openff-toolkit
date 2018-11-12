#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
An API for importing a ThermoML archive.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org> (original thermoml_parser)
* Levi N. Naden <levi.naden@choderalab.org> (original thermoml_parser)
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

from __future__ import with_statement

import re
import copy

from urllib.error import HTTPError
from urllib.request import urlopen

from enum import Enum, unique

from xml.etree import ElementTree

from openeye import oechem
from openeye import oeiupac

from simtk import unit

from openforcefield.thermodynamics import ThermodynamicState
# from openforcefield.measurements import MeasuredPhysicalProperty

from openforcefield.properties import PropertyPhase, PropertyType

from .property_dataset import PhysicalPropertyDataSet


# =============================================================================================
# ThermoMLConstraintType
# =============================================================================================

@unique
class ThermoMLConstraintType(Enum):

    Undefined = 0
    Temperature = 1
    Pressure = 2

    @staticmethod
    def from_string(string):

        constraint_type = ThermoMLConstraintType.Undefined

        if string == 'Temperature, K':
            constraint_type = ThermoMLConstraintType.Temperature
        elif string == 'Pressure, kPa':
            constraint_type = ThermoMLConstraintType.Pressure

        return constraint_type


# =============================================================================================
# ThermoMLConstraint
# =============================================================================================

class ThermoMLConstraint:

    def __init__(self):

        self.type = ThermoMLConstraintType.Undefined
        self.value = 0.0

    @classmethod
    def from_string(cls, constraint_type_string, unitless_value):

        return_value = cls()

        return_value.type = ThermoMLConstraintType.from_string(constraint_type_string)
        return_value.value = cls.add_unit_to_value(return_value.type, unitless_value)

        return return_value

    @classmethod
    def from_enum(cls, constraint_type_enum, unitless_value):

        return_value = cls()

        return_value.type = constraint_type_enum
        return_value.value = cls.add_unit_to_value(return_value.type, unitless_value)

        return return_value

    @staticmethod
    def add_unit_to_value(constraint_type, unitless_value):

        value = unitless_value

        if constraint_type is ThermoMLConstraintType.Temperature:
            value *= unit.kelvin
        elif constraint_type is ThermoMLConstraintType.Pressure:
            value *= unit.kilo * unit.pascal

        return value


# =============================================================================================
# ThermoMLVariableDefinition
# =============================================================================================

class ThermoMLVariableDefinition:

    def __init__(self):

        self.index = -1
        self.type = ThermoMLConstraintType.Undefined


# =============================================================================================
# ThermoMLComponent
# =============================================================================================

class ThermoMLComponent:

    def __init__(self):

        self.smiles = None
        self.index = -1

    # TODO: This functionality should really exist at the
    #       toolkit level in an OEChem independent way.
    @staticmethod
    def smiles_from_identifier(identifier):

        temp_molecule = oechem.OEMol()

        parse_smiles_options = oechem.OEParseSmilesOptions(quiet=True)

        if oechem.OEParseSmiles(temp_molecule, identifier, parse_smiles_options) is True:
            # Should make sure all smiles are OEChem consistent.
            return oechem.OEMolToSmiles(temp_molecule)

        if oeiupac.OEParseIUPACName(temp_molecule, identifier) is False:
            return None

        return oechem.OEMolToSmiles(temp_molecule)

    @classmethod
    def from_xml_node(cls, node, namespace):

        # Gather up all possible identifiers
        identifier_nodes = node.findall('ThermoML:sSmiles', namespace)

        # identifier_nodes.extend(node.findall('ThermoML:sStandardInChI', namespace))
        identifier_nodes.extend(node.findall('ThermoML:sIUPACName', namespace))
        identifier_nodes.extend(node.findall('ThermoML:sCommonName', namespace))

        if len(identifier_nodes) == 0:
            # convert common name to smiles
            raise RuntimeError('A ThermoML:Compound node does not have a proper identifier')

        identifier = cls.smiles_from_identifier(identifier_nodes[0].text)

        if identifier is None:

            print('WARNING: The component identifier ' + identifier_nodes[0].text +
                  ' could not be converted to a SMILES pattern and will be skipped')

            return None

        index_node = node.find('./ThermoML:RegNum//ThermoML:nOrgNum', namespace)

        if index_node is None or index_node.text.isdigit is False:
            raise RuntimeError('A ThermoML:Compound has an invalid (or non-existent) index')

        component_index = int(index_node.text)

        return_value = cls()

        return_value.smiles = identifier
        return_value.index = component_index

        return return_value


# =============================================================================================
# ThermoMLProperty
# =============================================================================================

class ThermoMLProperty(object):

    def __init__(self):

        self.index = None
        self.phase = PropertyPhase.Undefined

        self.type = PropertyType.Undefined
        self.method_name = None

        self.value = 0.0
        self.uncertainty = 0.0

        self.property_uncertainty_definitions = {}
        self.combined_uncertainty_definitions = {}

    @staticmethod
    def property_string_to_enum(string):

        string_split = string.split(',')
        return_value = PropertyType.Undefined

        property_string = string_split[0] if len(string_split) > 0 else None

        # Do we actually need to check the unit here? Probably not..
        # unit_string = string_split[0] if len(string_split) > 1 else None

        if property_string == 'Mass density':
            return_value = PropertyType.MassDensity

        return return_value

    @staticmethod
    def phase_string_to_enum(string):

        phase_string = string.lower()
        phase = PropertyPhase.Undefined

        if (phase_string.find('liquid') and not phase_string.find('crystal')) or \
                phase_string.find('fluid') or phase_string.find('solution'):

            phase = PropertyPhase.Liquid

        elif phase_string.find('crystal') and not phase_string.find('liquid'):

            phase = PropertyPhase.Solid

        elif phase_string.find('gas'):

            phase = PropertyPhase.Gas

        return phase

    @staticmethod
    def extract_uncertainty_definitions(node, namespace,
                                        property_uncertainty_definitions,
                                        combined_uncertainty_definitions):

        property_nodes = node.findall('ThermoML:CombinedUncertainty', namespace)

        for property_node in property_nodes:

            if property_node is None:
                continue

            uncertainty_definition = ThermoMLCombinedUncertainty.from_xml(property_node, namespace)

            if uncertainty_definition is None:
                continue

            combined_uncertainty_definitions[uncertainty_definition.index] = uncertainty_definition

        property_nodes = node.findall('ThermoML:PropUncertainty', namespace)

        for property_node in property_nodes:

            if property_node is None:
                continue

            uncertainty_definition = ThermoMLPropertyUncertainty.from_xml(property_node, namespace)

            if uncertainty_definition is None:
                continue

            property_uncertainty_definitions[uncertainty_definition.index] = uncertainty_definition

    @classmethod
    def from_xml_node(cls, node, namespace):

        # Gather up all possible identifiers
        index_node = node.find('ThermoML:nPropNumber', namespace)

        if index_node is None or index_node.text.isdigit is False:
            raise RuntimeError('A ThermoML:Property has an invalid (or non-existent) index')

        property_index = int(index_node.text)

        phase_node = node.find('./ThermoML:PropPhaseID//ThermoML:ePropPhase', namespace)

        if phase_node is None:

            raise RuntimeError('The ePropPhase property is missing.')

        phase = ThermoMLProperty.phase_string_to_enum(phase_node.text)

        if phase is PropertyPhase.Undefined:

            print('WARNING: A property was measured in an unsupported phase (' +
                  phase_node.text + ') and will be skipped.')

            return None

        # TODO: BEWARE - Property-MethodID also has a RegNum entry which describes
        # which component is referred to if the property is based on one of
        # the components... e.g. mass fraction of component 2
        #
        # This is ignored for now however...

        property_group_node = node.find('./ThermoML:Property-MethodID//ThermoML:PropertyGroup//*', namespace)

        property_name_node = property_group_node.find('./ThermoML:ePropName', namespace)
        method_name_node = property_group_node.find('./ThermoML:sMethodName', namespace)

        if method_name_node is None:

            method_name_node = property_group_node.find('./eMethodName', namespace)

        if method_name_node is None or property_name_node is None:

            raise RuntimeError('A property does not have a name / method entry.')

        property_type = ThermoMLProperty.property_string_to_enum(property_name_node.text)

        if property_type is PropertyType.Undefined:

            print('WARNING: An unsupported property was found (' + property_name_node.text +
                  ') and will be skipped.')

            return None

        return_value = cls()

        return_value.index = property_index
        return_value.phase = phase

        return_value.type = property_type
        return_value.method_name = method_name_node.text

        property_uncertainty_definitions = {}
        combined_uncertainty_definitions = {}

        cls.extract_uncertainty_definitions(node, namespace,
                                            property_uncertainty_definitions,
                                            combined_uncertainty_definitions)

        return_value.combined_uncertainty_definitions = combined_uncertainty_definitions
        return_value.property_uncertainty_definitions = property_uncertainty_definitions

        return return_value

    def set_value(self, value, uncertainty):

        if self.type is PropertyType.MassDensity:

            self.value = value / unit.meters ** 3 * unit.kilogram
            self.uncertainty = uncertainty / unit.meters ** 3 * unit.kilogram


# =============================================================================================
# ThermoMLPropertyUncertainty
# =============================================================================================

class ThermoMLPropertyUncertainty:

    prefix = ''

    def __init__(self):

        self.index = -1
        self.coverage_factor = None

    @classmethod
    def from_xml(cls, node, namespace):

        coverage_factor_node = node.find('ThermoML:n' + cls.prefix + 'CoverageFactor', namespace)
        confidence_node = node.find('ThermoML:n' + cls.prefix + 'UncertLevOfConfid', namespace)

        coverage_factor = None

        if coverage_factor_node is not None:
            coverage_factor = float(coverage_factor_node.text)
        elif confidence_node is not None and confidence_node.text == '95':
            # TODO: Is this actually correct?
            coverage_factor = 2
        else:
            return None

        index_node = node.find('ThermoML:n' + cls.prefix + 'UncertAssessNum', namespace)

        if index_node is None or not index_node.text.isdigit():
            raise RuntimeError('A ThermoML:PropUncert does not have a valid index')

        index = int(index_node.text)

        return_value = cls()

        return_value.coverage_factor = coverage_factor
        return_value.index = index

        return return_value


# =============================================================================================
# ThermoMLCombinedUncertainty
# =============================================================================================

class ThermoMLCombinedUncertainty(ThermoMLPropertyUncertainty):
    prefix = 'Comb'


# =============================================================================================
# ThermoMLPureOrMixtureData
# =============================================================================================

class ThermoMLPureOrMixtureData:

    @staticmethod
    def extract_component_indices(node, namespace, components):

        component_nodes = node.findall('ThermoML:Component', namespace)
        component_indices = []

        # Figure out which components are going to be associated with
        # the property entries.
        for component_node in component_nodes:

            index_node = component_node.find('./ThermoML:RegNum//ThermoML:nOrgNum', namespace)

            if index_node is None or index_node.text.isdigit() is False:

                raise RuntimeError('A ThermoML:Component entry within a ThermoML:PureOrMixtureData entry'
                                   ' has an invalid / non-existent index')

            component_index = int(index_node.text)

            if component_index not in components:

                print('WARNING: A PureOrMixtureData which depends on an unsupported compound has been skipped')
                return None

            if component_index in component_indices:
                raise RuntimeError('A ThermoML:PureOrMixtureData entry defines the same component twice')

            component_indices.append(component_index)

        return component_indices

    @staticmethod
    def extract_property_definitions(node, namespace):

        property_nodes = node.findall('ThermoML:Property', namespace)
        properties = {}

        for property_node in property_nodes:

            property_definition = ThermoMLProperty.from_xml_node(property_node, namespace)

            if property_definition is None:
                continue

            if property_definition.index in properties:

                raise RuntimeError('A ThermoML data set contains two '
                                   'properties with the same index')

            properties[property_definition.index] = property_definition

        return properties

    @staticmethod
    def extract_global_constraints(node, namespace):

        constraint_nodes = node.findall('ThermoML:Constraint', namespace)
        constraints = {}

        for constraint_node in constraint_nodes:

            type_node = constraint_node.find('.//ThermoML:ConstraintType/*', namespace)
            value = float(constraint_node.find('./ThermoML:nConstraintValue', namespace).text)

            constraint = ThermoMLConstraint.from_string(type_node.text, value)

            if constraint.type is ThermoMLConstraintType.Undefined:

                print('WARNING: An unsupported constraint exists upon the measurement (' +
                      type_node.text + ') - the offending PureOrMixtureData entry will be skipped')

                return None

            constraints[constraint.type] = constraint

        return constraints

    @staticmethod
    def extract_variable_definitions(node, namespace):

        variable_nodes = node.findall('ThermoML:Variable', namespace)
        variables = {}

        for variable_node in variable_nodes:

            index_node = variable_node.find('ThermoML:nVarNumber', namespace)
            type_node = variable_node.find('.//ThermoML:VariableType/*', namespace)

            variable = ThermoMLVariableDefinition()

            variable.index = int(index_node.text)
            variable.type = ThermoMLConstraintType.from_string(type_node.text)

            if variable.type is ThermoMLConstraintType.Undefined:

                print('WARNING: An unsupported variable exists upon the measurement (' +
                      type_node.text + ') and will be skipped')

                continue

            variables[variable.index] = variable

        return variables

    @staticmethod
    def calculate_uncertainty(node, namespace, property_definition):

        # Look for a standard uncertainty..
        uncertainty_node = node.find('.//ThermoML:nCombStdUncertValue', namespace)

        if uncertainty_node is None:
            uncertainty_node = node.find('.//ThermoML:nStdUncertValue', namespace)

        # We have found a std. uncertainty
        if uncertainty_node is not None:
            return float(uncertainty_node.text)

        # Try to calculate uncertainty from a coverage factor if present
        if len(property_definition.combined_uncertainty_definitions) == 0 and \
           len(property_definition.property_uncertainty_definitions) == 0:

            return None

        combined = len(property_definition.combined_uncertainty_definitions) > 0

        prefix = ThermoMLCombinedUncertainty.prefix if combined \
            else ThermoMLPropertyUncertainty.prefix

        index_node = node.find('.//ThermoML:n' + prefix + 'ExpandUncertValue', namespace)
        expanded_uncertainty_node = node.find('.//ThermoML:n' + prefix + 'ExpandUncertValue', namespace)

        if index_node is None or expanded_uncertainty_node is None:
            return None

        expanded_uncertainty = float(expanded_uncertainty_node.text)
        index = int(index_node.text)

        if combined and index not in property_definition.combined_uncertainty_definitions:
            return None

        if not combined and index not in property_definition.property_uncertainty_definitions:
            return None

        divisor = property_definition.combined_uncertainty_definitions[index].coverage_factor if combined \
            else property_definition.property_uncertainty_definitions[index].coverage_factor

        return expanded_uncertainty / divisor

    @staticmethod
    def extract_measured_properties(node, namespace,
                                    property_definitions,
                                    global_constraints,
                                    variable_definitions):

        value_nodes = node.findall('ThermoML:NumValues', namespace)

        measured_properties = []

        # Each value_node corresponds to one MeasuredProperty
        for value_node in value_nodes:

            constraints = copy.deepcopy(global_constraints)

            # First extract the values of any variable constraints
            variable_nodes = value_node.findall('ThermoML:VariableValue', namespace)

            skip_entry = False

            for variable_node in variable_nodes:

                variable_index = int(variable_node.find('./ThermoML:nVarNumber', namespace).text)

                if variable_index not in variable_definitions:

                    # The property was constrained by an unsupported variable and
                    # so will be skipped for now.
                    skip_entry = True
                    break

                variable_definition = variable_definitions[variable_index]
                variable_value = float(variable_node.find('ThermoML:nVarValue', namespace).text)

                if variable_definition.type in constraints:

                    raise RuntimeError('A ThermoML:NumValues entry attempts to define the value of ' +
                                       variable_definition.type + ' twice')

                "Convert the 'variable' into a full constraint entry"
                constraint = ThermoMLConstraint.from_enum(variable_definition.type, variable_value)
                constraints[constraint.type] = constraint

            if skip_entry:
                continue

            temperature = constraints[ThermoMLConstraintType.Temperature].value
            pressure = constraints[ThermoMLConstraintType.Pressure].value

            thermodynamic_state = ThermodynamicState(temperature, pressure)

            # Now extract the actual values of the measured properties, and their
            # uncertainties

            property_nodes = value_node.findall('ThermoML:PropertyValue', namespace)

            for property_node in property_nodes:

                property_index = int(property_node.find('./ThermoML:nPropNumber', namespace).text)

                if property_index not in property_definitions:

                    # Most likely the property was dropped earlier due to an unsupported phase / type
                    continue

                property_definition = property_definitions[property_index]

                uncertainty = ThermoMLPureOrMixtureData.calculate_uncertainty(property_node,
                                                                              namespace,
                                                                              property_definition)

                if uncertainty is None:

                    print('WARNING: A property (' + str(property_definition.type) +
                          ') without uncertainties was skipped')

                    continue

                measured_property = copy.deepcopy(property_definition)

                property_value_node = property_node.find('.//ThermoML:nPropValue', namespace)

                measured_property.set_value(float(property_value_node.text),
                                            float(uncertainty))

                measured_properties.append(measured_property)

        # By this point we now have the measured properties and the thermodynamic state
        # they were measured at - convert to the standardised classes.

        return measured_properties

    @staticmethod
    def from_xml_node(node, namespace, components):

        # Figure out which components are going to be associated with
        # the property entries.
        component_indices = ThermoMLPureOrMixtureData.extract_component_indices(node, namespace, components)

        if component_indices is None:
            # Most likely this entry depended on a non-parsable compound
            # and will be skipped entirely
            return None

        if len(component_indices) == 0:

            print('WARNING: A PureOrMixtureData entry with no components was skipped.')
            return None

        # Extract property definitions - values come later!
        property_definitions = ThermoMLPureOrMixtureData.extract_property_definitions(node, namespace)

        if len(property_definitions) == 0:

            print('WARNING: A PureOrMixtureData entry with no properties was skipped. ' +
                  'Most likely this entry only contained unsupported properties.')

            return None

        # Extract any constraints on the system e.g pressure, temperature
        global_constraints = ThermoMLPureOrMixtureData.extract_global_constraints(node, namespace)

        if global_constraints is None:
            return None

        # Extract any variables set on the system e.g pressure, temperature
        # Only the definition entry and not the value of the variable is extracted
        variable_definitions = ThermoMLPureOrMixtureData.extract_variable_definitions(node, namespace)

        if len(global_constraints) == 0 and len(variable_definitions) == 0:

            print('WARNING: A PureOrMixtureData entry with no constraints was skipped.')
            return None

        measured_properties = ThermoMLPureOrMixtureData.extract_measured_properties(node, namespace,
                                                                                    property_definitions,
                                                                                    global_constraints,
                                                                                    variable_definitions)

        # TODO: Construct thermodynamic state.
        # TODO: Construct construct Mixture... Need to speak to Levi about challenges..
        # TODO: Construct the final measured properties.

        return measured_properties


# =============================================================================================
# ThermoMLDataSet
# =============================================================================================

class ThermoMLDataSet(PhysicalPropertyDataSet):
    """A dataset of physical property measurements created from a ThermoML dataset.

    Examples
    --------

    For example, we can use the DOI `10.1016/j.jct.2005.03.012` as a key
    for retrieving the dataset from the ThermoML Archive:

    >>> dataset = ThermoMLDataSet('10.1016/j.jct.2005.03.012')

    You can also specify multiple ThermoML Archive keys to create a dataset from multiple ThermoML files:

    >>> thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
    >>> dataset = ThermoMLDataSet(thermoml_keys)

    You can see which DOIs contribute to the current `ThermoMLDataset` with the convenience functions:

    >>> thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
    >>> dataset = ThermoMLDataSet(thermoml_keys)

    """

    def __init__(self):

        super().__init__()

    @classmethod
    def from_doi(cls, doi):

        # E.g https://trc.nist.gov/ThermoML/10.1016/j.jct.2016.12.009.xml
        doi_url = 'https://trc.nist.gov/ThermoML/' + doi + '.xml'

        return cls.from_url(doi_url)

    @classmethod
    def from_url(cls, url):

        # TODO : Load file from URL
        data = ''

        try:
            with urlopen(url) as response:

                data = response.read()

        except HTTPError as error:

            print('WARNING: No ThermoML file could not be found at ' + url)
            return None

        return cls.from_xml(data)

    @classmethod
    def from_xml(cls, xml):

        root_node = ElementTree.fromstring(xml)

        if root_node is None:
            print('WARNING: The XML document could not be parsed.')
            return None

        # Extract the namespace that will prefix all type names
        namespace_string = re.search('\{.*\}', root_node.tag).group(0)[1:-1]
        namespace = {'ThermoML': namespace_string}

        return_value = ThermoMLDataSet()
        components = {}

        # Extract the base components present in the xml file
        for node in root_node.findall('ThermoML:Compound', namespace):

            component = ThermoMLComponent.from_xml_node(node, namespace)

            if component is None:
                continue

            if component.index in components:
                raise RuntimeError('A ThermoML data set contains two '
                                   'components with the same index')

            components[component.index] = component

        # Pull out any and all properties in the file.
        for node in root_node.findall('ThermoML:PureOrMixtureData', namespace):

            properties = ThermoMLPureOrMixtureData.from_xml_node(node, namespace, components)

            if properties is None or len(properties) == 0:
                continue

            return_value._measured_properties.extend(properties)

        return return_value
