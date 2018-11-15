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

import logging

# import json
# import ujson as json
import _pickle as cPickle

from urllib.error import HTTPError
from urllib.request import urlopen

from enum import IntEnum, unique

from xml.etree import ElementTree

from openeye import oechem
from openeye import oeiupac

from simtk import unit

from openforcefield.measurements import Source, MeasuredPhysicalProperty
from openforcefield.thermodynamics import ThermodynamicState
from openforcefield.substances import Mixture

from openforcefield.properties import PropertyPhase, PropertyType

from .property_dataset import PhysicalPropertyDataSet


# =============================================================================================
# Temporary Helper Methods
# =============================================================================================

def unit_from_thermoml_string(full_string):
    """A non-ideal way to convert a string to a simtk.unit.Unit"""

    full_string_split = full_string.split(',')

    unit_string = full_string_split[1] if len(full_string_split) > 1 else ''
    unit_string = unit_string.strip()

    if unit_string == 'K':
        return unit.kelvin
    elif unit_string == 'kPa':
        return unit.kilo * unit.pascal
    elif unit_string == 'kg/m3':
        return unit.kilogram / unit.meter**3
    elif unit_string == 'mol/kg':
        return unit.mole / unit.kilogram
    elif unit_string == 'mol/dm3':
        return unit.mole / unit.decimeter**3
    elif len(unit_string) == 0:
        return None
    else:
        raise NotImplementedError('The unit (' + unit_string + ') is not currently supported')


def phase_from_thermoml_string(string):

    phase_string = string.lower().strip()
    phase = PropertyPhase.Undefined

    if (phase_string == 'liquid' or
        phase_string.find('fluid') >= 0 or phase_string.find('solution') >= 0):

        phase = PropertyPhase.Liquid

    elif phase_string.find('crystal') >= 0 and not phase_string.find('liquid') >= 0:

        phase = PropertyPhase.Solid

    elif phase_string.find('gas') >= 0:

        phase = PropertyPhase.Gas

    return phase


# =============================================================================================
# ThermoMLConstraintType
# =============================================================================================

@unique
class ThermoMLConstraintType(IntEnum):

    """
    An enum containing the supported types of ThermoML constraints
    """

    Undefined            = 0x00
    Temperature          = 0x01
    Pressure             = 0x02
    ComponentComposition = 0x04
    SolventComposistion  = 0x08

    @staticmethod
    def from_node(node):

        """
            Converts either a ConstraintType or VariableType xml node to a ThermoMLConstraintType.
        """

        constraint_type = ThermoMLConstraintType.Undefined

        if node.tag.find('eTemperature') >= 0 and node.text == 'Temperature, K':
            constraint_type = ThermoMLConstraintType.Temperature
        elif node.tag.find('ePressure') >= 0 and node.text == 'Pressure, kPa':
            constraint_type = ThermoMLConstraintType.Pressure
        elif node.tag.find('eComponentComposition') >= 0 and node.text == 'Mole fraction':
            constraint_type = ThermoMLConstraintType.ComponentComposition
        elif node.tag.find('eSolventComposition') >= 0 and node.text == 'Mole fraction':
            constraint_type = ThermoMLConstraintType.SolventComposistion
        else:
            logging.warning(node.tag + '->' + node.text + ' is an unsupported constraint type.')

        return constraint_type


# =============================================================================================
# ThermoMLConstraint
# =============================================================================================

class ThermoMLConstraint:

    def __init__(self):

        self.type = ThermoMLConstraintType.Undefined
        self.value = 0.0

        self.solvents = []

        # Describes which compound the variable acts upon.
        self.compound_index = None

    @classmethod
    def from_node(cls, constraint_node, namespace):

        # Extract the xml nodes.
        type_node = constraint_node.find('.//ThermoML:ConstraintType/*', namespace)
        value_node = constraint_node.find('./ThermoML:nConstraintValue', namespace)

        solvent_index_nodes = constraint_node.find('./ThermoML:Solvent//nOrgNum', namespace)
        compound_index_node = constraint_node.find('./ThermoML:ConstraintID/ThermoML:RegNum/*', namespace)

        value = float(value_node.text)
        # Determine what the default unit for this variable should be.
        unit_type = unit_from_thermoml_string(type_node.text)

        return_value = cls()

        return_value.type = ThermoMLConstraintType.from_node(type_node)
        return_value.value = unit.Quantity(value, unit_type)

        if compound_index_node is not None:
            return_value.compound_index = int(compound_index_node.text)

        if solvent_index_nodes is not None:
            for solvent_index_node in solvent_index_nodes:
                return_value.solvents.append(int(solvent_index_node))

        return None if return_value.type is ThermoMLConstraintType.Undefined else return_value

    @classmethod
    def from_variable(cls, variable, value):

        return_value = cls()

        return_value.type = variable.type
        return_value.compound_index = variable.compound_index

        return_value.solvents.extend(variable.solvents)

        return_value.value = value

        return return_value


# =============================================================================================
# ThermoMLVariableDefinition
# =============================================================================================

class ThermoMLVariableDefinition:

    def __init__(self):

        self.index = -1
        self.type = ThermoMLConstraintType.Undefined

        self.solvents = []
        self.default_unit = None

        # Describes which compound the variable acts upon.
        self.compound_index = None

    @classmethod
    def from_node(cls, variable_node, namespace):

        # Extract the xml nodes.
        type_node = variable_node.find('.//ThermoML:VariableType/*', namespace)
        index_node = variable_node.find('ThermoML:nVarNumber', namespace)

        solvent_index_nodes = variable_node.find('./ThermoML:Solvent//nOrgNum', namespace)
        compound_index_node = variable_node.find('./ThermoML:VariableID/ThermoML:RegNum/*', namespace)

        return_value = cls()

        return_value.default_unit = unit_from_thermoml_string(type_node.text)

        return_value.index = int(index_node.text)
        return_value.type = ThermoMLConstraintType.from_node(type_node)

        if compound_index_node is not None:
            return_value.compound_index = int(compound_index_node.text)

        if solvent_index_nodes is not None:
            for solvent_index_node in solvent_index_nodes:
                return_value.solvents.append(int(solvent_index_node))

        return None if return_value.type is ThermoMLConstraintType.Undefined else return_value


# =============================================================================================
# ThermoMLCompound
# =============================================================================================

class ThermoMLCompound:

    def __init__(self):

        self.smiles = None

        self.index = -1

    # TODO: SMILES from name should exist at the toolkit level
    # in an OEChem independent way.
    @staticmethod
    def smiles_from_identifier(identifier):

        if identifier is None:
            return None

        temp_molecule = oechem.OEMol()
        parse_smiles_options = oechem.OEParseSmilesOptions(quiet=True)

        smiles = None

        if oechem.OEParseSmiles(temp_molecule, identifier, parse_smiles_options) is True:
            # Should make sure all smiles are OEChem consistent.
            smiles = oechem.OEMolToSmiles(temp_molecule)
        elif oeiupac.OEParseIUPACName(temp_molecule, identifier) is True:
            smiles = oechem.OEMolToSmiles(temp_molecule)

        if smiles is None:
            return None

        return smiles

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

        if identifier_nodes[0].text is None:
            # A pathological case where no name is actually given...
            return None

        identifier = cls.smiles_from_identifier(identifier_nodes[0].text)

        if identifier is None:

            logging.warning('The compound identifier ' + identifier_nodes[0].text +
                            ' could not be converted to a SMILES pattern and will be skipped')

            return None

        index_node = node.find('./ThermoML:RegNum/*', namespace)

        if index_node is None:
            raise RuntimeError('A ThermoML:Compound has a non-existent index')

        compound_index = int(index_node.text)

        return_value = cls()

        return_value.smiles = identifier

        return_value.index = compound_index

        return return_value


# =============================================================================================
# ThermoMLProperty
# =============================================================================================

class ThermoMLProperty(MeasuredPhysicalProperty):

    def __init__(self):

        super().__init__()

        self.index = None

        self.solvents = []

        self.property_uncertainty_definitions = {}
        self.combined_uncertainty_definitions = {}

        self.default_unit = None

    @staticmethod
    def property_string_to_enum(string):

        string_split = string.split(',')
        return_value = PropertyType.Undefined

        property_string = string_split[0] if len(string_split) > 0 else None

        if property_string == 'Mass density':
            return_value = PropertyType.MassDensity
        elif property_string == 'Relative permittivity at zero frequency':
            return_value = PropertyType.DielectricConstant

        return return_value

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

        property_index = int(index_node.text)

        phase_node = node.find('./ThermoML:PropPhaseID//ThermoML:ePropPhase', namespace)
        phase = PropertyPhase.Undefined

        if phase_node is not None:
            phase = phase_from_thermoml_string(phase_node.text)

        if phase == PropertyPhase.Undefined:

            # TODO: For now we just hope that the property defines the phase.
            # This needs to be better supported however.
            logging.warning('A property was measured in an unsupported phase (' +
                            phase_node.text + ') and will be skipped.')

            return None

        # TODO: Property->RegNum is currently ignored
        # Describes which compound is referred to if the property is based on one of
        # the compounds... e.g. mass fraction of compound 2

        property_group_node = node.find('./ThermoML:Property-MethodID//ThermoML:PropertyGroup//*', namespace)

        property_name_node = property_group_node.find('./ThermoML:ePropName', namespace)
        method_name_node = property_group_node.find('./ThermoML:eMethodName', namespace)

        if method_name_node is None:
            method_name_node = property_group_node.find('./ThermoML:sMethodName', namespace)

        if method_name_node is None or property_name_node is None:
            raise RuntimeError('A property does not have a name / method entry.')

        property_type = ThermoMLProperty.property_string_to_enum(property_name_node.text)

        if property_type == PropertyType.Undefined:

            logging.warning('An unsupported property was found (' + property_name_node.text +
                            ') and will be skipped.')

            return None

        return_value = cls()

        return_value.index = property_index
        return_value.phase = phase

        return_value.default_unit = unit_from_thermoml_string(property_name_node.text)

        return_value.type = property_type
        return_value.method_name = method_name_node.text

        property_uncertainty_definitions = {}
        combined_uncertainty_definitions = {}

        cls.extract_uncertainty_definitions(node, namespace,
                                            property_uncertainty_definitions,
                                            combined_uncertainty_definitions)

        return_value.combined_uncertainty_definitions = combined_uncertainty_definitions
        return_value.property_uncertainty_definitions = property_uncertainty_definitions

        solvent_index_nodes = node.find('./ThermoML:Solvent//nOrgNum', namespace)

        if solvent_index_nodes is not None:
            for solvent_index_node in solvent_index_nodes:
                return_value.solvents.append(int(solvent_index_node))

        return return_value

    def set_value(self, value, uncertainty):

        value_quantity = value
        uncertainty_quantity = uncertainty

        if not isinstance(value_quantity, unit.Quantity):
            value_quantity = unit.Quantity(value, self.default_unit)
        if not isinstance(uncertainty_quantity, unit.Quantity):
            uncertainty_quantity = unit.Quantity(uncertainty, self.default_unit)

        self.value = value_quantity
        self.uncertainty = uncertainty_quantity


# =============================================================================================
# ThermoMLPropertyUncertainty
# =============================================================================================

class ThermoMLPropertyUncertainty:

    # Reduce code redundancy by reusing this class for
    # both property and combined uncertainties.
    prefix = ''

    def __init__(self):

        self.index = -1
        self.coverage_factor = None

    @classmethod
    def from_xml(cls, node, namespace):

        coverage_factor_node = node.find('ThermoML:n' + cls.prefix + 'CoverageFactor', namespace)
        confidence_node = node.find('ThermoML:n' + cls.prefix + 'UncertLevOfConfid', namespace)

        coverage_factor = None

        # TODO: Does a confidence interval definitely equate to a coverage of 2?
        if coverage_factor_node is not None:
            coverage_factor = float(coverage_factor_node.text)
        elif confidence_node is not None and confidence_node.text == '95':
            coverage_factor = 2
        else:
            return None

        index_node = node.find('ThermoML:n' + cls.prefix + 'UncertAssessNum', namespace)
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
    def extract_compound_indices(node, namespace, compounds):

        component_nodes = node.findall('ThermoML:Component', namespace)
        compound_indices = []

        # Figure out which compounds are going to be associated with
        # the property entries.
        for component_node in component_nodes:

            index_node = component_node.find('./ThermoML:RegNum/*', namespace)

            compound_index = int(index_node.text)

            if compound_index not in compounds:

                logging.warning('A PureOrMixtureData entry depends on an '
                                'unsupported compound and has been ignored')

                return None

            if compound_index in compound_indices:

                raise RuntimeError('A ThermoML:PureOrMixtureData entry defines the same compound twice')

            compound_indices.append(compound_index)

        return compound_indices

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
    def extract_global_constraints(node, namespace, compounds):

        constraint_nodes = node.findall('ThermoML:Constraint', namespace)
        constraints = []

        for constraint_node in constraint_nodes:

            constraint = ThermoMLConstraint.from_node(constraint_node, namespace)

            if constraint is None or constraint.type is ThermoMLConstraintType.Undefined:
                logging.warning('An unsupported constraint has been ignored.')
                return None

            if constraint.compound_index is not None and \
               constraint.compound_index not in compounds:

                logging.warning('A constraint exists upon a non-existent compound and will be ignored.')
                return None

            if (constraint.type is ThermoMLConstraintType.SolventComposistion or
               constraint.type is ThermoMLConstraintType.ComponentComposition) and constraint.compound_index is None:

                logging.warning('An unsupported constraint has been ignored - composition constraints'
                                'need to have a corresponding compound_index.')

            constraints.append(constraint)

        return constraints

    @staticmethod
    def extract_variable_definitions(node, namespace, compounds):

        variable_nodes = node.findall('ThermoML:Variable', namespace)
        variables = {}

        for variable_node in variable_nodes:

            variable = ThermoMLVariableDefinition.from_node(variable_node, namespace)

            if variable is None or variable.type is ThermoMLConstraintType.Undefined:

                logging.warning('An unsupported variable has been ignored.')

                continue

            if variable.compound_index is not None and \
               variable.compound_index not in compounds:

                logging.warning('A constraint exists upon a non-existent compound and will be ignored.')

                continue

            if (variable.type is ThermoMLConstraintType.SolventComposistion or
               variable.type is ThermoMLConstraintType.ComponentComposition) and variable.compound_index is None:

                logging.warning('An unsupported variable has been ignored - composition variables'
                                'need to have a corresponding compound_index.')

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

        index_node = None

        if combined:
            index_node = node.find('./ThermoML:CombinedUncertainty/ThermoML:nCombUncertAssessNum', namespace)
        else:
            index_node = node.find('./ThermoML:PropUncertainty/ThermoML:nUncertAssessNum', namespace)

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
    def build_mixture(measured_property, constraints, compounds):

        if measured_property.phase != PropertyPhase.Liquid:

            logging.warning('Only properties measured in the liquid phase '
                            'are currently supported (' + str(measured_property.phase) + ').')

            return None

        mixture = Mixture()

        # Handle the easy case where the system has to
        # be pure
        if len(compounds) == 1:

            mixture.add_component(next(iter(compounds.values())).smiles, 1.0)
            return mixture

        mol_fractions = {}

        number_of_constraints = 0
        total_mol_fraction = 0.0

        for constraint in constraints:

            if constraint.type != ThermoMLConstraintType.ComponentComposition and \
               constraint.type != ThermoMLConstraintType.SolventComposistion:

                continue

            mol_fractions[constraint.compound_index] = constraint.value / unit.dimensionless

            total_mol_fraction += mol_fractions[constraint.compound_index]
            number_of_constraints += 1

        if number_of_constraints == len(compounds) and \
            abs(total_mol_fraction - 1.0) > 0.00001:

            raise RuntimeError('The total mol fraction does not add to 1.0')

        elif number_of_constraints > len(compounds):
            raise RuntimeError('There are more concentration constraints than componenents.')

        elif number_of_constraints < len(compounds) - 1:
            raise RuntimeError('There are too many unknown mole fractions.')

        elif number_of_constraints == len(compounds) - 1:

            for compound_index in compounds:

                if compound_index in mol_fractions:
                    continue

                mol_fractions[compound_index] = 1.0 - total_mol_fraction

        else:

            raise RuntimeError('Unexpected edge case..')

        for compound_index in compounds:

            if mol_fractions[compound_index] < 0.00001:
                continue

            compound = compounds[compound_index]
            mixture.add_component(compound.smiles, mol_fractions[compound_index])

        return mixture

    @staticmethod
    def extract_measured_properties(node, namespace,
                                    property_definitions,
                                    global_constraints,
                                    variable_definitions,
                                    compounds):

        value_nodes = node.findall('ThermoML:NumValues', namespace)

        measured_properties = []

        # Each value_node corresponds to one MeasuredProperty
        for value_node in value_nodes:

            constraints = []

            temperature_constraint = None
            pressure_constraint = None

            for global_constraint in global_constraints:

                # constraint = copy.deepcopy(global_constraint)
                # constraint = json.loads(json.dumps(global_constraint))
                constraint = cPickle.loads(cPickle.dumps(global_constraint, -1))
                constraints.append(constraint)

                if constraint.type == ThermoMLConstraintType.Temperature:
                    temperature_constraint = constraint
                elif constraint.type == ThermoMLConstraintType.Pressure:
                    pressure_constraint = constraint

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
                value_as_quantity = unit.Quantity(variable_value, variable_definition.default_unit)

                "Convert the 'variable' into a full constraint entry"
                constraint = ThermoMLConstraint.from_variable(variable_definition, value_as_quantity)
                constraints.append(constraint)

                if constraint.type == ThermoMLConstraintType.Temperature:
                    temperature_constraint = constraint
                elif constraint.type == ThermoMLConstraintType.Pressure:
                    pressure_constraint = constraint

            if skip_entry:
                continue

            # Extract the thermodynamic state that the property was measured at.
            if temperature_constraint is None or \
               pressure_constraint is None:

                logging.warning('A property did not state the T / p it was measured'
                                'at and will be ignored.')
                continue

            temperature = temperature_constraint.value
            pressure = pressure_constraint.value

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

                    logging.warning('A property (' + str(property_definition.type) +
                                    ') without uncertainties was ignored')

                    continue

                # measured_property = copy.deepcopy(property_definition)
                # measured_property = json.loads(json.dumps(property_definition))
                measured_property = cPickle.loads(cPickle.dumps(property_definition, -1))
                measured_property.thermodynamic_state = thermodynamic_state

                property_value_node = property_node.find('.//ThermoML:nPropValue', namespace)

                measured_property.set_value(float(property_value_node.text),
                                            float(uncertainty))

                mixture = ThermoMLPureOrMixtureData.build_mixture(measured_property,
                                                                  constraints,
                                                                  compounds)

                if mixture is None:

                    logging.warning('Could not build a mixture for a property (' +
                                    str(property_definition.type) + ').')

                    continue

                measured_property.substance = mixture

                measured_properties.append(measured_property)

        # By this point we now have the measured properties and the thermodynamic state
        # they were measured at converted to standardised classes.
        return measured_properties

    @staticmethod
    def from_xml_node(node, namespace, compounds):

        # Figure out which compounds are going to be associated with
        # the property entries.
        compound_indices = ThermoMLPureOrMixtureData.extract_compound_indices(node, namespace, compounds)

        if compound_indices is None:
            # Most likely this entry depended on a non-parsable compound
            # and will be skipped entirely
            return None

        if len(compound_indices) == 0:

            logging.warning('A PureOrMixtureData entry with no compounds was ignored.')
            return None

        # Extract property definitions - values come later!
        property_definitions = ThermoMLPureOrMixtureData.extract_property_definitions(node, namespace)

        if len(property_definitions) == 0:

            logging.warning('A PureOrMixtureData entry with no properties was ignored. ' +
                            'Most likely this entry only contained unsupported properties.')

            return None

        phase_nodes = node.findall('./ThermoML:PhaseID/ThermoML:ePhase', namespace)

        all_phases = None

        for phase_node in phase_nodes:

            phase = phase_from_thermoml_string(phase_node.text)

            if phase == PropertyPhase.Undefined:
                # TODO: For now we just hope that the property defines the phase.
                # This needs to be better supported however.
                logging.warning('A property was measured in an unsupported phase (' +
                                phase_node.text + ') and will be skipped.')

                return None

            all_phases = phase if all_phases is None else all_phases | phase

        for property_index in property_definitions:

            all_phases |= property_definitions[property_index].phase
            property_definitions[property_index].phase |= all_phases

        # Extract any constraints on the system e.g pressure, temperature
        global_constraints = ThermoMLPureOrMixtureData.extract_global_constraints(node, namespace, compounds)

        if global_constraints is None:
            return None

        # Extract any variables set on the system e.g pressure, temperature
        # Only the definition entry and not the value of the variable is extracted
        variable_definitions = ThermoMLPureOrMixtureData.extract_variable_definitions(node, namespace, compounds)

        if len(global_constraints) == 0 and len(variable_definitions) == 0:

            logging.warning('A PureOrMixtureData entry with no constraints was ignored.')
            return None

        used_compounds = {}

        for compound_index in compounds:

            if compound_index not in compound_indices:
                continue

            used_compounds[compound_index] = compounds[compound_index]

        measured_properties = ThermoMLPureOrMixtureData.extract_measured_properties(node, namespace,
                                                                                    property_definitions,
                                                                                    global_constraints,
                                                                                    variable_definitions,
                                                                                    used_compounds)

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

    >>> dataset = ThermoMLDataSet.from_doi_list('10.1016/j.jct.2005.03.012')

    You can also specify multiple ThermoML Archive keys to create a dataset from multiple ThermoML files:

    >>> thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
    >>> dataset = ThermoMLDataSet(*thermoml_keys)

    """

    def __init__(self):

        super().__init__()

    @classmethod
    def from_doi_list(cls, *doi_list):

        return_value = None

        for doi in doi_list:

            # E.g https://trc.nist.gov/ThermoML/10.1016/j.jct.2016.12.009.xml
            doi_url = 'https://trc.nist.gov/ThermoML/' + doi + '.xml'

            data_set = cls._from_url(doi_url, Source(doi=doi))

            if data_set is None or len(data_set.measured_properties) == 0:
                continue

            if return_value is None:
                return_value = data_set
            else:
                return_value.merge(data_set)

        return return_value

    @classmethod
    def from_url_list(cls, *url_list):

        return_value = None

        for url in url_list:

            data_set = cls._from_url(url)

            if data_set is None or len(data_set.measured_properties) == 0:
                continue

            if return_value is None:
                return_value = data_set
            else:
                return_value.merge(data_set)

        return return_value

    @classmethod
    def _from_url(cls, url, source = None):

        if source is None:
            source = Source(reference=url)

        return_value = None

        try:
            with urlopen(url) as response:

                return_value = cls.from_xml(response.read(), source)

        except HTTPError as error:

            logging.warning('WARNING: No ThermoML file could not be found at ' + url)

        return return_value

    @classmethod
    def from_file_list(cls, *file_list):

        return_value = None
        counter = 0

        for file in file_list:

            data_set = cls._from_file(file)

            print('Reading file ', counter, ' of ', len(file_list), ' (', file, ')')
            counter += 1

            if data_set is None or len(data_set.measured_properties) == 0:
                continue

            if return_value is None:
                return_value = data_set
            else:
                return_value.merge(data_set)

        return return_value

    @classmethod
    def _from_file(cls, path):

        source = Source(reference=path)
        return_value = None

        try:

            with open(path) as file:

                return_value = ThermoMLDataSet.from_xml(file.read(), source)

        except FileNotFoundError as error:

            logging.warning('No ThermoML file could not be found at ' + path)

        return return_value

    @classmethod
    def from_xml(cls, xml, source):

        root_node = ElementTree.fromstring(xml)

        if root_node is None:
            logging.warning('The ThermoML XML document could not be parsed.')
            return None

        if root_node.tag.find('DataReport') < 0:
            logging.warning('The ThermoML XML document does not contain the expected root node.')
            return None

        # Extract the namespace that will prefix all type names
        namespace_string = re.search('{.*\}', root_node.tag).group(0)[1:-1]
        namespace = {'ThermoML': namespace_string}

        return_value = ThermoMLDataSet()
        compounds = {}

        # Extract the base compounds present in the xml file
        for node in root_node.findall('ThermoML:Compound', namespace):

            compound = ThermoMLCompound.from_xml_node(node, namespace)

            if compound is None:
                continue

            if compound.index in compounds:
                raise RuntimeError('A ThermoML data set contains two '
                                   'compounds with the same index')

            compounds[compound.index] = compound

        # Pull out any and all properties in the file.
        for node in root_node.findall('ThermoML:PureOrMixtureData', namespace):

            properties = ThermoMLPureOrMixtureData.from_xml_node(node, namespace, compounds)

            if properties is None or len(properties) == 0:
                continue

            return_value._measured_properties.extend(properties)

        for measured_property in return_value._measured_properties:
            # Set the source of the data.
            measured_property.source = source

        return_value._sources.append(source)

        return return_value
