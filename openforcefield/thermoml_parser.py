import re
import copy
import numpy as np
import simtk.unit as unit
from .data.properties import thermoml_schema
from .substances import Substance, Mixture
from .thermodynamics import ThermodynamicState

from typing import List

def determine_method_name(e_name, s_name):
    if e_name is None or ("other" in e_name.lower() and s_name is not None):
        return s_name
    elif e_name is not None and s_name is not None:
        print("Ambiguous Output, {} -- {}".format(e_name, s_name))
    return e_name


class Parser(object):

    _QUANTITY_KEYWORDS = ('Mole fraction',)
    _STATE_MAP = {'Pressure, kPa': ("pressure", unit.kilopascals),
                  'Temperature, K': ("temperature", unit.kelvin)}

    def __init__(self, filename):
        """Create a parser object from an XML filename."""
        self.filename = filename
        self.root = thermoml_schema.CreateFromDocument(open(self.filename).read())
        self.compound_num_to_name = {}
        self.compound_name_to_formula = {}
        self.store_compounds()
        self.doi = self.root.Citation.sDOI
    
    def store_compounds(self):
        """Extract and store compounds from a thermoml XML file."""

        for Compound in self.root.Compound:
            nOrgNum = Compound.RegNum.nOrgNum
            sCommonName = Compound.sCommonName[0]
            sFormulaMolec = Compound.sFormulaMolec

            self.compound_num_to_name[nOrgNum] = sCommonName
            self.compound_name_to_formula[sCommonName] = sFormulaMolec

    def _compile_measurement(self):
        output = dict()
        output['doi'] = self.doi
        output['filename'] = self.filename
        return output

    @staticmethod
    def _build_components(pure_or_mix_data):
        # Build the list of compounds
        components = []
        for Component in pure_or_mix_data.Component:
            n_org_num = Component.RegNum.nOrgNum
            components.append(n_org_num)
        return components

    def _component_string_from_ids(self, component_ids: List[int]):
        string_list = []
        for comp_id in component_ids:
            s_common_name = self.compound_num_to_name[comp_id]
            string_list.append(s_common_name)
        return string_list

    def _build_proprty_dict(self, pure_or_mixture_data, components_string, only_emethod=False):
        # Build the list of properties (things we are observing)
        property_dict = {}
        for Property in pure_or_mixture_data.Property:
            n_prop_number = Property.nPropNumber  # Property ID number
            property_content = Property.Property_MethodID.PropertyGroup.orderedContent()[0].value
            e_prop_name = property_content.ePropName  # ASSUMING LENGTH 1
            e_method_name = property_content.eMethodName
            s_method_name = property_content.sMethodName
            # Compound ID number

            compound = self._parse_compounds(Property.Property_MethodID)
            # if n_org_number is None:
            #     compound_name = components_string
            # else:
            #     compound_name = self.compound_name_to_formula[n_org_number]
            if not only_emethod:
                true_method_name = determine_method_name(e_method_name, s_method_name)
            else:
                true_method_name = e_method_name
            e_prop_phase = Property.PropPhaseID[0].ePropPhase  # ASSUMING LENGTH 1
            single_property_dict = {'method': true_method_name,
                                    'type': e_prop_name,
                                    'compound': compound,
                                    'phase': e_prop_phase}
            property_dict[n_prop_number] = single_property_dict
        return property_dict

    @staticmethod
    def _build_constraints(pure_or_mixture_data):
        constraint_dict = {}

        for index, Constraint in enumerate(pure_or_mixture_data.Constraint):
            n_constraint_value = Constraint.nConstraintValue  # Actual value of constraint
            constraint_type = Constraint.ConstraintID.ConstraintType  # Type dict
            assert len(constraint_type.orderedContent()) == 1
            # Content of the constraint type (true constraint name)
            try:
                reg_num = Constraint.SoConstraintIDlvent.RegNum
            except AttributeError:
                reg_num = None
            if reg_num is None:
                compound = None
            else:
                compound = [x.nOrgNum for x in reg_num]
            constraint_type = constraint_type.orderedContent()[0].value
            constraint_data = {'type': constraint_type,
                               'value': n_constraint_value,
                               'compound': compound}
            constraint_dict[index] = constraint_data
        return constraint_dict

    def _build_variables(self, pure_or_mixture_data):
        variable_dict = {}
        for Variable in pure_or_mixture_data.Variable:
            n_var_number = Variable.nVarNumber  # Variable ID
            variable_type = Variable.VariableID.VariableType
            assert len(variable_type.orderedContent()) == 1
            # Assume length 1, haven't found counterexample yet.
            variable_type = variable_type.orderedContent()[0].value
            compound = self._parse_compounds(Variable)
            variable_data = {'type': variable_type,
                             'compound': compound}
            variable_dict[n_var_number] = variable_data
        return variable_dict

    @staticmethod
    def _parse_compounds(node):
        try:
            reg_num = node.RegNum
        except AttributeError:
            reg_num = None
        if reg_num is None:
            compound = None
        else:
            try:
                compound = [x.nOrgNum for x in reg_num]
            except TypeError:
                compound = reg_num.nOrgNum
        return compound

    def _update_state(self, state, *data):
        for d in data:
            for dk, dv in d.items():
                if dv['type'] in self._STATE_MAP:
                    state_variable, correct_unit = self._STATE_MAP[dv['type']]
                    setattr(state, state_variable, dv['value'] * correct_unit)

    def _update_composition(self, composition, *data):
        for d in data:
            for dk, dv in d.items():
                compound = dv['compound']
                # Ensure compound is not none
                if compound is not None and any([(x in dv['type']) for x in self._QUANTITY_KEYWORDS]):
                    if isinstance(compound, int):
                        compound = [compound]
                    for c in compound:
                        composition[c]['Mole fraction'] = dv['value']

    @staticmethod
    def _validate_data(state, composition):
        state_good = False
        if state.pressure is not None and state.pressure is not None:
            state_good = True
        composition_good = False
        net_composition = 0.0
        count_invalid_mol_fracs = 0
        missing_mol_frac = None
        for comp, comp_val in composition.items():
            mole_frac = comp_val['Mole fraction']
            if mole_frac is not None:
                net_composition += mole_frac
            else:
                count_invalid_mol_fracs += 1
                missing_mol_frac = comp_val
        # Missing 1 item means we can deduce the last mole fraction
        if count_invalid_mol_fracs == 1:
            missing_mol_frac['Mole fraction'] = 1-net_composition
        if net_composition == 1.0 and all([(value['Mole fraction'] is not None) for value in composition.values()]):
            # This will only catch all but the MOST rigid data.
            # E.g. solubility in binary fluid: je500991b.xml
            # The data reported are the solubility of the solid in the binary fluid, but only the mole fraction
            # of one of the fluid components is reported. We have to deduce the other one.
            # Also, every component in the given PureOrMixtureData does not have to contribute to the same phase.
            # For instance: in this example, the Variable is on Acetonitrile, but the Methanol is left to be
            # derived for the Solvent purity, pre-measurement, so Xa + Xm = 1, but the urea component is also listed
            # so going only on the component list will lead to undefined mixture.
            # If we have to do context-specific learning, this gets MUCH harder.
            composition_good = True
        if state_good and composition_good:
            return True
        return False

    def parse(self, only_emethod=False):
        """Parse the current XML filename and return a list of measurements."""
        alldata = []
        """
        - Nab the source
        - For each Data entry
          - Build the list of compounds
          - Build the list of properties (measurements)
          - Cycle through constraints (constants)
          - Cycle through Variable (Controlled variables)
          - Cycle through values
            - Get Source
            - Get target compound
            - Get get the value
            - Get the uncertainty
            - Get the State
            - Get the method
            - Get the observable
        """
        # For each collection of data
        for PureOrMixtureData in self.root.PureOrMixtureData:
            components = self._build_components(PureOrMixtureData)
            # Assume we don't know mole_components until further notice

            components_string = "__".join(self._component_string_from_ids(components))
            property_dict = self._build_proprty_dict(PureOrMixtureData, components_string, only_emethod=only_emethod)
            # Temperature and Pressure States determined at value value runtime
            constraint_dict = self._build_constraints(PureOrMixtureData)

            variable_dict = self._build_variables(PureOrMixtureData)

            # Populate empty objects
            constant_state = ThermodynamicState(temperature=None, pressure=None)
            constant_composition = {index: {'Mole fraction': None,
                                            'name': self.compound_num_to_name[index]}
                                    for index in components}

            self._update_state(constant_state, constraint_dict)
            self._update_composition(constant_composition, constraint_dict)

            for NumValues in PureOrMixtureData.NumValues:
                local_composition = copy.deepcopy(constant_composition)
                local_state = copy.deepcopy(constant_state)
                for VariableValue in NumValues.VariableValue:
                    n_var_number = VariableValue.nVarNumber
                    n_var_value = VariableValue.nVarValue
                    variable_data = variable_dict[n_var_number]
                    # Append value to the data
                    variable_data['value'] = n_var_value
                property_collection = {}
                for PropertyValue in NumValues.PropertyValue:
                    # Get property
                    n_prop_number = PropertyValue.nPropNumber
                    n_prop_value = PropertyValue.nPropValue
                    try:
                        uncertainty = PropertyValue.PropUncertainty[0].nStdUncertValue
                    except IndexError:
                        uncertainty = None
                    single_property = {'value': n_prop_value,
                                       'uncertainty': uncertainty,
                                       **property_dict[n_prop_number]}
                    property_collection[n_prop_number] = single_property
                self._update_state(local_state, variable_dict, property_collection)
                self._update_composition(local_composition, variable_dict, property_collection)
                if self._validate_data(local_state, local_composition):
                    output_state = local_state
                    output_mixture = Mixture()
                    for comp_value in local_composition.values():
                        output_mixture.add_component(comp_value['name'],
                                                     mole_fraction=comp_value['Mole fraction'])
                    for single_property in property_collection.values():
                        alldata.append((single_property, output_state, output_mixture))

        return alldata


def count_atoms(formula_string):
    """Parse a chemical formula and return the total number of atoms."""
    element_counts = formula_to_element_counts(formula_string)
    return sum(val for key, val in element_counts.items())


def count_atoms_in_set(formula_string, which_atoms):
    """Parse a chemical formula and return the number of atoms in a set of atoms."""
    element_counts = formula_to_element_counts(formula_string)
    return sum(val for key, val in element_counts.items() if key in which_atoms)


def get_first_entry(cas):
    """If cirpy returns several CAS results, extracts the first one."""
    if type(cas) == type([]):
        cas = cas[0]
    return cas


def formula_to_element_counts(formula_string):
    """Transform a chemical formula into a dictionary of (element, number) pairs."""
    pattern = r'([A-Z][a-z]{0,2}\d*)'
    pieces = re.split(pattern, formula_string)
    data = pieces[1::2]
    _ = filter(None, pieces[0::2])
    pattern2 = r'([A-Z][a-z]{0,2})'
    results = {}
    for piece in data:
        element, number = re.split(pattern2, piece)[1:]
        try:
            number = int(number)
        except ValueError:
            number = 1
        results[element] = number

    return results
