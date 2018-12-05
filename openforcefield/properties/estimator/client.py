#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Property estimator client side API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import logging

from openforcefield.properties.estimator.runner import PropertyCalculationRunner


# =============================================================================================
# Property Estimator
# =============================================================================================

class PropertyEstimator(object):
    """
    The object responsible for requesting a set of properties
    be calculated by the low-level property calculation backend,
    and for analysing the performance of the parameters.
    """

    @staticmethod
    def compute_properties(data_set, parameter_set, worker_threads=1):
        """
        Submit the property and parameter set for calculation.

        Parameters
        ----------
        data_set : PropertyDataSet
            The set of properties to attempt to compute.
        parameter_set : ParameterSet
            The OpenFF parameter set to use for the calculations.
        worker_threads : int
            The number of worker threads to calculate the properties on.
        """

        if data_set is None or parameter_set is None:

            raise ValueError('Both a data set and parameter set must be '
                             'present to compute physical properties.')

        # In principle this method would simply push all properties and
        # to the backend which will decide what to do with them.

        # For now, just create the backend manually on the local device.
        calculation_runner = PropertyCalculationRunner(worker_threads)

        # In practice such a complicated runner will need to report back
        # detailed diagnostics of what ran and why, and what if anything
        # went wrong.
        calculated_properties = calculation_runner.run(data_set, parameter_set)

        return calculated_properties

    @staticmethod
    def _store_properties_in_hierarchy(original_set):
        """Refactor a property list into a hierarchy of substance->state->type.

        Parameters
        ----------
        original_set : dict(str, list(PhysicalProperty))
            The set of properties to refactor.
        """
        property_hierarchy = {}

        for substance_tag in original_set:

            for calculated_property in original_set[substance_tag]:

                if substance_tag not in property_hierarchy:
                    property_hierarchy[substance_tag] = {}

                state_tag = hash(calculated_property.thermodynamic_state)

                if state_tag not in property_hierarchy[substance_tag]:
                    property_hierarchy[substance_tag][state_tag] = {}

                if calculated_property.type not in property_hierarchy[substance_tag][state_tag]:
                    property_hierarchy[substance_tag][state_tag][calculated_property.type] = {}

                property_hierarchy[substance_tag][state_tag][calculated_property.type] = calculated_property

        return property_hierarchy

    @staticmethod
    def produce_calculation_report(measured_data_set, calculated_data_set):
        """
        Produce a report detailing how well a measured and calculated data
        set match.

        Parameters
        ----------
        measured_data_set : PhysicalPropertyDataSet
            The set of measured properties to compare against.
        calculated_data_set : CalculatedPropertySet
            The set of calculated properties to analyse.
        """
        measured_properties = PropertyEstimator._store_properties_in_hierarchy(
            measured_data_set.properties)

        calculated_properties = PropertyEstimator._store_properties_in_hierarchy(
            calculated_data_set.properties)

        for substance in calculated_properties:

            for state in calculated_properties[substance]:

                if len(calculated_properties[substance][state]) <= 0:
                    continue

                state_string = next(iter(calculated_properties[substance][state].values())).thermodynamic_state

                logging.info('PROPERTIES FOR ' + substance + ' AT ' + str(state_string))

                for property_type in calculated_properties[substance][state]:

                    measured_property = measured_properties[substance][state][property_type]
                    calculated_property = calculated_properties[substance][state][property_type]

                    logging.info('Property: ' + str(property_type) +
                                 ' Measured: ' + str(measured_property.value) +
                                 '(' + str(measured_property.uncertainty) + ')' +
                                 ' Calculated: ' + str(calculated_property.value) +
                                 '(' + str(calculated_property.uncertainty) + ')')

        return
