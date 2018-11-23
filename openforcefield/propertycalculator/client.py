#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Property calculator client side API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import logging

from openforcefield.propertycalculator.runner import PropertyCalculationRunner


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
    def compute_properties(data_set, parameter_set, worker_threads = 1):
        """
        Submit the property and parameter set for calculation.

        Parameters
        ----------
        data_set : PropertyDataSet
            The set of properties to attempt to compute.
        parameter_set : ParameterSet
            The OpenFF parameter set to use for the calculations.
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
    def produce_calculation_report(measured_data_set, calculated_data_set):
        """
        Produce a report detailing how well a measured and calculated data
        set match.

        Parameters
        ----------
        measured_data_set : PropertyDataSet
            The set of measured properties to compare against.
        calculated_data_set : ParameterSet
            The set of calculated properties to analyse.
        """

        # TODO: The way properties are stored needs to be refactored to be uniform.
        measured_properties = {}

        for measured_property in measured_data_set.measured_properties:

            substance_tag = measured_property.substance.to_tag()

            if substance_tag not in measured_properties:
                measured_properties[substance_tag] = {}

            state_tag = measured_property.thermodynamic_state.to_tag()

            if state_tag not in measured_properties[substance_tag]:
                measured_properties[substance_tag][state_tag] = {}

            measured_properties[substance_tag][state_tag][measured_property.type] = measured_property

        calculated_properties = {}

        for substance_tag in calculated_data_set.properties:

            for calculated_property in calculated_data_set.properties[substance_tag]:

                if substance_tag not in calculated_properties:
                    calculated_properties[substance_tag] = {}

                state_tag = calculated_property.thermodynamic_state.to_tag()

                if state_tag not in calculated_properties[substance_tag]:
                    calculated_properties[substance_tag][state_tag] = {}

                calculated_properties[substance_tag][state_tag][calculated_property.type] = calculated_property

        for substance in calculated_properties:

            for state in calculated_properties[substance]:

                logging.info('PROPERTIES FOR ' + substance + ' AT ' + state)

                for property_type in calculated_properties[substance][state]:

                    measured_property = measured_properties[substance][state][property_type]
                    calculated_property = calculated_properties[substance][state][property_type]

                    logging.info('Property: ' + str(property_type) +
                                 ' Measured: ' + str(measured_property.value) +
                                 '(' + str(measured_property.uncertainty) + ')' +
                                 ' Calculated: ' + str(calculated_property.value) +
                                 '(' + str(calculated_property.uncertainty) + ')')

        return
