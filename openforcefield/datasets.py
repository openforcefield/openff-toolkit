#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Substances API.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>
* Levi N. Naden <levi.naden@choderalab.org>

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import copy
import json
import inspect
import pandas as pd
from openforcefield.measurements import MeasuredPhysicalProperty, Source
from openforcefield.thermodynamics import ThermodynamicState
import openforcefield.measurements as measurements

from typing import List, Tuple, Optional, Union

measurement_members = (member for member in inspect.getmembers(measurements, inspect.isclass)
                       if issubclass(member, measurements.MeasurementMethod))


# =============================================================================================
# DATASET
# =============================================================================================

class PhysicalPropertyDataset(object):
    """A dataset of physical property measurements.

    Implements the container API.

    Properties
    ----------

    """
    def __init__(self, dataset: Optional[List[MeasuredPhysicalProperty]]=None):
        """Create a physical property dataset.

        Parameters
        ----------
        dataset : iterable of MeasuredPhysicalProperty's, optional, default=None
            If a dataset is specified, it is copied to create a new dataset.

        """
        self._measurements = list()

        if dataset is not None:
            for measurement in dataset:
                self._measurements.append(copy.deepcopy(measurement))

    def __len__(self):
        return len(self._measurements)

    def __delitem__(self, doi):
        """Delete all items with the given DOI.
        """
        measurement = self[doi]
        self._measurements.remove(measurement)

    def __iter__(self):
        for measurement in self._measurements:
            yield measurement

    def __contains__(self, measurement):
        for measurement2 in self._measurements:
            if measurement == measurement2:
                return True
        return False

    def filter(self, **kwargs):
        """
        Filter the dataset down based on key-word parameters in place

        This is an irreversible filter

        Parameters
        ----------
        kwargs : TBD

        """
        pass

    def to_pandas(self) -> pd.DataFrame:
        """
        Converts the data set to a Pandas DataFrame object

        Returns
        -------
        data_frame : pandas.DataFrame
            Data Frame representation of the current
        """
        return pd.DataFrame

    @classmethod
    def from_pandas(cls, data_frame) -> 'PhysicalPropertyDataset':
        """
        Creates a new PhysicalPropertyDataset from a Pandas DataFrame.
        This is the reverse of the ``to_pandas`` function

        Parameters
        ----------
        data_frame : pandas.DataFrame
            DataFrame representation of the measured property dataset.

        Returns
        -------
        data_set : PhysicalPropertyDataset
            New PhysicalPropertyDataset constructed from the Pandas dataframe
        """
        return cls


    @property
    def dois(self) -> List[str]:
        """Returns a list of the DOI's which compose this data set derived from sources"""
        return ["placeholder"]

    @property
    def dois(self) -> List[str]:
        """Returns a list of the references which compose this data set derived from sources"""
        return ["placeholder"]

    def _collect_sources(self):
        """Gather sources from __init__"""
        pass


# =============================================================================================
# THERMOML DATASET
# =============================================================================================

class ThermoMLDataset(PhysicalPropertyDataset):
    """A dataset of physical property measurements created from a ThermoML dataset.

    Examples
    --------

    For example, we can use the DOI `10.1016/j.jct.2005.03.012` as a key for retrieving the dataset from the ThermoML Archive:

    >>> dataset = ThermoMLDataset('10.1016/j.jct.2005.03.012')

    You can also specify multiple ThermoML Archive keys to create a dataset from multiple ThermoML files:

    >>> thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
    >>> dataset = ThermoMLDataset(thermoml_keys)

    You can see which DOIs contribute to the current `ThermoMLDataset` with the convenience functions:

    >>> thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
    >>> dataset = ThermoMLDataset(thermoml_keys)

    """
    def __init__(self, doi=None, url=None):
        """Retrieve ThermoML files matching specified keys from the specified URL.

        Parameters
        ----------
        doi : str or list of str, optional, default=None
            If specified, ThermoML files with the specified DOI keys will be retrieved
        url : str, optional, default=None
            If specified, this URL (which may point to a local filesystem) will be used instead of the ThermoML Archive.

        """
        self._sources = []  # TODO: Fix this for __getitem__
        observables = []
        compiled = []
        cycle_doi = []
        cycle_url = []
        if isinstance(doi, str):
            doi = [doi]
        elif doi is None:
            doi = []
        cycle_doi.extend(doi)
        if isinstance(url, str):
            url = [url]
        elif url is None:
            url = []
        cycle_url.extend(url)
        for do in cycle_doi:
            self._sources.append(do)
            observables.append(self._process_doi(do))
        for ur in cycle_url:
            self._sources.append(ur)
            observables.extend(self._process_url(ur))
        for observable in observables:
            substance, temperature, pressure, property_name, method_name, value, uncertainty, source = observable
            measurement_method = self._find_method(property_name, method_name)
            source = self._make_source(source)
            if measurement_method is None:
                # No known method found
                continue
            state = ThermodynamicState(temperature, pressure)
            measurement = MeasuredPhysicalProperty(substance,
                                                   state,
                                                   measurement_method,
                                                   value=value,
                                                   uncertainty=uncertainty)
            compiled.append(measurement)
        super().__init__(compiled)

    def __getitem__(self, doi):
        """Return all measurements with the given DOI.
        """
        measurements = list()
        for measurement in self._measurements:
            if measurement.doi == doi:
                measurements.append(measurement)
        if len(measurements) == 0:
            return KeyError("DOI (%s) not found in dataset." % doi)
        return measurements

    def retrieve(self,
                 doi: Optional[Union[str, List[str]]] = None,
                 url: Optional[Union[str, List[str]]] = None,
                 ) -> List[MeasuredPhysicalProperty]:
        """
        Fetch a subset of parameters from a set of DOI(s) and/or URL(s)

        Parameters
        ----------
        doi : str or list of str or None, optional, default None
            DOI source(s) to fetch data from
        url : str or list of str or None, optional, default None

        Returns
        -------
        properties : list of MeasuredPhysicalProperty
            Output measured observables

        """
        # Placeholder return
        return []

    def apply_nist_uncertainties(self,
                                 nist_uncertainties: dict,
                                 adjust_uncertainties: Optional[bool] = True,
                                 discard_outliers: Optional[bool] = True):
        """
        Apply uncertainty corrections from a NIST set. This is an in-place update and relies on the source
        DOI and/or references for accurate corrections.

        NIST has compiled data are in JSON format which must be converted to dict before being passed into this function

        Parameters
        ----------
        nist_uncertainties : dict
            Collection of NIST uncertainty corrections converted from json
        adjust_uncertainties : bool, optional, default True
            TODO: Figure out what this does separate from the function its in (it part of spec)
        discard_outliers: bool, optional, default True
            Discard observables which are statistical outliers from the corrected uncertainties
        """
        pass

    @staticmethod
    def _process_doi(doi: str) -> List[Tuple]:
        """
        Fetch and process the data from a DOI source

        Parameters
        ----------
        doi

        Returns
        -------
        observables : list of tuple
            Observables compiled from source, fetches them all in the source
            Returns a form
                Substance,
                Temperature,
                Pressure,
                PropertyName,
                MethodName,
                Value,
                Uncertainty
                Source,
        """
        pass

    @staticmethod
    def _process_url(doi: str) -> List[Tuple]:
        """
        Fetch and process the data from a URL source

        Parameters
        ----------
        url

        Returns
        -------
        observables : dict of {(Property, MethodName) : (Value, Uncertainty)}
            Observables compiled from source, fetches them all in the source
            Returns a form
                Substance,
                Temperature,
                Pressure,
                PropertyName,
                MethodName,
                Value,
                Uncertainty,
                Source
        """
        pass

    @staticmethod
    def _find_method(property_name, method_name):
        output = None
        for member in measurement_members:
            try:
                x = member()
                if hasattr(x, property_name) and hasattr(x, method_name):
                    output = x
                    break
            except TypeError:
                # Is an abstract class since `member()` will fail
                continue
        return output

    @staticmethod
    def _make_source(source) -> Source:
        """Given a source from DOI/URL, return a source class"""
        return Source(doi="placeholder")

