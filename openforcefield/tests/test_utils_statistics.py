"""
Units tests for openforcefield.utils.statistics
"""

import logging

import numpy as np
from pymbar import timeseries

from openforcefield.utils import statistics


def test_statistical_inefficiency():
    """Test the statistical inefficiency calculation utility."""

    data_size = 200000

    random_array = np.random.rand(data_size)
    numpy_vector_array = []

    for i in range(data_size):
        numpy_vector_array.append([random_array[i]])

    A = np.array(numpy_vector_array)

    statistical_inefficiency = statistics.calculate_statistical_inefficiency(A, minimum_samples=3)
    pymbar_statistical_inefficiency = timeseries.statisticalInefficiency(A, mintime=3)

    print('utils: {}, pymbar: {}', statistical_inefficiency, pymbar_statistical_inefficiency)

    assert abs(statistical_inefficiency - pymbar_statistical_inefficiency) < 0.00001
