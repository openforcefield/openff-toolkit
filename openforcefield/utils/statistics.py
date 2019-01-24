# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
A set of utilities for performing statistical analysis on a time series.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org> (original pymbar implementation)
* Simon Boothroyd <simon.boothroyd@choderalab.org> (port to higher dimensions)
"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import math

import numpy as np
from pymbar.utils import ParameterError


# =============================================================================================
# Utilities
# =============================================================================================

def calculate_statistical_inefficiency(time_series, minimum_samples=3):
    """Calculates the statistical inefficiency of a time series.

    Notes
    -----
    The statistical inefficiency g, is related to the autocorrelation time
    by g = 1+2*tau

    This method is based on the paper by J. D. Chodera [1], and the implementation at
    https://github.com/choderalab/pymbar - extending the code to support multidimensional data.

    References
    ----------
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
    histogram analysis method for the analysis of simulated and parallel tempering simulations.
    JCTC 3(1):26-41, 2007.

    Parameters
    ----------
    time_series: np.ndarray, shape=(num_frames, num_dimensions), dtype=float
        The time series to calculate the statistical inefficiency of.
    minimum_samples: int
        The minimum number of data points to consider in the calculation.

    Returns
    -------
    float:
        The statistical inefficiency.
    """

    number_of_timesteps = time_series.shape[0]
    time_series_dimension = 1 if len(time_series.shape) == 1 else time_series.shape[1]

    time_series_mean = time_series.mean(0)

    shifted_data = time_series.astype(np.float32) - time_series_mean

    sigma_squared_array = np.zeros(number_of_timesteps)

    for i in range(number_of_timesteps):

        if time_series_dimension > 1:
            sigma_squared_array[i] = shifted_data[i].dot(shifted_data[i])
        else:
            sigma_squared_array[i] = shifted_data[i] * shifted_data[i]

    sigma_squared = sigma_squared_array.mean()

    if sigma_squared == 0:
        raise ParameterError('Sample covariance sigma_AB^2 = 0 -- cannot compute statistical inefficiency')

    current_timestep = 1
    statistical_inefficiency = 1.0

    while current_timestep < number_of_timesteps - 1:

        autocorrelation_array = np.zeros([number_of_timesteps - current_timestep, time_series_dimension])

        for i in range(number_of_timesteps - current_timestep):

            if time_series_dimension > 1:
                autocorrelation_array[i] = shifted_data[i].dot(shifted_data[i + current_timestep])
            else:
                autocorrelation_array[i] = shifted_data[i] * shifted_data[i + current_timestep]

        autocorrelation_function = autocorrelation_array.mean() / sigma_squared

        if autocorrelation_function <= 0.0 and current_timestep > minimum_samples:
            break

        statistical_inefficiency += (2.0 * autocorrelation_function *
                                     (1.0 - float(current_timestep) / float(number_of_timesteps)))

        current_timestep += 1

    # Enforce a minimum autocorrelation time of 0.
    if statistical_inefficiency < 1.0:
        statistical_inefficiency = 1.0

    return statistical_inefficiency


def calculate_autocorrelation_time(time_series, minimum_samples=3):
    """Calculates the autocorrelation time of a time series via its
    statistical inefficiency.

    Parameters
    ----------
    time_series: np.ndarray, shape=(num_frames, num_dimensions), dtype=float
        The time series to calculate the autocorrelation time of.
    minimum_samples: int
        The minimum number of data points to consider in the calculation.

    Returns
    -------
    float:
        The autocorrelation time.
    """

    statistical_inefficiency = calculate_statistical_inefficiency(time_series, minimum_samples)
    return (statistical_inefficiency - 1.0) / 2.0


def detect_equilibration(time_series, minimum_samples=3):
    """Detect when a time series set has effectively become stationary (i.e has reached equilibrium).

    Notes
    -----
    This method is based on the paper by J. D. Chodera [1], and the implementation at
    https://github.com/choderalab/pymbar - extending the code to support multidimensional data.

    References
    ----------
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
        histogram analysis method for the analysis of simulated and parallel tempering simulations.
        JCTC 3(1):26-41, 2007.

    Parameters
    ----------
    time_series: np.ndarray, shape=(num_frames, num_dimensions), dtype=float
        The time series to analyse.
    minimum_samples: int
        The minimum number of data points to consider in the calculation.

    Returns
    -------
    float:
        The time at which the data has reached equilibrium.
    float:
        The statistical inefficiency of the data.
    int:
        The effective number of uncorrelated samples.
    """

    number_of_timesteps = time_series.shape[0]
    statistical_inefficiency_array = np.ones([number_of_timesteps - 1], np.float32)

    # Special case if the time series is constant.
    if time_series.std() == 0.0:
        return 0, 1, 1

    effect_samples_array = np.ones([number_of_timesteps - 1], np.float32)

    current_timestep = 0

    for current_timestep in range(0, number_of_timesteps - 1):

        try:
            statistical_inefficiency_array[current_timestep] = calculate_statistical_inefficiency(
                time_series[current_timestep:number_of_timesteps], minimum_samples)
        except ParameterError:  # Fix for issue https://github.com/choderalab/pymbar/issues/122
            statistical_inefficiency_array[current_timestep] = (number_of_timesteps - current_timestep + 1)

        effect_samples_array[current_timestep] = (number_of_timesteps - current_timestep + 1) / \
                                                 statistical_inefficiency_array[current_timestep]

    maximum_effective_samples = effect_samples_array.max()
    equilibration_time = effect_samples_array.argmax()
    statistical_inefficiency = statistical_inefficiency_array[current_timestep]

    return equilibration_time, statistical_inefficiency, maximum_effective_samples


def decorrelate_time_series(time_series):
    """Extracts an uncorrelated sub-time series from a possibly correlated one.

    Parameters
    ----------
    time_series : np.ndarray, shape=(num_frames, num_dimensions), dtype=float
        The possibly correlated time series.

    Returns
    -------
    np.ndarray, shape=(num_frames, num_dimensions), dtype=float
        The uncorrelated time series from which the average was calculated.
    int
        The index after which the data is considered well equilibrated.
    float
        The statistical inefficiency of the original time series.
    """

    # Compute the indices of the uncorrelated time series
    [equilibration_index, inefficiency, effective_samples] = detect_equilibration(time_series)
    equilibrated_data = time_series[equilibration_index:]

    # Extract a set of uncorrelated data points.
    indices = get_uncorrelated_indices(equilibrated_data.shape[0], inefficiency)
    uncorrelated_time_series = equilibrated_data[indices]

    return uncorrelated_time_series, equilibration_index, inefficiency


def get_uncorrelated_indices(time_series_length, statistical_inefficiency):
    """Returns the indices of the uncorrelated frames of a time series.

        Parameters
        ----------
        time_series_length : int
            The length of the time series to extract frames from.
        statistical_inefficiency: float
            The statistical inefficiency of the time series.

        Returns
        -------
        list of int
            The indices of the uncorrelated frames.
        """

    # Extract a set of uncorrelated data points.
    stride = int(math.ceil(statistical_inefficiency))
    return range(0, time_series_length, stride)
