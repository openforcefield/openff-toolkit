# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Base Property Estimator Backend

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import math
import logging
import multiprocessing
import os


# =============================================================================================
# Base Backend Definition
# =============================================================================================

class PropertyEstimatorBackend:

    def __init__(self, number_of_workers=1, threads_per_worker=None):
        """Constructs a new PropertyEstimatorBackend object.

        Parameters
        ----------
        number_of_workers : int
            The number of works to run the calculations on. One worker
            can perform a single task (e.g run a simulation) at once.
        threads_per_worker : int, optional
            The number of threads available to each worker. This enables
            multi-threaded tasks (e.g run a simulation on more than one core).

            If None, the workers will split any unused threads between
            themselves.
        """

        self._number_of_workers = number_of_workers
        self._threads_per_worker = threads_per_worker

        self._calculate_number_of_simulation_threads()

    def _calculate_number_of_simulation_threads(self):
        """Determines how many threads will be used per simulation
        if no value is given by the user. The default option is
        to use as many threads as are available.
        """
        maximum_threads = multiprocessing.cpu_count()

        if self._threads_per_worker is None:
            self._threads_per_worker = math.floor(maximum_threads / self._number_of_workers)

        total_threads = self._number_of_workers * self._threads_per_worker

        if total_threads > maximum_threads:

            raise ValueError('The total number of requested threads ({}) must be less '
                             'than the available {}.'.format(total_threads, maximum_threads))

        logging.info(str(self._threads_per_worker) + ' threads will be used per worker')

        os.environ["OPENMM_NUM_THREADS"] = str(self._threads_per_worker)

    def start(self):
        """TODO: Missing docstrings."""
        pass

    def submit_task(self, function, *args, **kwargs):
        """Submit a task to the compute resources
        managed by this backend.

        Parameters
        ----------
        function: function
            The function to run.
        """
        pass
