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


# =============================================================================================
# Base Backend Definition
# =============================================================================================

class BackendResources:
    """An object which stores how many of each type of computational resource
    (threads or gpu's) is available to a calculation task."""

    @property
    def number_of_threads(self):
        return self._number_of_threads

    @property
    def number_of_gpus(self):
        return self._number_of_gpus

    def __init__(self, number_of_threads=1, number_of_gpus=0):
        """Constructs a new BackendResources object.

        Parameters
        ----------
        number_of_threads: int
            The number of the threads available.
        number_of_gpus
            The number of the gpu's available.
        """

        self._number_of_threads = number_of_threads
        self._number_of_gpus = number_of_gpus
        
        assert self._number_of_threads >= 0
        assert self._number_of_gpus >= 0

        assert self._number_of_threads > 0 or self._number_of_gpus > 0

    def dict(self):
        return self.__getstate__()

    def __getstate__(self):
        return {
            'number_of_threads': self.number_of_threads,
            'number_of_gpus': self.number_of_gpus
        }

    def __setstate__(self, state):

        self._number_of_threads = state['number_of_threads']
        self._number_of_gpus = state['number_of_gpus']

    def __eq__(self, other):
        return self.number_of_threads == other.number_of_threads and \
               self.number_of_gpus == other.number_of_gpus

    def __ne__(self, other):
        return not self.__eq__(other)


class PropertyEstimatorBackend:
    """An abstract base representation of a property estimator backend.

    A backend will be responsible for coordinating and running calculations
    on the available hardware.

    Notes
    -----
    All estimator backend classes must inherit from this class, and must implement the
    `start`, `stop`, and `submit_task` method.
    """

    def __init__(self, number_of_workers=1, threads_per_worker=None, resources_per_task=BackendResources()):
        """Constructs a new PropertyEstimatorBackend object.

        Parameters
        ----------
        number_of_workers : int
            The number of works to run the calculations on. One worker
            can perform a single task (e.g run a simulation) at once.
        threads_per_worker : int, optional
            The number of threads per each launched worker.
        resources_per_task: BackendResources
            The number of resources available to each calculation task.
        """

        self._number_of_workers = number_of_workers
        self._threads_per_worker = threads_per_worker

        self._resources_per_task = resources_per_task

    def start(self):
        """Start the calculation backend."""
        pass

    def stop(self):
        """Stop the calculation backend."""
        pass

    def submit_task(self, function, *args, **kwargs):
        """Submit a task to the compute resources
        managed by this backend.

        Parameters
        ----------
        function: function
            The function to run.

        Returns
        -------
        Future
            Returns a future object which will eventually point to the results
            of the submitted task.
        """
        pass
