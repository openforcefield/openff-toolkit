# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Dask Property Estimator Backend

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import multiprocessing

from dask import distributed

from openforcefield.properties.estimator.backends.base import PropertyEstimatorBackend, BackendResources


# =============================================================================================
# Base Backend Definition
# =============================================================================================

class DaskLocalClusterBackend(PropertyEstimatorBackend):
    """A property estimator backend which uses a dask `LocalCluster` to
    run calculations.
    """

    def __init__(self, number_of_workers=1, threads_per_worker=None, resources_per_task=BackendResources()):
        """Constructs a new DaskLocalClusterBackend"""

        super().__init__(number_of_workers, threads_per_worker, resources_per_task)

        maximum_threads = multiprocessing.cpu_count()
        requested_threads = number_of_workers * threads_per_worker * resources_per_task.number_of_threads

        if requested_threads > maximum_threads:

            raise ValueError('The total number of requested threads ({})is greater than is available on the'
                             'machine ({})'.format(requested_threads, maximum_threads))

        # TODO: Check GPUs

        self._cluster = None
        self._client = None

    def start(self):

        self._cluster = distributed.LocalCluster(self._number_of_workers,
                                                 self._threads_per_worker,
                                                 processes=False,
                                                 resources=self._resources_per_task.dict())

        self._client = distributed.Client(self._cluster,
                                          processes=False)

    def stop(self):

        self._client.close()
        self._cluster.close()

    def submit_task(self, function, *args):

        return self._client.submit(function, *args, resources=self._resources_per_task.dict(),
                                   available_resources=self._resources_per_task)
