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

from openforcefield.properties.estimator.backends.base import PropertyEstimatorBackend

from dask import distributed


# =============================================================================================
# Base Backend Definition
# =============================================================================================

class DaskLocalClusterBackend(PropertyEstimatorBackend):

    def __init__(self, number_of_workers=1, threads_per_worker=None):
        """Constructs a new DaskLocalClusterBackend"""
        super().__init__(number_of_workers, threads_per_worker)

        self._cluster = distributed.LocalCluster(number_of_workers, 1, processes=False)
        self._client = distributed.Client(self._cluster, processes=False)

    def submit_task(self, function, *args):
        return self._client.submit(function, *args)
