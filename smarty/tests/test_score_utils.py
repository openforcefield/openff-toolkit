import unittest
import smarty
from smarty import score_utils
import sys

class TestScoreUtils(unittest.TestCase):
    def test_processing_trajfile(self):
        """
        Test ability to load and process trajectory files from smarty test
        """
        # load_trajectory file: converts csv to dictionary
        timeseries = score_utils.load_trajectory('test_smirky.csv')
        # scores_vs_time: converts numerator/denominator to fractional scores
        scores_vs_time = score_utils.scores_vs_time(timeseries)
        # Test plotting function
        score_utils.create_plot_file('test_smirky.csv', 'test_smirky.pdf')

