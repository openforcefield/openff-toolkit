import os
import numpy
from numpy import array
import commands as c
import re
from smarty.score_utils import *

iterations = 5000
times = 3
temperatures = [0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1]

output = open('testing_summary_2.txt', 'w')
header = "%5s %10s %10s %10s %10s\n" % ('run', "temperature", "Max Score", "Final Score", "First 100%")
sep = "="*70
#for frag in ['VdW', 'Bond', 'Angle', 'Torsion']:
for frag in ['VdW','Bond','Angle']:
    os.chdir('%s/' % frag)
    output.writelines("%s\n%s\n%s\n" % (sep, frag, header))
    for temp in temperatures:
        for idx in range(times):
            out = "%s_%i_%.2e" % (frag, idx, temp)
            timeseries = load_trajectory("%s.csv" % out)
            time_fractions = scores_vs_time(timeseries, 'fractionmatched')
            max_score = max(time_fractions['all'])
            if max_score < 1.0:
                iteration = str(-1)
            else:
                iteration = str(list(time_fractions['all']).index(1.0))
            max_score = "%.3f" % max_score
            final = "%.3f" % time_fractions['all'][-1]
            write = "%5s %10s %10s %10s %10s\n" % (idx, str(temp), max_score, final, iteration)
            print write
            output.writelines(write)
    os.chdir('../')

output.close()
