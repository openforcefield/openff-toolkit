import os
import numpy
from numpy import array
import commands as c
import re
from smarty.score_utils import *

iterations = 5000
times = 3
temperatures = [0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1]

output = open('testing_summary.txt', 'w')
output.writelines("%5s %15s %15s %15s\n" % ('run', "temperature", "Max Score", "Final Score"))
sep = "="*50
#for frag in ['VdW', 'Bond', 'Angle', 'Torsion']:
for frag in ['VdW', 'Angle', 'Torsion']:
    print frag
    os.chdir('%s/' % frag)
    output.writelines("%s\n%s\n\n" % (sep, frag))
    for temp in temperatures:
        for idx in range(times):
            out = "%s_%i_%.2e" % (frag, idx, temp)
            command = "./run_smirky.bash %i %f %s" % (iterations, temp, out)
            print c.getoutput(command)
            timeseries = load_trajectory("%s.csv" % out)
            time_fractions = scores_vs_time(timeseries, 'fractionmatched')
            max_score = "%.3f" % max(time_fractions['all'])
            final = "%.3f" % time_fractions['all'][-1]
            write = "%5s %15s %15s %15s\n" % (idx, str(temp), max_score, final)
            print write
            output.writelines(write)
    os.chdir('../')

output.close()
