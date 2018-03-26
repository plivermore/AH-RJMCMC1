import numpy as np
import os
import sys

""" Script to remove outlying data
    Requires three files to be in the current directory, all outputs from the RJ-MCMC codes:
    1) data.dat
    2) credible_lower.dat
    3) credible_upper.dat

    The output is to the file:   data_within_credible_bounds.dat
"""

if not ( os.path.exists('credible_lower.dat') and os.path.exists('credible_upper.dat') and os.path.exists('data.dat') ):
    print('One or more of the required files are missing from the current directory')
    sys.exit(0)

lx, ly = np.loadtxt('credible_lower.dat', unpack=True)
ux, uy = np.loadtxt('credible_upper.dat', unpack=True)
x, x_err, y, y_err, strat = np.loadtxt('data.dat', unpack=True)

count = 0
f = open('data_within_credible_bounds.dat','w')
for i in range(len(x)):
    if y[i] >= np.interp(x[i],lx,ly) and y[i] <= np.interp(x[i],ux,uy):
        f.writelines("%10.3f %10.3f %10.3f %10.3f %i \n" % (x[i], x_err[i], y[i], y_err[i], strat[i] ) )
        count += 1

f.close()
print('Of the ' + str(len(x)) + ' data read, ' + str(count) + ' data items fell within the credible interval')
