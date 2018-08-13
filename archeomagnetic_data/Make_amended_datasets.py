# Makes other datasets that are simple amendments of Paris-700
import numpy as np

Paris = np.loadtxt('Paris700.txt')

Paris[:,3] = 0
np.savetxt('Paris700_no_age_errors.txt',Paris,fmt='%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %i')

New_Paris = Paris.copy()
New_Paris[:,5] = New_Paris[:,5] * 2
np.savetxt('Paris700_no_age_errors_twice_F_error.txt',New_Paris,fmt='%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %i')


for i in range(0,len(Paris[:,5])):
    Paris[i,5] = max(5,Paris[i,5])
np.savetxt('Paris700_no_age_errors_min_5muT.txt',Paris,fmt='%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %i')
