#  Collection of fit models
#
#

import math
import numpy as np

#
#  Staggered fit model for two states
#
def stagg_2_state(t, a0, m0, a1, m1):
     ''' a0*np.exp(-m0*t)  + (-1)**t*a1*np.exp(-m1*t)   '''
     ss =(-1)**t
     ans = a0*np.exp(-m0*t)  + ss*a1*np.exp(-m1*t)  
     return ans
