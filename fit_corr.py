#!/usr/local/bin/python3


#  parameters of the fit
t_start = 2
t_end   = 10
corr_tag = "onemp.ll"

##ainv_gev = 0.1973 / 0.15 
ainv_gev = 1.0

##
##
##

import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize, report_fit

# my modules
import sys
sys.path.append('src')

import util
import models

fff = "t0_onempHy_m0.450.txt"

fit_model = models.stagg_2_state

nt = 96
no_config = 1008




corr = util.load_data(fff, nt,no_config, corr_tag  )
corr *= -1 

print("Computing jackknife correlators")
tt, corr_mean, corr_err = util.calc_corr(corr, nt,no_config, 0.0) 


plt.errorbar(tt, corr_mean , corr_err  ,  fmt= 'ro' )


##
##  setup plots
##

plt.title("Hybrid meson correlator at charm")
##plt.legend()

plt.yscale('log')

plt.xlabel('t')
plt.ylabel('corr(t)')
plt.xlim(0,10.5)
##plt.ylim(0,8)

##plt.text(0.3,1.4 , "Ensemble 1648f211b580m013m065m838a (200 configs)")
##plt.text(0.3,1.2 , "Charm mass m = 0.838")

###plt.savefig("onemp_corr.pdf")



print("Starting the fits")




##def model_1(t, a0, m0, a1, m1):
##     ss =(-1)**t
##     ans = a0*np.exp(-m0*t)  + ss*a1*np.exp(-m1*t)  
##     return ans




#
#  err_global**2 is the chi**2/
#
def err_global(p, x, data, data_err):
    # p is now a_1, b, c_1, a_2, c_2, with b shared between the two
    a0 = p['a0']
    a1 = p['a1']

    m0 = p['m0']
    m1 = p['m1']

##    err1 = (data - model_1(x, a0, m0, a1, m1) ) / data_err
    err1 = (data - fit_model(x, a0, m0, a1, m1) ) / data_err

    return err1


fit_params = Parameters()
fit_params.add('a0', value=1)
fit_params.add('a1', value=1)

fit_params.add('m0', value=3)
fit_params.add('m1', value=3)



out = minimize(err_global, fit_params, args=(tt[t_start:t_end], corr_mean[t_start:t_end], corr_err[t_start:t_end] ))


##report_fit(out.params)



# 
print("Channel " , corr_tag)
print("Time fit region " , t_start , t_end )
print("Fit model: " , fit_model.__doc__)

##print("chi**2 =  {:.2e} " . format( out.chisqr ))
print("chi**2/dof =  {:.2e} " .format(out.redchi ))

a0best = out.params["a0"].value
m0best = out.params["m0"].value

a1best = out.params["a1"].value
m1best = out.params["m1"].value

print("Mass, m0 = {:.2e}" .format(m0best))
print("Mass, m1 = {:.2e}" .format(m1best))

print("Amplitude, a0 = {:.2e}" .format(a0best) )
print("Amplitude, a1 = {:.2e}" .format( a1best) )


##
##   plot the fit model of the correlator
##
y_fit = [] 
t_fit = [] 

for t,c_,ce_ in zip(tt[t_start:t_end],corr_mean[t_start:t_end] , corr_err[t_start:t_end]  ) :
   fit_v =  fit_model(t, a0best, m0best, a1best, m1best)
#   print (t, "fit=" , fit_v, "data ", c_, ce_)
   t_fit.append(t)
   y_fit.append(fit_v)


plt.plot(t_fit, y_fit, "g") 

plt.show()

