import numpy as np
import math

##
##  Common jackknie routines
##


def jackknife(jloop,nsample):
  jmean = 0.0
  for isample in  range(0,nsample) :
    jmean += jloop[isample]

  jmean  /= nsample

  jtot = 0.0
  for isample in  range(0,nsample) :
    ttt   = (jloop[isample] - jmean )
    jtot +=  ttt * ttt

  return math.sqrt( (nsample - 1.0)/(1.0*nsample) * jtot     )



def load_data(filename_in, nt,no_config, tag, scale = 1.0 ):

  corr = np.zeros(( nt,no_config ))

  print ("Loading correlators from " , filename_in )
  try:
    f = open(filename_in, "r")
  except:
    print ("Error opening " , filename_in)
    sys.exit(1)

  iconfig = 0 
  for line in f:
    ll  = line.rstrip('\n')
    tmp = ll.split()
    tag_l = tmp.pop(0)
    if tag_l == tag :
       t = 0 
       for xx in tmp:
#          print (xx)
          corr[t,iconfig] = scale*float(xx)
          t = t + 1
       iconfig = iconfig + 1
#
  return corr
#


def calc_meff(corr, nt,no_config, ss = 0, ainv = 1) :

    tmp  = np.zeros(( no_config ))

    tmax = int(nt / 2)  - 1
    tt     = np.zeros(( tmax ))
    mm     = np.zeros(( tmax ))
    mm_err = np.zeros(( tmax ))
  
    for t in range(0,tmax):
       tt[t] = t + ss
       now = 0.0
       inc = 0.0 
       for i in range(0, no_config  ):
          ok = True

          now = now + corr[t,i] 
          inc = inc + corr[t+1,i] 

          #  jackknife analysis
          jnow = 0
          jinc = 0
          for j in range(0, no_config  ):
            if i != j :
              jnow = jnow + corr[t  ,j] 
              jinc = jinc + corr[t+1,j] 

          jnow = jnow / (no_config - 1) 
          jinc = jinc / (no_config - 1) 
          if jnow > 0 and jinc > 0 :
             tmp[i] = math.log(jnow / jinc)
          else:
             ok = False

       now = now / no_config
       inc = inc / no_config

       if (now > 0 and inc > 0) and ok   :
          meff = ainv * math.log(now/inc)
          jerr = ainv * jackknife(tmp,no_config)
       else:
          meff = -1
          jerr = 0.0

       mm[t]     = meff
       mm_err[t] = jerr

#       print(t, meff, jerr)


    return tt, mm, mm_err




def calc_corr(corr, nt,no_config, ss = 0) :

    tmp  = np.zeros(( no_config ))

    tmax =   nt 
    tt     = np.zeros(( tmax ))
    mm     = np.zeros(( tmax ))
    mm_err = np.zeros(( tmax ))
  
    for t in range(0,tmax):
       tt[t] = t + ss
       now = 0.0
       for i in range(0, no_config  ):
          ok = True

          now = now + corr[t,i] 
          #  jackknife analysis
          jnow = 0
          for j in range(0, no_config  ):
            if i != j :
              jnow = jnow + corr[t  ,j] 

          tmp[i] = jnow / (no_config - 1) 

       mm[t]     = now / no_config
       mm_err[t] = jackknife(tmp,no_config)

#       print(t, meff, jerr)


    return tt, mm, mm_err
