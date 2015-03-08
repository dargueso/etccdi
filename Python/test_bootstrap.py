#!/usr/bin/env python

""" test_boostrap.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Thu Feb 12 16:56:19 AEDT 2015

"""
import netCDF4 as nc
import numpy as np
import datetime as dt
import calendar

import pdb

def calc_thres(varin,nyears):
    
    varinp=np.ones((365,),dtype=np.float64)
    doy_s=np.tile(np.arange(365),nyears)
    use_doy=np.zeros(doy_s.shape,dtype=np.int)
    for d in range(0,365):
      use_doy[:]=1
      init_d=(d-2)%365
      end_d=(d+2)%365
      if (init_d>end_d):
        use_doy[(doy_s<init_d) & (doy_s>=end_d)]=0
      else:
        use_doy[(doy_s<init_d) | (doy_s>=end_d)]=0
      varinp[d]=np.asarray(np.percentile(varin[use_doy==1],90))
      
    return varinp

i=765 #longitude
j=193 #latitude
fin=nc.Dataset("tmax.1961-1990.nc","r")

tmax=fin.variables['tmax'][:,j,i]



dates = nc.num2date(fin.variables['time'][:],units=fin.variables['time'].units,calendar='standard')
years_all = np.asarray([dates[i].year for i in xrange(len(dates))])
months_all= np.asarray([dates[i].month for i in xrange(len(dates))])
days_all=np.asarray([dates[i].day for i in xrange(len(dates))])

#Removing leap days
dates=dates[((months_all==2) & (days_all==29))==False]
years=np.asarray([dates[i].year for i in xrange(len(dates))])
months=np.asarray([dates[i].month for i in xrange(len(dates))])  
days=np.asarray([dates[i].day for i in xrange(len(dates))])

tmax=tmax[((months_all==2) & (days_all==29))==False,:,:]

nyears=years[-1]-years[0]+1
bsyear=years[0]
beyear=years[-1]

tmaxpbs=np.ones((nyears,nyears-1,365),dtype=np.float64)
tmax_nobs=np.ones((365,),dtype=np.float64)
tmax_29yr=np.ones((365,nyears),dtype=np.float64)


print nyears
for yr in range(nyears):
  tmax_boot=np.resize(np.repeat(tmax,29),(len(tmax),29))
  
  yr_n=0
  print "year: %s" %(yr)
  for yr_iter in range(nyears):
    print yr_iter+bsyear
    if yr!=yr_iter:
      tmax_boot[years==yr+bsyear,yr_n]=tmax[years==yr_iter+bsyear]
      tmaxpbs[yr,yr_n,:]=calc_thres(tmax_boot[:,yr_n],nyears)
      yr_n=yr_n+1
      

  tmax_29yr[:,yr]=calc_thres(tmax[years!=yr+bsyear],nyears-1)
tmax_nobs=calc_thres(np.mean(tmax_boot,axis=1),nyears)
  

TXp=np.zeros((nyears,nyears-1),dtype=np.int)
TXp_nobs=np.zeros((nyears,),dtype=np.int)
TXp_29yr=np.zeros((nyears,),dtype=np.int)
for yr in range(nyears):
  yr_n=0
  for yr_iter in range(nyears):
    if yr_iter!=yr:
      TXp[yr,yr_n]=np.sum((tmax[(years==yr+bsyear)]>tmaxpbs[yr,yr_n,:]))
      yr_n=yr_n+1
      
  TXp_nobs[yr]=np.sum((tmax[(years==yr+bsyear)]>tmax_nobs[:]))
  TXp_29yr[yr]=np.sum((tmax[(years==yr+bsyear)]>tmax_29yr[:,yr]))
sys.exit()



import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(range(30),TXp,color='0.9',linestyle='--')
ax.plot(range(30),np.mean(TXp,axis=1),color='0.6')
ax.plot(range(30),TXp_nobs,color='k')
ax.plot(range(30),TXp_29yr,color='r')
plt.show
