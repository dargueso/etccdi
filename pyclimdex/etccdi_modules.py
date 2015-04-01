#!/usr/bin/env python

""" etccdi_modules.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Thu Feb  5 16:30:24 AEDT 2015

"""
import re
import sys
import os
import glob
import calendar
from constants import const
import varinfo as vinf
import datetime as dt
import numpy as np
import netCDF4 as nc

import pdb

varinfo=vinf.VariablesInfo()


def read_input(filename):
  """
  Namelist with input arguments
  """
  
  filein = open(filename,"r")
  
  lines = filein.readlines()
  
  inputinf = {}
  entryname = []
  entryvalue = []
  
  for line in lines:
    line=re.sub('\s+',' ',line)
    line=line.replace(" ", "")
    line=line.replace('"', '')
    li=line.strip()
    #Ignore empty lines
    if li:
      #Ignore commented lines
      if not li.startswith("#"): 
        values=li.split('=')
        entryname.append(values[0])
        entryvalue.append(values[1])
        
  filein.close()
  
  for ii in xrange(len(entryname)):
    inputinf[entryname[ii]]=entryvalue[ii]
    
  return inputinf

###############################################
###############################################

  
def calc_dates(syear,eyear):
  """
  Function to calculate vectors with dates, years, months
  """
  d1 = dt.datetime(syear,1,1,0,0)
  d2 = dt.datetime(eyear,12,31,0,0)
  diff = d2 - d1
  dates=np.asarray([d1 + dt.timedelta(i) for i in range(diff.days+1)])
  years = np.asarray([dates[i].year for i in xrange(len(dates))])
  months = np.asarray([dates[i].month for i in xrange(len(dates))])

  return dates,years,months

###############################################
###############################################

def calc_qualitymask(var,years,inputinf):
  """
  Function to mask out years without enough data according to NoMissingThreshold
  """
  missthres=float(inputinf['NoMissingThreshold'])
  
  nyears=years[-1]-years[0]+1
  quality_mask=np.zeros((nyears,)+var.shape[1:],dtype=np.bool)
  for yr in range(nyears):
    year=years[0]+yr
    quality_mask[yr,:,:]=(np.sum(var[years==year].mask,axis=0)/float(np.sum(years==year)))>(1.-missthres)
  return quality_mask


###############################################
###############################################

def calc_otime(years,mode):
	"""
	Function to calculate the time vector to write the netCDF
	"""
	
	if mode == 'years':
		otime = np.asarray([dt.datetime(x,6,1,0,0) for x in xrange(years[0],years[-1]+1)])
	if mode == 'months':
		otime=[]
		for yr in range(years[0],years[-1]+1):
			for months in range(1,13):
				otime.append(dt.datetime(yr,months,15,0,0))
		otime=np.asarray(otime)
	
	return otime
				
###############################################
###############################################

  
def calc_FD(tmax,tmin,years):
  """
  Function to calculate extremes defined using thresholds (FD,SU,ID,TN)
  
  """
  nyears=years[-1]-years[0]+1
  FD=np.ones((nyears,)+tmax.shape[1:],dtype=np.float)*const.missingval
  ID=np.ones((nyears,)+tmax.shape[1:],dtype=np.float)*const.missingval
  SU=np.ones((nyears,)+tmax.shape[1:],dtype=np.float)*const.missingval
  TR=np.ones((nyears,)+tmax.shape[1:],dtype=np.float)*const.missingval
  
  print "Procesing FD, ID, SU, TR..."
  
  for yr in range(nyears):
    year=years[0]+yr
    #Frosting days
    FD[yr,:,:]=np.ma.sum(tmin[years==year,:,:]<0.,axis=0)
    #Icing days
    ID[yr,:,:]=np.ma.sum(tmax[years==year,:,:]<0.,axis=0)
    #Summer days
    SU[yr,:,:]=np.ma.sum(tmax[years==year,:,:]>25.,axis=0)
    #Tropical nights
    TR[yr,:,:]=np.ma.sum(tmin[years==year,:,:]>20.,axis=0)
  
       
  
  return FD,ID,SU,TR
  

def calc_gsl(tmax,tmin,years,months):
    """
    Function to calculate extremes defined using thresholds (FD,SU,ID,TN)

    """
###############################################
###############################################

def calc_TXx(tmax,tmin,years,months):
  """
   Function to calculate absolute extremes (TXx,TXn,TNn,TNx,DTR)
   Extremes are caculated at monthly and annual timescales
  """
  print "Processing TXx,TXn,TNn,TNx,DTR ..."
  
  nyears=years[-1]-years[0]+1
  TXx = np.ones((nyears*12,)+tmax.shape[1:],dtype=np.float)*const.missingval
  TXn = np.ones((nyears*12,)+tmax.shape[1:],dtype=np.float)*const.missingval
  TNn = np.ones((nyears*12,)+tmax.shape[1:],dtype=np.float)*const.missingval
  TNx = np.ones((nyears*12,)+tmax.shape[1:],dtype=np.float)*const.missingval
  DTR = np.ones((nyears*12,)+tmax.shape[1:],dtype=np.float)*const.missingval
  for yr in range(nyears):
    year=years[0]+yr
    for month in range(1,13):
      index_ym=yr*12+month-1
      TXx[index_ym,:,:] = np.ma.max(tmax[(years==year) & (months==month),:,:],axis=0)
      TNx[index_ym,:,:] = np.ma.max(tmin[(years==year) & (months==month),:,:],axis=0)
      TXn[index_ym,:,:] = np.ma.min(tmax[(years==year) & (months==month),:,:],axis=0)
      TNn[index_ym,:,:] = np.ma.min(tmin[(years==year) & (months==month),:,:],axis=0)
      DTR[index_ym,:,:] = np.ma.mean(tmax[(years==year) & (months==month),:,:]-tmin[(years==year) & (months==month),:,:],axis=0)
      

  return TXx,TXn,TNn,TNx,DTR

###############################################
###############################################
  
def calc_R10mm(prec,years,Rnnmm_t):
  """
   Function to calculate threshold extremes for precipitation (R10mm, R20mm, Rnnmm,SDII)

  """
  
  print "Processing R10mm, R20mm, Rnnmm,SDII ..."
  
  nyears=years[-1]-years[0]+1
  R10mm=np.ones((nyears,)+prec.shape[1:],dtype=np.float)*const.missingval
  R20mm=np.ones((nyears,)+prec.shape[1:],dtype=np.float)*const.missingval
  Rnnmm=np.ones((nyears,)+prec.shape[1:],dtype=np.float)*const.missingval
  SDII=np.ones((nyears,)+prec.shape[1:],dtype=np.float)*const.missingval
  
  for yr in range(nyears):
    year=years[0]+yr
    
    #R10mm
    R10mm[yr,:,:]=np.ma.sum(prec[years==year,:,:]>10.,axis=0)
    #R20mm
    R20mm[yr,:,:]=np.ma.sum(prec[years==year,:,:]>20.,axis=0)
    #Rnnmm
    Rnnmm[yr,:,:]=np.ma.sum(prec[years==year,:,:]>Rnnmm_t,axis=0)
    #SDII
    prec1mm=np.ma.sum(prec[prec[years==year,:,:]>1.],axis=0)
    days1mm=np.ma.sum(prec[years==year,:,:]>1.,axis=0)
    SDII[yr,:,:]=np.ma.sum(prec1mm,axis=0)/days1mm.astype('float')
    
    
    
  return R10mm, R20mm, Rnnmm,SDII

###############################################
###############################################

def calc_Rx5day(prec,years,months):
  """
   Function to calculate absolute extremes (Rx1day,Rx5day)
   Extremes are caculated at monthly and annual timescales
  """
  print "Processing Rx1day,Rx5day ..."
  
  nyears=years[-1]-years[0]+1
  Rx1day = np.ones((nyears*12,)+prec.shape[1:],dtype=np.float)*const.missingval
  Rx5day = np.ones((nyears*12,)+prec.shape[1:],dtype=np.float)*const.missingval

  prec5day=prec.copy()
  prec5day[4:,:,:]=prec[4:,:,:]+prec[3:-1,:,:]+prec[2:-2,:,:]+prec[1:-3,:,:]+prec[0:-4,:,:]
  for yr in range(nyears):
    year=years[0]+yr
    for month in range(1,13):
      index_ym=yr*12+month-1
      Rx5day[index_ym,:,:] = np.ma.max(prec5day[(years==year) & (months==month),:,:],axis=0)
      Rx1day[index_ym,:,:] = np.ma.max(prec[(years==year) & (months==month),:,:],axis=0)
      
  return Rx1day,Rx5day
  
###############################################
###############################################


def calc_R95p(prec,years,inputinf):
  """
   Function to calculate percntile extremes for precipitation (R95p, R99p, PRCPtot)

  """
  
  is_thresfile=int(inputinf['is_thresfile'])
  syear = int(inputinf['syear'])
  eyear = int(inputinf['eyear'])
  bsyear = int(inputinf['basesyear'])
  beyear = int(inputinf['baseeyear'])
  byrs = beyear-bsyear+1
  print "Processing R95p, R99p, PRCPtot ..."
  
  nyears=years[-1]-years[0]+1
  R95p=np.ones((nyears,)+prec.shape[1:],dtype=np.float)*const.missingval
  R99p=np.ones((nyears,)+prec.shape[1:],dtype=np.float)*const.missingval
  PRCPtot=np.ones((nyears,)+prec.shape[1:],dtype=np.float)*const.missingval
  prec95=np.ones(prec.shape[1:],dtype=np.float)*const.missingval
  prec99=np.ones(prec.shape[1:],dtype=np.float)*const.missingval
  
  if is_thresfile==0:
    print "No threshold file required. Percentiles will be calculated (precip)"
    prec_base=prec[(years>=bsyear) & (years<=beyear),:,:]
    for i in range(prec.shape[1]):
      for j in range(prec.shape[2]):
        prec_base=prec[(years>=bsyear) & (years<=beyear),i,j]
        prec_base_rain=np.ma.masked_less_equal(prec_base,1.)
        if np.count_nonzero(~prec_base_rain.mask)!=0:
          prec95[i,j]=np.percentile(prec_base_rain[~prec_base_rain.mask],95,axis=0)
          prec99[i,j]=np.percentile(prec_base_rain[~prec_base_rain.mask],99,axis=0)
    
    
  elif is_thresfile==1:
    print "Threshold file will be used. No percentiles will be calculated (precip)"
    thresfile=nc.Dataset(inputinf['thres_filename'],'r')
    prec95=thresfile.variables['prec95'][:]
    prec99=thresfile.variables['prec95'][:]
    thresfile.close()
    
    
  else:
    sys.exit("ERROR:Wrong entry for is_thresfile (0/1)")
  
  for yr in range(nyears):
    year=years[0]+yr

    for i in range(prec.shape[1]):
      for j in range(prec.shape[2]):
        aux=prec[years==year,i,j]
        aux95=aux>prec95[i,j]
        aux99=aux>prec99[i,j]
        R95p[yr,i,j] = np.ma.sum(aux[aux95],axis=0)
        R99p[yr,i,j] = np.ma.sum(aux[aux99],axis=0)
        PRCPtot[yr,i,j] = np.ma.sum(aux[aux>1.],axis=0)
  
  
    
  return R95p,R99p,PRCPtot,prec95,prec99
###############################################
###############################################
def calc_TX10p(tmax,tmin,dates,inputinf):
  """
   Function to calculate percntile extremes for temperature (TN10p,TN50p,TN90p,TX10p,TX50p,TX90p)

  """
  is_thresfile=int(inputinf['is_thresfile'])
  syear = int(inputinf['syear'])
  eyear = int(inputinf['eyear'])
  bsyear = int(inputinf['basesyear'])
  beyear = int(inputinf['baseeyear'])
  byrs = beyear-bsyear+1
  version = inputinf['thres_version']
  print version
  
  print "Processing TN10p,TN50p,TN90p,TX10p,TX50p,TX90p ..."
  years_all=np.asarray([dates[i].year for i in xrange(len(dates))])
  months_all=np.asarray([dates[i].month for i in xrange(len(dates))])
  days_all=np.asarray([dates[i].day for i in xrange(len(dates))])
  
  #Removing leap days
  dates=dates[((months_all==2) & (days_all==29))==False]
  years=np.asarray([dates[i].year for i in xrange(len(dates))])
  months=np.asarray([dates[i].month for i in xrange(len(dates))])  
  days=np.asarray([dates[i].day for i in xrange(len(dates))])
  
  months_clim=months[:365]
  
  tmax=tmax[((months_all==2) & (days_all==29))==False,:,:]
  tmin=tmin[((months_all==2) & (days_all==29))==False,:,:]
  
  tmax_base=tmax[(years>=bsyear) & (years<=beyear),:,:]
  tmin_base=tmin[(years>=bsyear) & (years<=beyear),:,:]
  years_base=years[(years>=bsyear) & (years<=beyear)]
  
  nyears=years[-1]-years[0]+1
  
  
  TXp = np.zeros((3,nyears*12,)+tmax.shape[1:],dtype=np.float)
  TNp = np.zeros((3,nyears*12,)+tmax.shape[1:],dtype=np.float)
  tminp=np.ones((365,3)+tmin.shape[1:],dtype=np.float)*const.missingval
  tmaxp=np.ones((365,3)+tmin.shape[1:],dtype=np.float)*const.missingval
  tminpbs=np.ones((byrs,365,3)+tmin.shape[1:],dtype=np.float)*const.missingval
  tmaxpbs=np.ones((byrs,365,3)+tmin.shape[1:],dtype=np.float)*const.missingval
 
  

  
  if is_thresfile==0:
    if version=='bootstrap':
      print "No threshold file required. Percentiles will be calculated (temp)"
      print "A bootstrap will be carried out. This part take most of the time"
    
    
      tminp=calc_thres(tmin_base,byrs)
      tmaxp=calc_thres(tmax_base,byrs)
    
      for yr in range(byrs):
        print "year: %s" %(yr+bsyear)
        tmax_boot=tmax_base.copy()
        tmin_boot=tmin_base.copy()
        for yr_iter in range(byrs):
          if (yr_iter!=yr):
            tmax_boot[years_base==yr+bsyear,:,:]=tmax_boot[years_base==yr_iter+bsyear,:,:]
            tmin_boot[years_base==yr+bsyear,:,:]=tmin_boot[years_base==yr_iter+bsyear,:,:]
      
          tminpbs[yr,:,:,:,:]=calc_thres(tmin_boot,byrs)
          tmaxpbs[yr,:,:,:,:]=calc_thres(tmax_boot,byrs)
    
    elif version=='all_years':
      print "No threshold file required. Percentiles will be calculated (temp)"
      print "All years wihtin the base period will be used for all years (no boostrap or other method selected)"
      
      tminp=calc_thres(tmin_base,byrs)
      tmaxp=calc_thres(tmax_base,byrs)
      
    elif version=='exclude_year':
      print "No threshold file required. Percentiles will be calculated (temp)"
      print "Each of the base period years will be removed at once and the percentile calculated"
      print "Similar to bootstrap but faster, although not as robust."
      
      tminp=calc_thres(tmin_base,byrs)
      tmaxp=calc_thres(tmax_base,byrs)
      
      for yr in range(byrs):
        print "year: %s" %(yr+bsyear)
        tmax_boot=tmax_base[years_base!=yr+bsyear,:,:].copy()
        tmin_boot=tmin_base[years_base!=yr+bsyear,:,:].copy()
        
        tminpbs[yr,:,:,:,:]=calc_thres(tmin_boot,byrs-1)
        tmaxpbs[yr,:,:,:,:]=calc_thres(tmax_boot,byrs-1)
    
    else:
      
      print "Method to calculate the thresholds not supported. It must be one of the following: "
      print "bootstrap all_years exclude_year"
      sys.exit("ERROR: No valid threshold method chosen")
      
      
      
    

    for yr in range(nyears):
      year=years[0]+yr
      
      if version=="all_years":
        for month in range(1,13):
          index_ym=yr*12+month-1
          TNp[0,index_ym,:,:]=np.ma.sum(tmin[(years==year) & (months==month),:,:]<tminp[months_clim==month,0,:,:],axis=0)
          TNp[1,index_ym,:,:]=np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminp[months_clim==month,1,:,:],axis=0)
          TNp[2,index_ym,:,:]=np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminp[months_clim==month,2,:,:],axis=0)

          TXp[0,index_ym,:,:]=np.ma.sum(tmax[(years==year) & (months==month),:,:]<tmaxp[months_clim==month,0,:,:],axis=0)
          TXp[1,index_ym,:,:]=np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxp[months_clim==month,1,:,:],axis=0)
          TXp[2,index_ym,:,:]=np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxp[months_clim==month,2,:,:],axis=0)
      
      
      if version in ['bootstrap','exclude_year']: 
      
        if (year<bsyear) or (year>beyear):
          #Out of the base period
          for month in range(1,13):
            index_ym=yr*12+month-1
            #pdb.set_trace()
            TNp[0,index_ym,:,:]=np.ma.sum(tmin[(years==year) & (months==month),:,:]<tminp[months_clim==month,0,:,:],axis=0)
            TNp[1,index_ym,:,:]=np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminp[months_clim==month,1,:,:],axis=0)
            TNp[2,index_ym,:,:]=np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminp[months_clim==month,2,:,:],axis=0)
          
            TXp[0,index_ym,:,:]=np.ma.sum(tmax[(years==year) & (months==month),:,:]<tmaxp[months_clim==month,0,:,:],axis=0)
            TXp[1,index_ym,:,:]=np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxp[months_clim==month,1,:,:],axis=0)
            TXp[2,index_ym,:,:]=np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxp[months_clim==month,2,:,:],axis=0)
        
        elif (year>=bsyear) and (year<=beyear):
          #Within the base period
          for month in range(1,13):
            index_ym=yr*12+month-1
            for yr_iter in range(byrs):
              if yr_iter!=yr:
                index_ym=yr*12+month-1
  
                TNp[0,index_ym,:,:]=TNp[0,index_ym,:,:]+np.ma.sum(tmin[(years==year) & (months==month),:,:]<tminpbs[yr_iter,months_clim==month,0,:,:],axis=0)
                TNp[1,index_ym,:,:]=TNp[1,index_ym,:,:]+np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminpbs[yr_iter,months_clim==month,1,:,:],axis=0)
                TNp[2,index_ym,:,:]=TNp[2,index_ym,:,:]+np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminpbs[yr_iter,months_clim==month,2,:,:],axis=0)
          
                TXp[0,index_ym,:,:]=TXp[0,index_ym,:,:]+np.ma.sum(tmax[(years==year) & (months==month),:,:]<tmaxpbs[yr_iter,months_clim==month,0,:,:],axis=0)
                TXp[1,index_ym,:,:]=TXp[1,index_ym,:,:]+np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxpbs[yr_iter,months_clim==month,1,:,:],axis=0)
                TXp[2,index_ym,:,:]=TXp[2,index_ym,:,:]+np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxpbs[yr_iter,months_clim==month,2,:,:],axis=0)
          
        
            TNp[:,index_ym,:,:]=TNp[:,index_ym,:,:]/float(byrs-1)                            
            TXp[:,index_ym,:,:]=TXp[:,index_ym,:,:]/float(byrs-1)
        


    
    
  elif is_thresfile==1:
    print "Threshold file will be used. No percentiles will be calculated (temp)"
    thresfile=nc.Dataset(inputinf['thres_filename'],'r')
    tminp=thresfile.variables['tminp'][:]
    tmaxp=thresfile.variables['tmaxp'][:]
    
    tminpbs=thresfile.variables['tminpbs'][:]
    tmaxpbs=thresfile.variables['tmaxpbs'][:]

    thresfile.close()
    
    for yr in range(nyears):
      year=years[0]+yr
      
      if (year<bsyear) or (year>beyear):
        #Out of the base period
        for month in range(1,13):
          index_ym=yr*12+month-1
          TNp[0,index_ym,:,:]=np.ma.sum(tmin[(years==year) & (months==month),:,:]<tminp[months_clim==month,0,:,:],axis=0)
          TNp[1,index_ym,:,:]=np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminp[months_clim==month,1,:,:],axis=0)
          TNp[2,index_ym,:,:]=np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminp[months_clim==month,2,:,:],axis=0)
          
          TXp[0,index_ym,:,:]=np.ma.sum(tmax[(years==year) & (months==month),:,:]<tmaxp[months_clim==month,0,:,:],axis=0)
          TXp[1,index_ym,:,:]=np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxp[months_clim==month,1,:,:],axis=0)
          TXp[2,index_ym,:,:]=np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxp[months_clim==month,2,:,:],axis=0)
        
      else:
        #Within the base period
        for yr_iter in range(byrs):
          if yr_iter!=yr:
              TNp[0,index_ym,:,:]=TNp[0,index_ym,:,:]+np.ma.sum(tmin[(years==year) & (months==month),:,:]<tminpbs[yr_iter,months_clim==month,0,:,:],axis=0) #CHECK ALL THIS (THE MONTHS STUFF WILL GIVE AN ERROR! BECAUSE OF THE SIZE)
              TNp[1,index_ym,:,:]=TNp[1,index_ym,:,:]+np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminpbs[yr_iter,months_clim==month,1,:,:],axis=0)
              TNp[2,index_ym,:,:]=TNp[2,index_ym,:,:]+np.ma.sum(tmin[(years==year) & (months==month),:,:]>tminpbs[yr_iter,months_clim==month,2,:,:],axis=0)
          
              TXp[0,index_ym,:,:]=TXp[0,index_ym,:,:]+np.ma.sum(tmax[(years==year) & (months==month),:,:]<tmaxpbs[yr_iter,months_clim==month,0,:,:],axis=0)
              TXp[1,index_ym,:,:]=TXp[1,index_ym,:,:]+np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxpbs[yr_iter,months_clim==month,1,:,:],axis=0)
              TXp[2,index_ym,:,:]=TXp[2,index_ym,:,:]+np.ma.sum(tmax[(years==year) & (months==month),:,:]>tmaxpbs[yr_iter,months_clim==month,2,:,:],axis=0)
          
        
        TNp[:,index_ym,:,:]=TNp[:,index_ym,:,:]/float(byrs-1)                    
        TXp[:,index_ym,:,:]=TX10p[:,index_ym,:,:]/float(byrs-1)


    
  else:
    sys.exit("ERROR:Wrong entry for is_thresfile (0/1)") 
  return TNp,TXp,tminp,tmaxp,tminpbs,tmaxpbs
  
###############################################
###############################################

def calc_thres(varin,byrs):
    
    
    varinp=np.ones((365,3)+varin.shape[1:],dtype=np.float64)
    doy_s=np.tile(np.arange(365),byrs)
    use_doy=np.zeros(doy_s.shape,dtype=np.int)
    for d in range(0,365):
      use_doy[:]=1
      init_d=(d-2)%365 
      end_d=(d+2)%365
      if (init_d>end_d):
        use_doy[(doy_s<init_d) & (doy_s>=end_d)]=0
      else:
        use_doy[(doy_s<init_d) | (doy_s>=end_d)]=0
      varinp[d,:,:,:]=np.asarray(np.percentile(varin[use_doy==1,:,:],[10,50,90],axis=0))
      
      
    return varinp
        
###############################################
###############################################

def roughly_split(a, n):
  "Function to divide a list into roughly equal chunks - to be used in joblib"
  chunks=np.zeros((2,n))
  k, m = len(a) / n, len(a) % n
  chunks[0,:]=np.array(list( i * k + min(i, m) for i in xrange(n)))
  chunks[1,:]=np.array(list((i + 1) * k + min(i + 1, m) for i in xrange(n)))
  return chunks       

###############################################
###############################################



def write_fileout(ovar,varname,otime,out_file,lat,lon,inputinf):
  
  """
  Function to write out a netcdf file with the extreme variables

  """
  
  print "Writing out %s file..." %(varname)
  outfile=nc.Dataset(out_file,mode="w")
  outfile.createDimension('time',None)
  outfile.createDimension('lat',ovar.shape[1])
  outfile.createDimension('lon',ovar.shape[2])
  
  outvar=outfile.createVariable(varname,'f4',('time','lat','lon'),fill_value=const.missingval)
  outlat=outfile.createVariable('lat','f4',('lat','lon'),fill_value=const.missingval)
  outlon=outfile.createVariable('lon','f4',('lat','lon'),fill_value=const.missingval)
  outtime=outfile.createVariable('time','f4',('time'),fill_value=const.missingval)
  
  outvar[:]=ovar[:]
  outlat[:]=lat[:]
  outlon[:]=lon[:]

  outtime[:]=nc.date2num(otime,units='days since %s' %(nc.datetime.strftime(dt.datetime(1949,01,01,00), '%Y-%m-%d %H:%M:%S')),calendar='standard')
  
  setattr(outvar,"units",varinfo.get_units(varname))
  setattr(outvar,"long_name",varinfo.get_longname(varname))
  if varname=='Rnmm':
    setattr(outvar,"long_name","Days with rainfall larger than %smm" %(inputinf['Rnmm']))
  
  setattr(outlat,"standard_name","latitude")
  setattr(outlat,"long_name","latitude")
  setattr(outlat,"units","degrees_north")
  setattr(outlat,"axis","Y")
  
  setattr(outlon,"standard_name","longitude")
  setattr(outlon,"long_name","longitude")
  setattr(outlon,"units","degrees_east")
  setattr(outlon,"axis","X")
  
  setattr(outtime,"standard_name","time")
  setattr(outtime,"long_name","Time")
  setattr(outtime,"units","days since %s" %(nc.datetime.strftime(dt.datetime(1949,01,01,00), '%Y-%m-%d %H:%M:%S')))
  setattr(outtime,"calendar","standard")

  setattr(outfile,'date',dt.date.today().strftime('%Y-%m-%d'))
  setattr(outfile,'author','Daniel Argueso @CCRC UNSW')
  setattr(outfile,'contact','d.argueso@unsw.edu.au')
  setattr(outfile,'comments','Base period: %s-%s' %(inputinf['basesyear'],inputinf['baseeyear']))
  setattr(outfile,'method', 'Method to calculate thresholds: %s' %(inputinf['thres_version']))
  
  outfile.close()
  
###############################################
###############################################


def write_thresfile(tminp,tmaxp,tminpbs,tmaxpbs,prec95,prec99,lat,lon,outpath,inputinf):
  
  """
  Function to write out a netcdf file with the extreme variables

  """
  #global tmax10bs,tmax50bs,tmax90bs,tmin10bs,tmin50bs,tmin90bs,tmax10,tmax50,tmax90,tmax10,tmax50,tmax90,prec95,prec99
  print "Writing out thresholds file..." 
  outfile=nc.Dataset("%sthresholds.nc" %(outpath),mode="w")
  outfile.createDimension('time',None)
  outfile.createDimension('DoY',365)
  outfile.createDimension('perc',3)
  outfile.createDimension('lat',lat.shape[0])
  outfile.createDimension('lon',lat.shape[1])
  
  outvar_xpbs=outfile.createVariable('tmaxpbs','f4',('time','DoY','perc','lat','lon'),fill_value=const.missingval)
  outvar_npbs=outfile.createVariable('tminpbs','f4',('time','DoY','perc','lat','lon'),fill_value=const.missingval)
  

  
  outvar_x=outfile.createVariable('tmaxp','f4',('DoY','perc','lat','lon'),fill_value=const.missingval)
  outvar_n=outfile.createVariable('tminp','f4',('DoY','perc','lat','lon'),fill_value=const.missingval) 

  
  outvar_p95=outfile.createVariable('prec95','f4',('lat','lon'),fill_value=const.missingval)
  outvar_p99=outfile.createVariable('prec99','f4',('lat','lon'),fill_value=const.missingval)
  
  outlat=outfile.createVariable('lat','f4',('lat','lon'),fill_value=const.missingval)
  outlon=outfile.createVariable('lon','f4',('lat','lon'),fill_value=const.missingval)
  outtime=outfile.createVariable('time','f4',('time'),fill_value=const.missingval)
  
  outvar_xpbs[:]=tmaxpbs[:]
  outvar_npbs[:]=tminpbs[:]
  outvar_x[:]=tmaxp[:]
  outvar_n[:]=tminp[:]
  
  outvar_p95[:]=prec95[:]
  outvar_p99[:]=prec99[:]
  
  
  outlat[:]=lat[:]
  outlon[:]=lon[:]

  bsyear = int(inputinf['basesyear'])
  beyear = int(inputinf['baseeyear'])
  byrs = beyear-bsyear+1
  
  otime= [dt.datetime(bsyear+x,06,01,00) for x in range(0,byrs)]
  outtime[:]=nc.date2num(otime,units='days since %s' %(nc.datetime.strftime(dt.datetime(1949,01,01,00), '%Y-%m-%d %H:%M:%S')),calendar='standard')
  
  setattr(outlat,"standard_name","latitude")
  setattr(outlat,"long_name","latitude")
  setattr(outlat,"units","degrees_north")
  setattr(outlat,"axis","Y")
  
  setattr(outlon,"standard_name","longitude")
  setattr(outlon,"long_name","longitude")
  setattr(outlon,"units","degrees_east")
  setattr(outlon,"axis","X")
  
  setattr(outtime,"standard_name","time")
  setattr(outtime,"long_name","Time")
  setattr(outtime,"units","days since %s" %(nc.datetime.strftime(dt.datetime(1949,01,01,00), '%Y-%m-%d %H:%M:%S')))
  setattr(outtime,"calendar","standard")

  setattr(outfile,'date',dt.date.today().strftime('%Y-%m-%d'))
  setattr(outfile,'author','Daniel Argueso @CCRC UNSW')
  setattr(outfile,'contact','d.argueso@unsw.edu.au')
  setattr(outfile,'comments','Base period: %s-%s' %(inputinf['basesyear'],inputinf['baseeyear']))
  
  outfile.close()