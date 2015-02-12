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
  
  if is_thresfile==0:
    print "No threshold file required. Percentiles will be calculated (precip)"
    prec_base=prec[(years>=bsyear) & (years<=beyear),:,:]
    for i in range(prec.shape[1]):
      for j in range(prec.shape[2]):
        prec_base=prec[(years>=bsyear) & (years<=beyear),i,j]
        prec_base_1mm=np.ma.masked_less_equal(prec_base,1.)
        
        prec95[i,j]=numpy.percentiles(prec_base_rain[~prec_base_rain.mask],95,axis=0)
        prec99[i,j]=numpy.percentiles(prec_base_rain[~prec_base_rain.mask],99,axis=0)
    
    
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
    R95p[yr,:,:] = np.ma.sum(prec[prec[years==year,:,:]>prec95],axis=0)
    R99p[yr,:,:] = np.ma.sum(prec[prec[years==year,:,:]>prec99],axis=0)
    PRCPtot[yr,:,:] = np.ma.sum(prec[prec[years==year,:,:]>1.],axis=0)
  
  
    
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
  
  print "Processing TN10p,TN50p,TN90p,TX10p,TX50p,TX90p ..."
  years_all=np.asarray([dates[i].year for i in xrange(len(dates))])
  months_all=np.asarray([dates[i].month for i in xrange(len(dates))])
  days_all=np.asarray([dates[i].day for i in xrange(len(dates))])
  
  #Removing leap days
  dates=dates[((months_all==2) & (days_all==29))==False]
  years=np.asarray([dates[i].year for i in xrange(len(dates))])
  months=np.asarray([dates[i].month for i in xrange(len(dates))])  
  days=np.asarray([dates[i].day for i in xrange(len(dates))])
  
  tmax=tmax[((months_all==2) & (days_all==29))==False,:,:]
  tmin=tmin[((months_all==2) & (days_all==29))==False,:,:]
  
  tmax_base=tmax[(years>=bsyear) & (years<=beyear),:,:]
  tmin_base=tmin[(years>=bsyear) & (years<=beyear),:,:]
  
  nyears=years[-1]-years[0]+1
  
  
  TXp = np.zeros((3,nyears*12,)+tmax.shape[1:],dtype=np.float)
  TNp = np.zeros((3,nyears*12,)+tmax.shape[1:],dtype=np.float)
  tminpbs=np.ones((byrs,365,3)+tmin.shape[1:],dtype=np.float)*const.missingval
  tmaxpbs=np.ones((byrs,365,3)+tmin.shape[1:],dtype=np.float)*const.missingval
 
  

  
  if is_thresfile==0:
    print "No threshold file required. Percentiles will be calculated (temp)"
    print "A bootstrap will be carried out. This part take most of the time"
    
    
    tminp=calc_thres(tmin_base,years,bsyear,beyear)
    tmaxp=calc_thres(tmax_base,years,bsyear,beyear)
    
    for yr in range(byrs+1):
      print "year: %s" %(yr+bsyear)
      tmax_boot=tmax_base.copy()
      tmin_boot=tmin_base.copy()
      for yr_iter in range(byrs+1):
        if (yr_iter!=yr):
          tmax_boot[years==yr+bsyear,:,:]=tmax_boot[years==yr_iter+bsyear,:,:]
          tmin_boot[years==yr+bsyear,:,:]=tmin_boot[years==yr_iter+bsyear,:,:]
        
        tminpbs[yr,:,:,:,:]=calc_thres(tmin_boot,years,bsyear,beyear)
        tmaxpbs[yr,:,:,:,:]=calc_thres(tmax_boot,years,bsyear,beyear)
    

    for yr in range(nyears):
      year=years[0]+yr
      
      if (year<bsyear) or (year>beyear):
        #Out of the base period
        for month in range(1,13):
          index_ym=yr*12+month-1
          TNp[0,index_ym,:,:]=np.ma.sum(tmin[(years==year & months==month),:,:]<tminp[:,0,:,:],axis=0)
          TNp[1,index_ym,:,:]=np.ma.sum(tmin[(years==year & months==month),:,:]>tminp[:,1,:,:],axis=0)
          TNp[2,index_ym,:,:]=np.ma.sum(tmin[(years==year & months==month),:,:]>tminp[:,2,:,:],axis=0)
          
          TXp[0,index_ym,:,:]=np.ma.sum(tmin[(years==year & months==month),:,:]<tmaxp[:,0,:,:],axis=0)
          TXp[1,index_ym,:,:]=np.ma.sum(tmin[(years==year & months==month),:,:]>tmaxp[:,1,:,:],axis=0)
          TXp[2,index_ym,:,:]=np.ma.sum(tmin[(years==year & months==month),:,:]>tmaxp[:,2,:,:],axis=0)
        
      else:
        #Within the base period
        for yr_iter in range(byrs):
          if yr_iter!=yr:
            TNp[0,index_ym,:,:]=TNp[0,index_ym,:,:]+np.ma.sum(tmin[(years==year & months==month),:,:]<tminpbs[yr_iter,:,0,:,:],axis=0)
            TNp[1,index_ym,:,:]=TNp[1,index_ym,:,:]+np.ma.sum(tmin[(years==year & months==month),:,:]>tminpbs[yr_iter,:,1,:,:],axis=0)
            TNp[2,index_ym,:,:]=TNp[2,index_ym,:,:]+np.ma.sum(tmin[(years==year & months==month),:,:]>tminpbs[yr_iter,:,2,:,:],axis=0)
          
            TXp[0,index_ym,:,:]=TXp[0,index_ym,:,:]+np.ma.sum(tmax[(years==year & months==month),:,:]<tmaxpbs[yr_iter,:,0,:,:],axis=0)
            TXp[1,index_ym,:,:]=TXp[1,index_ym,:,:]+np.ma.sum(tmax[(years==year & months==month),:,:]>tmaxpbs[yr_iter,:,1,:,:],axis=0)
            TXp[2,index_ym,:,:]=TXp[2,index_ym,:,:]+np.ma.sum(tmax[(years==year & months==month),:,:]>tmaxpbs[yr_iter,:,2,:,:],axis=0)
          
        
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
          TNp[0,index_ym,:,:]=np.ma.sum(tmin[(years==year & months==month),:,:]<tminp[:,0,:,:],axis=0)
          TNp[1,index_ym,:,:]=np.ma.sum(tmin[(years==year & months==month),:,:]>tminp[:,1,:,:],axis=0)
          TNp[2,index_ym,:,:]=np.ma.sum(tmin[(years==year & months==month),:,:]>tminp[:,2,:,:],axis=0)
          
          TXp[0,index_ym,:,:]=np.ma.sum(tmax[(years==year & months==month),:,:]<tmaxp[:,0,:,:],axis=0)
          TXp[1,index_ym,:,:]=np.ma.sum(tmax[(years==year & months==month),:,:]>tmaxp[:,1,:,:],axis=0)
          TXp[2,index_ym,:,:]=np.ma.sum(tmax[(years==year & months==month),:,:]>tmaxp[:,2,:,:],axis=0)
        
      else:
        #Within the base period
        for yr_iter in range(byrs):
          if yr_iter!=yr:
              TNp[0,index_ym,:,:]=TNp[0,index_ym,:,:]+np.ma.sum(tmin[(years==year & months==month),:,:]<tminpbs[yr_iter,:,0,:,:],axis=0) #CHECK ALL THIS (THE MONTHS STUFF WILL GIVE AN ERROR! BECAUSE OF THE SIZE)
              TNp[1,index_ym,:,:]=TNp[1,index_ym,:,:]+np.ma.sum(tmin[(years==year & months==month),:,:]>tminpbs[yr_iter,:,1,:,:],axis=0)
              TNp[2,index_ym,:,:]=TNp[2,index_ym,:,:]+np.ma.sum(tmin[(years==year & months==month),:,:]>tminpbs[yr_iter,:,2,:,:],axis=0)
          
              TXp[0,index_ym,:,:]=TXp[0,index_ym,:,:]+np.ma.sum(tmax[(years==year & months==month),:,:]<tmaxpbs[yr_iter,:,0,:,:],axis=0)
              TXp[1,index_ym,:,:]=TXp[1,index_ym,:,:]+np.ma.sum(tmax[(years==year & months==month),:,:]>tmaxpbs[yr_iter,:,1,:,:],axis=0)
              TXp[2,index_ym,:,:]=TXp[2,index_ym,:,:]+np.ma.sum(tmax[(years==year & months==month),:,:]>tmaxpbs[yr_iter,:,2,:,:],axis=0)
          
        
        TNp[:,index_ym,:,:]=TNp[:,index_ym,:,:]/float(byrs-1)                    
        TXp[:,index_ym,:,:]=TX10p[:,index_ym,:,:]/float(byrs-1)


    
  else:
    sys.exit("ERROR:Wrong entry for is_thresfile (0/1)") 
    
  return TNp,TXp,tminp,tmaxp,tminpbs,tmaxpbs
  
###############################################
###############################################

def calc_thres(varin,years,bsyear,beyear):
    
    
    var_base=varin[(years>=bsyear) & (years<=beyear),:,:]
    byrs=beyear-bsyear+1
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
      varinp[d,:,:,:]=np.asarray(np.percentile(var_base[use_doy==1,:,:],[10,50,90],axis=0))
      
      
    return varinp
        
        

###############################################
###############################################



def write_fileout(ovar,varname,otime,out_file,lat,lon,inputinf):
  
  """
  Function to write out a netcdf file with the extreme variables

  """
  
  print "Writing out %s file..." %(varname)
  outfile=nc.Dataset(out_file,mode="w")
  outfile.createDimension('time',None)
  outfile.createDimension('y',ovar.shape[1])
  outfile.createDimension('x',ovar.shape[2])
  
  outvar=outfile.createVariable(varname,'f4',('time','y','x'),fill_value=const.missingval)
  outlat=outfile.createVariable('lat','f4',('y','x'),fill_value=const.missingval)
  outlon=outfile.createVariable('lon','f4',('y','x'),fill_value=const.missingval)
  outtime=outfile.createVariable('time','f4',('time'),fill_value=const.missingval)
  
  outvar[:]=ovar[:]
  outlat[:]=lat[:]
  outlon[:]=lon[:]

  outtime[:]=nc.date2num(otime,units='days since %s' %(nc.datetime.strftime(dt.datetime(1949,01,01,00), '%Y-%m-%d_%H:%M:%S')),calendar='standard')
  
  setattr(outvar,"units",varinfo.get_units(varname))
  setattr(outvar,"long_name",varinfo.get_longname(varname))
  if varname=='Rnmm':
    setattr(outvar,"long_name","Days with rainfall larger than %smm" %(inputinf['Rnmm']))
  
  setattr(outlat,"standard_name","latitude")
  setattr(outlat,"long_name","Latitude")
  setattr(outlat,"units","degrees_north")
  setattr(outlat,"_CoordinateAxisType","Lat")
  setattr(outlat,"axis","Y")
  
  setattr(outlon,"standard_name","longitude")
  setattr(outlon,"long_name","Longitude")
  setattr(outlon,"units","degrees_east")
  setattr(outlon,"_CoordinateAxisType","Lon")
  setattr(outlon,"axis","X")
  
  setattr(outtime,"standard_name","time")
  setattr(outtime,"long_name","Time")
  setattr(outtime,"units","days since %s" %(nc.datetime.strftime(dt.datetime(1949,01,01,00), '%Y-%m-%d_%H:%M:%S')))
  setattr(outtime,"calendar","standard")

  setattr(outfile,'date',dt.date.today().strftime('%Y-%m-%d'))
  setattr(outfile,'author','Daniel Argueso @CCRC UNSW')
  setattr(outfile,'contact','d.argueso@unsw.edu.au')
  setattr(outfile,'comments','Base period: %s-%s' %(inputinf['basesyear'],inputinf['baseeyear']))
  
  outfile.close()
  
###############################################
###############################################


def write_thresfilet(outpath,otime,inputinf):
  
  """
  Function to write out a netcdf file with the extreme variables

  """
  #global tmax10bs,tmax50bs,tmax90bs,tmin10bs,tmin50bs,tmin90bs,tmax10,tmax50,tmax90,tmax10,tmax50,tmax90,prec95,prec99
  print "Writing out %s file..." %(varname)
  outfile=nc.Dataset("%s_thresholds.nc" %(outpath),mode="w")
  outfile.createDimension('time',None)
  outfile.createDimennsion('DoY',365)
  outfile.createDimennsion('perc',3)
  outfile.createDimension('y',ovar.shape[1])
  outfile.createDimension('x',ovar.shape[2])
  
  outvar_xpbs=outfile.createVariable('tmaxpbs','f4',('time','Doy','perc','y','x'),fill_value=const.missingval)
  outvar_npbs=outfile.createVariable('tminpbs','f4',('time','Doy','perc','y','x'),fill_value=const.missingval)
  

  
  outvar_x=outfile.createVariable('tmaxp','f4',('Doy','perc','y','x'),fill_value=const.missingval)
  outvar_n=outfile.createVariable('tminp','f4',('Doy','perc','y','x'),fill_value=const.missingval) 

  
  outvar_p95=outfile.createVariable('prec95','f4',('y','x'),fill_value=const.missingval)
  outvar_p99=outfile.createVariable('prec99','f4',('y','x'),fill_value=const.missingval)
  
  outlat=outfile.createVariable('lat','f4',('y','x'),fill_value=const.missingval)
  outlon=outfile.createVariable('lon','f4',('y','x'),fill_value=const.missingval)
  outtime=outfile.createVariable('time','f4',('time'),fill_value=const.missingval)
  
  outvar_xpbs[:]=tmaxpbs[:]
  outvar_npbs[:]=tminpbs[:]
  outvar_x[:]=tmaxp[:]
  outvar_n[:]=tminp[:]
  
  outvar_p95[:]=prec95[:]
  outvar_p99[:]=prec99[:]
  
  
  outlat[:]=lat[:]
  outlon[:]=lon[:]

  outtime[:]=nc.date2num(otime,units='days since %s' %(nc.datetime.strftime(dt.datetime(1949,01,01,00), '%Y-%m-%d_%H:%M:%S')),calendar='standard')
  
  setattr(outlat,"standard_name","latitude")
  setattr(outlat,"long_name","Latitude")
  setattr(outlat,"units","degrees_north")
  setattr(outlat,"_CoordinateAxisType","Lat")
  setattr(outlat,"axis","Y")
  
  setattr(outlon,"standard_name","longitude")
  setattr(outlon,"long_name","Longitude")
  setattr(outlon,"units","degrees_east")
  setattr(outlon,"_CoordinateAxisType","Lon")
  setattr(outlon,"axis","X")
  
  setattr(outtime,"standard_name","time")
  setattr(outtime,"long_name","Time")
  setattr(outtime,"units","days since %s" %(nc.datetime.strftime(dt.datetime(1949,01,01,00), '%Y-%m-%d_%H:%M:%S')))
  setattr(outtime,"calendar","standard")

  setattr(outfile,'date',dt.date.today().strftime('%Y-%m-%d'))
  setattr(outfile,'author','Daniel Argueso @CCRC UNSW')
  setattr(outfile,'contact','d.argueso@unsw.edu.au')
  setattr(outfile,'comments','Base period: %s-%s' %(inputinf['basesyear'],inputinf['baseeyear']))
  
  outfile.close()