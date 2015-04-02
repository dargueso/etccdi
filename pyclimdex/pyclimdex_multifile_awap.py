#!/usr/bin/env python

""" fclimdex.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Thu Feb  5 16:24:00 AEDT 2015

Modified: Thu Apr  2 11:27:33 AEDT 2015 - from pyclimdex.py It incorporates the possibility to have multiple input files for each variable

"""

import netCDF4 as nc
import numpy as np
from optparse import OptionParser
import etccdi_modules as em
import datetime as dt
from constants import const
import sys
from joblib import Parallel, delayed
import glob

import pdb


############### INPUT ##############

parser = OptionParser()
parser.add_option("-i","--infile",dest = "infile",
help = "namelist with input arguments", metavar = "namelist")

(opts,args) = parser.parse_args()

#####################################

inputinf=em.read_input(opts.infile)
syear = int(inputinf['syear'])
eyear = int(inputinf['eyear'])
bsyear = int(inputinf['basesyear'])
beyear = int(inputinf['baseeyear'])
byrs = beyear-bsyear+1
fullpathout="%s/%s%s-%s_" %(inputinf['outpath'],inputinf['outname'],syear,eyear)
dates,years,months = em.calc_dates(syear,eyear)
otime_y = em.calc_otime(years,"years")
otime_m = em.calc_otime(years,"months")

###############################################
###############################################

## Calculating precipitation extremes ##

print "Loading precipitation files..."

ifiles_pr_all = sorted(glob.glob("%s/%s"%(inputinf['inpath_prec'],inputinf['file_prec'])))
ifiles_pr= em.sel_files_postprocess(ifiles_pr_all,inputinf['inpattern'],syear,eyear)
print ifiles_pr 

filepr = nc.MFDataset(ifiles_pr)
print "Loading precip files between %s and %s" %(syear,eyear)


prec = filepr.variables[inputinf['prec_name']][:]
time = filepr.variables[inputinf['time_username']][:]
lon  = filepr.variables[inputinf['lon_username']][:]
lat  = filepr.variables[inputinf['lat_username']][:]

if ((len(lon.shape)==1) & (len(lat.shape)==1)):
  
  if lat[0]>lat[-1]:
    lat=lat[::-1]
    prec=prec[:,::-1,:]
  
  newlon=lon.repeat(lat.shape[0])
  newlon=np.resize(newlon,(lon.shape[0],lat.shape[0])).T
  newlat=lat.repeat(lon.shape[0])
  newlat=np.resize(newlat,(lat.shape[0],lon.shape[0]))
  lon=newlon.copy()
  lat=newlat.copy()
  
print "Preciptiation dimensions are: ", prec.shape[:]
dates_input = nc.num2date(filepr.variables[inputinf['time_username']][:],units=filepr.variables['time'].units,calendar='standard')
years_input = np.asarray([dates_input[i].year for i in xrange(len(dates_input))])

if int(inputinf['is_thresfile'])==0:
  if (bsyear<years_input[0]) or (beyear>years_input[-1]):
    print bsyear,years_input[0],beyear,years_input[-1]
    print years_input[0]-bsyear
    print beyear-years_input[-1]
    sys.exit('ERROR: Input files do not cover the base period. Need a longer input file or different base period')
  else:
    years_base = years_input.copy()
    years_base[(years>=bsyear) & ((years<=beyear))] = -999

# #Reduce variable to analysis period
# prec=prec[(years_input>=syear) & (years_input<=eyear),:,:]
# 
# ## Calculate qualitymask for precip
# prec_mask=em.calc_qualitymask(prec,years,inputinf)
# 
# 
# 
# R10mm,R20mm,Rnnmm,SDII=em.calc_R10mm(prec,years,float(inputinf['Rnnmm']))
# 
# R10mm[prec_mask==True]=const.missingval
# R20mm[prec_mask==True]=const.missingval
# Rnnmm[prec_mask==True]=const.missingval
# SDII[prec_mask==True]=const.missingval
# 
# em.write_fileout(R10mm,"R10mm",otime_y,"%s%s.nc"%(fullpathout,"R10mm"),lat,lon,inputinf)
# em.write_fileout(R20mm,"R20mm",otime_y,"%s%s.nc"%(fullpathout,"R20mm"),lat,lon,inputinf)
# em.write_fileout(Rnnmm,"Rnnmm",otime_y,"%s%s.nc"%(fullpathout,"Rnnmm"),lat,lon,inputinf)
# em.write_fileout(SDII ,"SDII",otime_y, "%s%s.nc"%(fullpathout, "SDII"),lat,lon,inputinf)
# 
# Rx1day,Rx5day=em.calc_Rx5day(prec,years,months)
# 
# Rx1day[prec_mask==True]=const.missingval
# Rx5day[prec_mask==True]=const.missingval
# 
# em.write_fileout(Rx1day,"Rx1day",otime_m,"%s%s.nc"%(fullpathout,"Rx1day"),lat,lon,inputinf)
# em.write_fileout(Rx5day,"Rx5day",otime_m,"%s%s.nc"%(fullpathout,"Rx5day"),lat,lon,inputinf)
# 
# R95p,R99p,PRCPtot,prec95,prec99=em.calc_R95p(prec,years,inputinf)
# 
# R95p[prec_mask==True]=const.missingval
# R99p[prec_mask==True]=const.missingval
# PRCPtot[prec_mask==True]=const.missingval
# 
# em.write_fileout(R95p,"R95p",otime_y,"%s%s.nc"%(fullpathout,"R95p"),lat,lon,inputinf)
# em.write_fileout(R99p,"R99p",otime_y,"%s%s.nc"%(fullpathout,"R99p"),lat,lon,inputinf)
# em.write_fileout(PRCPtot,"PRCPtot",otime_y,"%s%s.nc"%(fullpathout,"PRCPtot"),lat,lon,inputinf)
# 
# filepr.close()
###############################################
###############################################

## Calculating Temp extremes ##


print "Loading temp files..."
print "%s/%s"%(inputinf['inpath_temp'],inputinf['file_tmax'])
ifiles_tx_all = sorted(glob.glob("%s/%s"%(inputinf['inpath_temp'],inputinf['file_tmax'])))
ifiles_tn_all = sorted(glob.glob("%s/%s"%(inputinf['inpath_temp'],inputinf['file_tmin'])))

ifiles_tx=em.sel_files_postprocess(ifiles_tx_all,inputinf['inpattern'],syear,eyear)
ifiles_tn=em.sel_files_postprocess(ifiles_tn_all,inputinf['inpattern'],syear,eyear)

print ifiles_tx
print ifiles_tn


print "Loading tmax files between %s and %s" %(syear,eyear)
filetx=nc.MFDataset(ifiles_tx)
tmax = filetx.variables[inputinf['tmax_name']][:]
tmax_units=filetx.variables[inputinf['tmax_name']].units

print "Loading tmin files between %s and %s" %(syear,eyear)
filetn=nc.MFDataset(ifiles_tn)
tmin = filetn.variables[inputinf['tmin_name']][:]
tmin_units=filetn.variables[inputinf['tmin_name']].units



print "Tmax dimensions are: ", tmax.shape[:]
print "Tmin dimensions are: ", tmin.shape[:]
### Check dimensions ###

time = filetx.variables[inputinf['time_username']][:]
lon  = filetx.variables[inputinf['lon_username']][:]
lat  = filetx.variables[inputinf['lat_username']][:]



if ((len(lon.shape)==1) & (len(lat.shape)==1)):
  
  if lat[0]>lat[-1]:
    lat=lat[::-1]
    tmax=tmax[:,::-1,:]
    tmin=tmin[:,::-1,:]
    
  newlon=lon.repeat(lat.shape[0])
  newlon=np.resize(newlon,(lon.shape[0],lat.shape[0])).T
  newlat=lat.repeat(lon.shape[0])
  newlat=np.resize(newlat,(lat.shape[0],lon.shape[0]))
  lon=newlon.copy()
  lat=newlat.copy()
print lat.shape
print lon.shape
# if (tmax.shape[1:]!=lon.shape) or (tmin.shape[1:]!=lon.shape) or (lat.shape!=lon.shape):
#     sys.exit("ERROR: Coordinate dimensions in temperature do not coincide with lat/lon size (from prec file). Check dimensions for all variables")
# if (tmax.shape[0]!=len(time)):
#     sys.exit("ERROR: Time dimensions in temperature are not equal to time variable (from prec). Check first dimension for all variables")



#Reduce variables to analysis period
tmax=tmax[(years_input>=syear) & (years_input<=eyear),:,:]
tmin=tmin[(years_input>=syear) & (years_input<=eyear),:,:]

# Convert to Celsius
if tmax_units=='K':
    tmax=tmax-const.tkelvin

if tmin_units=='K':
    tmin=tmin-const.tkelvin

## Calculate qualitymask for tmax, tmin
if isinstance(tmax,np.ma.core.MaskedArray):
  tmax_mask=em.calc_qualitymask(tmax,years,inputinf)
else:
  tmax_mask=np.zeros((eyear-syear+1,)+tmax.shape[1:],dtype=np.bool)
  
if isinstance(tmin,np.ma.core.MaskedArray):
  tmin_mask=em.calc_qualitymask(tmin,years,inputinf)
else:
  tmin_mask=np.zeros((eyear-syear+1,)+tmin.shape[1:],dtype=np.bool)
  
  

filetx.close()
filetn.close()
FD,ID,SU,TR = em.calc_FD(tmax,tmin,years)

FD[tmin_mask==True]=const.missingval
ID[tmax_mask==True]=const.missingval
SU[tmax_mask==True]=const.missingval
TR[tmin_mask==True]=const.missingval

em.write_fileout(FD,"FD",otime_y,"%s%s.nc"%(fullpathout,"FD"),lat,lon,inputinf)
em.write_fileout(ID,"ID",otime_y,"%s%s.nc"%(fullpathout,"ID"),lat,lon,inputinf)
em.write_fileout(SU,"SU",otime_y,"%s%s.nc"%(fullpathout,"SU"),lat,lon,inputinf)
em.write_fileout(TR,"TR",otime_y,"%s%s.nc"%(fullpathout,"TR"),lat,lon,inputinf)

TXx,TXn,TNn,TNx,DTR = em.calc_TXx(tmax,tmin,years,months)

TXx[tmax_mask==True]=const.missingval
TXn[tmax_mask==True]=const.missingval
TNn[tmin_mask==True]=const.missingval
TNx[tmin_mask==True]=const.missingval
DTR[((tmin_mask==True) | (tmax_mask==True))]=const.missingval

em.write_fileout(TXx,"TXx",otime_m,"%s%s.nc"%(fullpathout,"TXx"),lat,lon,inputinf)
em.write_fileout(TXn,"TXn",otime_m,"%s%s.nc"%(fullpathout,"TXn"),lat,lon,inputinf)
em.write_fileout(TNx,"TNx",otime_m,"%s%s.nc"%(fullpathout,"TNx"),lat,lon,inputinf)
em.write_fileout(TNn,"TNn",otime_m,"%s%s.nc"%(fullpathout,"TNn"),lat,lon,inputinf)
em.write_fileout(DTR,"DTR",otime_m,"%s%s.nc"%(fullpathout,"DTR"),lat,lon,inputinf)






njobs=10
patches=em.roughly_split(range(tmax.shape[1]),njobs)


if njobs>1:
  
  var_in=[(tmax[:,patches[0,0]:patches[1,0],:],tmin[:,patches[0,0]:patches[1,0],:],dates,inputinf)]
  for proc in xrange(1,njobs):
    var_in.append((tmax[:,patches[0,proc]:patches[1,proc],:],tmin[:,patches[0,proc]:patches[1,proc],:],dates,inputinf))
  
  all_var=Parallel(n_jobs=njobs)(delayed(em.calc_TX10p)(*var_in[i]) for i in xrange(njobs))
  
  TNp,TXp,tminp,tmaxp,tminpbs,tmaxpbs=zip(*all_var)
  TNp=np.concatenate(TNp,axis=-2)
  TXp=np.concatenate(TXp,axis=-2)
  tminp=np.concatenate(tminp,axis=-2)
  tmaxp=np.concatenate(tmaxp,axis=-2)
  tminpbs=np.concatenate(tminpbs,axis=-2)
  tmaxpbs=np.concatenate(tmaxpbs,axis=-2)

else:
  TNp,TXp,tminp,tmaxp,tminpbs,tmaxpbs=em.calc_TX10p(tmax,tmin,dates,inputinf)
TNp[:,tmin_mask==True]=const.missingval
TXp[:,tmax_mask==True]=const.missingval

em.write_fileout(TNp[0,:,:,:],"TN10p",otime_m,"%s%s.nc"%(fullpathout,"TN10p"),lat,lon,inputinf)
em.write_fileout(TNp[1,:,:,:],"TN50p",otime_m,"%s%s.nc"%(fullpathout,"TN50p"),lat,lon,inputinf)
em.write_fileout(TNp[2,:,:,:],"TN90p",otime_m,"%s%s.nc"%(fullpathout,"TN90p"),lat,lon,inputinf)
                                                                             
em.write_fileout(TXp[0,:,:,:],"TX10p",otime_m,"%s%s.nc"%(fullpathout,"TX10p"),lat,lon,inputinf)
em.write_fileout(TXp[1,:,:,:],"TX50p",otime_m,"%s%s.nc"%(fullpathout,"TX50p"),lat,lon,inputinf)
em.write_fileout(TXp[2,:,:,:],"TX90p",otime_m,"%s%s.nc"%(fullpathout,"TX90p"),lat,lon,inputinf)





###########################
#Checking missing values###
#pdb.set_trace()
# prec_missing=(np.sum(prec.mask,axis=0)/float(len(time))>inputinf['NoMissingThreshold'])
# tmax_missing=(np.sum(tmax.mask,axis=0)/float(len(time))>inputinf['NoMissingThreshold'])
# tmin_missing=(np.sum(tmin.mask,axis=0)/float(len(time))>inputinf['NoMissingThreshold'])


# em.calc_gsl(tmax,tmin)
# 

# em.calc_consecext(prec)





if int(inputinf['save_thresholds'])==1:
    em.write_thresfile(tminp,tmaxp,tminpbs,tmaxpbs,prec95,prec99,lat,lon,fullpathout,inputinf)
  
  



