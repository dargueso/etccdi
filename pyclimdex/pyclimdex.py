#!/usr/bin/env python

""" fclimdex.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Thu Feb  5 16:24:00 AEDT 2015

"""

import netCDF4 as nc
import numpy as np
from optparse import OptionParser
import etccdi_modules as em
import datetime as dt
from constants import const
import sys
from joblib import Parallel, delayed

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
fullpathin="%s/%s%s-%s_" %(inputinf['outpath'],inputinf['outname'],syear,eyear)
dates,years,months = em.calc_dates(syear,eyear)
otime_y = em.calc_otime(years,"years")
otime_m = em.calc_otime(years,"months")

###############################################
###############################################

## Calculating precipitation extremes ##

print "Loading precipitation file..."


filepr = nc.Dataset("%s/%s"%(inputinf['inpath'],inputinf['file_prec']),"r")

print "Loading ... %s/%s"%(inputinf['inpath'],inputinf['file_prec'])
prec = filepr.variables[inputinf['prec_name']][:]

time = filepr.variables[inputinf['time_username']][:]
lon  = filepr.variables[inputinf['lon_username']][:]
lat  = filepr.variables[inputinf['lat_username']][:]

if ((len(lon.shape)==1) & (len(lat.shape)==1)):
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

#Reduce variable to analysis period
prec=prec[(years_input>=syear) & (years_input<=eyear),:,:]

# R10mm,R20mm,Rnnmm,SDII=em.calc_R10mm(prec,years,float(inputinf['Rnnmm']))
# em.write_fileout(R10mm,"R10mm",otime_y,"%s%s.nc"%(fullpathin,"R10mm"),lat,lon,inputinf)
# em.write_fileout(R20mm,"R20mm",otime_y,"%s%s.nc"%(fullpathin,"R20mm"),lat,lon,inputinf)
# em.write_fileout(Rnnmm,"Rnnmm",otime_y,"%s%s.nc"%(fullpathin,"Rnnmm"),lat,lon,inputinf)
# em.write_fileout(SDII ,"SDII",otime_y, "%s%s.nc"%(fullpathin, "SDII"),lat,lon,inputinf)

# Rx1day,Rx5day=em.calc_Rx5day(prec,years,months)
# em.write_fileout(Rx1day,"Rx1day",otime_m,"%s%s.nc"%(fullpathin,"Rx1day"),lat,lon,inputinf)
# em.write_fileout(Rx5day,"Rx5day",otime_m,"%s%s.nc"%(fullpathin,"Rx5day"),lat,lon,inputinf)

# R95p,R99p,PRCPtot,prec95,prec99=em.calc_R95p(prec,years,inputinf)
# em.write_fileout(R95p,"R95p",otime_y,"%s%s.nc"%(fullpathin,"R95p"),lat,lon,inputinf)
# em.write_fileout(R99p,"R99p",otime_y,"%s%s.nc"%(fullpathin,"R99p"),lat,lon,inputinf)
# em.write_fileout(PRCPtot,"PRCPtot",otime_y,"%s%s.nc"%(fullpathin,"PRCPtot"),lat,lon,inputinf)


filepr.close()
###############################################
###############################################

## Calculating Temp extremes ##


print "Loading temp files..."
filetx = nc.Dataset("%s/%s"%(inputinf['inpath'],inputinf['file_tmax']),"r")
filetn = nc.Dataset("%s/%s"%(inputinf['inpath'],inputinf['file_tmin']),"r")


print "Loading ... %s/%s"%(inputinf['inpath'],inputinf['file_tmax'])
tmax = filetx.variables[inputinf['tmax_name']][:]
tmax_units=filetx.variables[inputinf['tmax_name']].units

print "Loading ... %s/%s"%(inputinf['inpath'],inputinf['file_tmin'])
tmin = filetn.variables[inputinf['tmin_name']][:]
tmin_units=filetn.variables[inputinf['tmin_name']].units


print "Tmax dimensions are: ", tmax.shape[:]
print "Tmin dimensions are: ", tmin.shape[:]
### Check dimensions ###

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

filetx.close()
filetn.close()
# FD,ID,SU,TR = em.calc_SU(tmax,years)
# em.write_fileout(FD,"FD",otime_y,"%s%s.nc"%(fullpathin,"FD"),lat,lon,inputinf)
# em.write_fileout(ID,"ID",otime_y,"%s%s.nc"%(fullpathin,"ID"),lat,lon,inputinf)
# em.write_fileout(SU,"SU",otime_y,"%s%s.nc"%(fullpathin,"SU"),lat,lon,inputinf)
# em.write_fileout(TR,"TR",otime_y,"%s%s.nc"%(fullpathin,"TR"),lat,lon,inputinf)

# TXx,TXn,TNn,TNx,DTR = em.calc_TXx(tmax,tmin,years,months)
# em.write_fileout(TXx,"TXx",otime_m,"%s%s.nc"%(fullpathin,"TXx"),lat,lon,inputinf)
# em.write_fileout(TXn,"TXn",otime_m,"%s%s.nc"%(fullpathin,"TXn"),lat,lon,inputinf)
# em.write_fileout(TNx,"TNx",otime_m,"%s%s.nc"%(fullpathin,"TNx"),lat,lon,inputinf)
# em.write_fileout(TNn,"TNn",otime_m,"%s%s.nc"%(fullpathin,"TNn"),lat,lon,inputinf)
# em.write_fileout(DTR,"DTR",otime_m,"%s%s.nc"%(fullpathin,"DTR"),lat,lon,inputinf)




njobs=1

nlen=np.ceil(tmax.shape[1]/np.float(njobs))
patches=np.zeros((2,njobs))
patches[0,:]=np.arange(njobs)*nlen
patches[1,:]=np.arange(1,njobs+1)*nlen
patches[1,-1]=lon.shape[0]

TNp,TXp,tminp,tmaxp,tminpbs,tmaxpbs=Parallel(n_jobs=njobs)(delayed(em.calc_TX10p)(tmax[:,patches[0,i]:patches[1,i],:],tmin[:,patches[0,i]:patches[1,i],:],dates,inputinf) for i in xrange(njobs))

em.write_fileout(TNp[0,:,:,:],"TN10p",otime_m,"%s%s.nc"%(fullpathin,"TN10p"),lat,lon,inputinf)
em.write_fileout(TNp[1,:,:,:],"TN50p",otime_m,"%s%s.nc"%(fullpathin,"TN50p"),lat,lon,inputinf)
em.write_fileout(TNp[2,:,:,:],"TN90p",otime_m,"%s%s.nc"%(fullpathin,"TN90p"),lat,lon,inputinf)
                                                                             
em.write_fileout(TXp[0,:,:,:],"TX10p",otime_m,"%s%s.nc"%(fullpathin,"TX10p"),lat,lon,inputinf)
em.write_fileout(TXp[1,:,:,:],"TX50p",otime_m,"%s%s.nc"%(fullpathin,"TX50p"),lat,lon,inputinf)
em.write_fileout(TXp[2,:,:,:],"TX90p",otime_m,"%s%s.nc"%(fullpathin,"TX90p"),lat,lon,inputinf)





###########################
#Checking missing values###
#pdb.set_trace()
# prec_missing=(np.sum(prec.mask,axis=0)/float(len(time))>inputinf['NoMissingThreshold'])
# tmax_missing=(np.sum(tmax.mask,axis=0)/float(len(time))>inputinf['NoMissingThreshold'])
# tmin_missing=(np.sum(tmin.mask,axis=0)/float(len(time))>inputinf['NoMissingThreshold'])


# em.calc_gsl(tmax,tmin)
# 

# em.calc_consecext(prec)





if int(inputinf['save_thres'])==1:
    em.write_thresfile(fullpathin)
  
  



