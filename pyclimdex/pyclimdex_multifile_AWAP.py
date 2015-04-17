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
nyears=eyear-syear+1
bsyear = int(inputinf['basesyear'])
beyear = int(inputinf['baseeyear'])
byrs = beyear-bsyear+1
fullpathout="%s/%s%s-%s_" %(inputinf['outpath'],inputinf['outname'],syear,eyear)
dates,years,months = em.calc_dates(syear,eyear)
otime_y = em.calc_otime(years,"years")
otime_m = em.calc_otime(years,"months")
nsplit=20

###############################################
###############################################

# Calculating precipitation extremes ##

print "Loading precipitation files..."

ifiles_pr_all = sorted(glob.glob("%s/%s"%(inputinf['inpath_prec'],inputinf['file_prec'])))

if inputinf['data_source']=='NARCLIM':
  
  ifiles_pr= em.sel_files_postprocess(ifiles_pr_all,inputinf['inpattern'],syear,eyear)
elif inputinf['data_source']=='AWAP':
  ifiles_pr= em.sel_files_awap(ifiles_pr_all,syear,eyear)
else:
  sys.exit('ERROR: Data source not supported for multifile')  

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

#Reduce variable to analysis period
prec=prec[(years_input>=syear) & (years_input<=eyear),:,:]

## Calculate qualitymask for precip
if isinstance(prec,np.ma.core.MaskedArray):
  prec_mask=em.calc_qualitymask(prec,years,inputinf)
else:
  prec_mask=np.zeros((eyear-syear+1,)+prec.shape[1:],dtype=np.bool)




R10mm,R20mm,Rnnmm,SDII=em.calc_R10mm(prec,years,float(inputinf['Rnnmm']))

R10mm[prec_mask==True]=const.missingval
R20mm[prec_mask==True]=const.missingval
Rnnmm[prec_mask==True]=const.missingval
SDII[prec_mask==True]=const.missingval

em.write_fileout(R10mm,"R10mm",otime_y,"%s%s.nc"%(fullpathout,"R10mm"),lat,lon,inputinf)
em.write_fileout(R20mm,"R20mm",otime_y,"%s%s.nc"%(fullpathout,"R20mm"),lat,lon,inputinf)
em.write_fileout(Rnnmm,"Rnnmm",otime_y,"%s%s.nc"%(fullpathout,"Rnnmm"),lat,lon,inputinf)
em.write_fileout(SDII ,"SDII",otime_y, "%s%s.nc"%(fullpathout, "SDII"),lat,lon,inputinf)

Rx1day,Rx5day=em.calc_Rx5day(prec,years,months)

Rx1day[prec_mask==True]=const.missingval
Rx5day[prec_mask==True]=const.missingval

em.write_fileout(Rx1day,"Rx1day",otime_m,"%s%s.nc"%(fullpathout,"Rx1day"),lat,lon,inputinf)
em.write_fileout(Rx5day,"Rx5day",otime_m,"%s%s.nc"%(fullpathout,"Rx5day"),lat,lon,inputinf)

R95p,R99p,PRCPtot,prec95,prec99=em.calc_R95p(prec,years,inputinf)

R95p[prec_mask==True]=const.missingval
R99p[prec_mask==True]=const.missingval
PRCPtot[prec_mask==True]=const.missingval

em.write_fileout(R95p,"R95p",otime_y,"%s%s.nc"%(fullpathout,"R95p"),lat,lon,inputinf)
em.write_fileout(R99p,"R99p",otime_y,"%s%s.nc"%(fullpathout,"R99p"),lat,lon,inputinf)
em.write_fileout(PRCPtot,"PRCPtot",otime_y,"%s%s.nc"%(fullpathout,"PRCPtot"),lat,lon,inputinf)




filepr.close()
##############################################
##############################################

## Calculating Temp extremes ##


print "Loading temp files..."
print "%s/%s"%(inputinf['inpath_temp'],inputinf['file_tmax'])
ifiles_tx_all = sorted(glob.glob("%s/%s"%(inputinf['inpath_temp'],inputinf['file_tmax'])))
ifiles_tn_all = sorted(glob.glob("%s/%s"%(inputinf['inpath_temp'],inputinf['file_tmin'])))

if inputinf['data_source']=='NARCLIM':
  ifiles_tx=em.sel_files_postprocess(ifiles_tx_all,inputinf['inpattern'],syear,eyear)
  ifiles_tn=em.sel_files_postprocess(ifiles_tn_all,inputinf['inpattern'],syear,eyear)
elif inputinf['data_source']=='AWAP':
  ifiles_tx=em.sel_files_awap(ifiles_tx_all,syear,eyear)
  ifiles_tn=em.sel_files_awap(ifiles_tn_all,syear,eyear)
else:
  sys.exit('ERROR: Data source not supported for multifile')



print "Loading tmax files between %s and %s" %(syear,eyear)
filetx=nc.MFDataset(ifiles_tx)

print "Loading tmin files between %s and %s" %(syear,eyear)
filetn=nc.MFDataset(ifiles_tn)

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





####### ARRAY IS VERY LARGE - REQUIRES SPLITTING #############

patches=em.roughly_split(range(lat.shape[0]),nsplit)

em.create_fileout("TN10p",otime_m,"%s%s.nc"%(fullpathout,"TN10p"),lat,lon,inputinf)
em.create_fileout("TN50p",otime_m,"%s%s.nc"%(fullpathout,"TN50p"),lat,lon,inputinf)
em.create_fileout("TN90p",otime_m,"%s%s.nc"%(fullpathout,"TN90p"),lat,lon,inputinf)
em.create_fileout("TX10p",otime_m,"%s%s.nc"%(fullpathout,"TX10p"),lat,lon,inputinf)
em.create_fileout("TX50p",otime_m,"%s%s.nc"%(fullpathout,"TX50p"),lat,lon,inputinf)
em.create_fileout("TX90p",otime_m,"%s%s.nc"%(fullpathout,"TX90p"),lat,lon,inputinf)

em.create_fileout("FD",otime_y,"%s%s.nc"%(fullpathout,"FD"),lat,lon,inputinf)
em.create_fileout("ID",otime_y,"%s%s.nc"%(fullpathout,"ID"),lat,lon,inputinf)
em.create_fileout("SU",otime_y,"%s%s.nc"%(fullpathout,"SU"),lat,lon,inputinf)
em.create_fileout("TR",otime_y,"%s%s.nc"%(fullpathout,"TR"),lat,lon,inputinf)

em.create_fileout("TXx",otime_m,"%s%s.nc"%(fullpathout,"TXx"),lat,lon,inputinf)
em.create_fileout("TXn",otime_m,"%s%s.nc"%(fullpathout,"TXn"),lat,lon,inputinf)
em.create_fileout("TNx",otime_m,"%s%s.nc"%(fullpathout,"TNx"),lat,lon,inputinf)
em.create_fileout("TNn",otime_m,"%s%s.nc"%(fullpathout,"TNn"),lat,lon,inputinf)
em.create_fileout("DTR",otime_m,"%s%s.nc"%(fullpathout,"DTR"),lat,lon,inputinf)

for p in range(patches.shape[1]):
  
  
  print "###############################################"  
  print "LOADING PATCH %s : from %s to %s" %(p,patches[0,p],patches[1,p])
  print "###############################################"
  
  tmax = filetx.variables[inputinf['tmax_name']][:,patches[0,p]:patches[1,p],:]
  tmax_units=filetx.variables[inputinf['tmax_name']].units
  
  tmin = filetn.variables[inputinf['tmin_name']][:,patches[0,p]:patches[1,p],:]
  tmin_units=filetn.variables[inputinf['tmin_name']].units

  print "Tmax dimensions are: ", tmax.shape[:]
  print "Tmin dimensions are: ", tmin.shape[:]
  ### Check dimensions ###



  dates_input = nc.num2date(filetx.variables[inputinf['time_username']][:],units=filetx.variables['time'].units,calendar='standard')
  years_input = np.asarray([dates_input[i].year for i in xrange(len(dates_input))])



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
  
  FD,ID,SU,TR = em.calc_FD(tmax,tmin,years)
  FD[tmin_mask==True]=const.missingval
  ID[tmax_mask==True]=const.missingval
  SU[tmax_mask==True]=const.missingval
  TR[tmin_mask==True]=const.missingval

  
  em.put_variable(FD,"FD","%s%s.nc"%(fullpathout,"FD"),patches[:,p])
  em.put_variable(ID,"ID","%s%s.nc"%(fullpathout,"ID"),patches[:,p])
  em.put_variable(SU,"SU","%s%s.nc"%(fullpathout,"SU"),patches[:,p])
  em.put_variable(TR,"TR","%s%s.nc"%(fullpathout,"TR"),patches[:,p])

  TXx,TXn,TNn,TNx,DTR = em.calc_TXx(tmax,tmin,years,months)

  TXx[tmax_mask==True]=const.missingval
  TXn[tmax_mask==True]=const.missingval
  TNn[tmin_mask==True]=const.missingval
  TNx[tmin_mask==True]=const.missingval
  DTR[((tmin_mask==True) | (tmax_mask==True))]=const.missingval

  em.put_variable(TXx,"TXx","%s%s.nc"%(fullpathout,"TXx"),patches[:,p])
  em.put_variable(TXn,"TXn","%s%s.nc"%(fullpathout,"TXn"),patches[:,p])
  em.put_variable(TNx,"TNx","%s%s.nc"%(fullpathout,"TNx"),patches[:,p])
  em.put_variable(TNn,"TNn","%s%s.nc"%(fullpathout,"TNn"),patches[:,p])
  em.put_variable(DTR,"DTR","%s%s.nc"%(fullpathout,"DTR"),patches[:,p])
  
  
  
  # ###PARALLEL ####
  # if int(inputinf['is_thresfile'])==0:
  #   njobs=16
  # else:
  #   njobs=1
  # p_parallel=em.roughly_split(range(lat.shape[1]),njobs)
  # 
  # 
  # if njobs>1:
  # 
  #   var_in=[(tmax[:,:,p_parallel[0,0]:p_parallel[1,0]],tmin[:,:,p_parallel[0,0]:p_parallel[1,0]],dates,inputinf)]
  #   
  #   
  #   for proc in xrange(1,njobs):
  #     var_in.append((tmax[:,:,p_parallel[0,proc]:p_parallel[1,proc]],tmin[:,:,p_parallel[0,proc]:p_parallel[1,proc]],dates,inputinf))
  # 
  #   all_var=Parallel(n_jobs=njobs)(delayed(em.calc_TX10p)(*var_in[j]) for j in xrange(njobs))
  # 
  #   TNp,TXp,tminp,tmaxp,tminpbs,tmaxpbs=zip(*all_var)
  #   TNp=np.concatenate(TNp,axis=-1)
  #   TXp=np.concatenate(TXp,axis=-1)
  #   tminp=np.concatenate(tminp,axis=-1)
  #   tmaxp=np.concatenate(tmaxp,axis=-1)
  #   tminpbs=np.concatenate(tminpbs,axis=-1)
  #   tmaxpbs=np.concatenate(tmaxpbs,axis=-1)
  # 
  # else:
  #   TNp,TXp,tminp,tmaxp,tminpbs,tmaxpbs=em.calc_TX10p(tmax,tmin,dates,inputinf)
  # TNp[:,tmin_mask==True]=const.missingval
  # TXp[:,tmax_mask==True]=const.missingval
  # 
  # 
  # 
  # 
  # 
  # em.put_variable(TNp[0,:],"TN10p","%s%s.nc"%(fullpathout,"TN10p"),patches[:,p])
  # em.put_variable(TNp[1,:],"TN50p","%s%s.nc"%(fullpathout,"TN50p"),patches[:,p])
  # em.put_variable(TNp[2,:],"TN90p","%s%s.nc"%(fullpathout,"TN90p"),patches[:,p])
  # em.put_variable(TXp[0,:],"TX10p","%s%s.nc"%(fullpathout,"TX10p"),patches[:,p])
  # em.put_variable(TXp[1,:],"TX50p","%s%s.nc"%(fullpathout,"TX50p"),patches[:,p])
  # em.put_variable(TXp[2,:],"TX90p","%s%s.nc"%(fullpathout,"TX90p"),patches[:,p])
  # 
  # 
  # if int(inputinf['save_thresholds'])==1:
  #   em.create_thresfile(lat,lon,fullpathout,inputinf)
  #   em.put_variable_thfile(tminp,'tminp',"%sthresholds.nc" %(fullpathout),patches[:,p])
  #   em.put_variable_thfile(tminpbs,'tminpbs',"%sthresholds.nc" %(fullpathout),patches[:,p])
  #   em.put_variable_thfile(tmaxp,'tmaxp',"%sthresholds.nc" %(fullpathout),patches[:,p])
  #   em.put_variable_thfile(tmaxpbs,'tmaxpbs',"%sthresholds.nc" %(fullpathout),patches[:,p])
  # 
  # print "###############################################"  
  # print "FINISHED PATCH %s : from %s to %s" %(i,patches[0,p],patches[1,p])
  # print "###############################################" 


filetx.close()
filetn.close()

if int(inputinf['save_thresholds'])==1:
  if not os.path.exists("%sthresholds.nc" %(fullpathout)):
    em.create_thresfile(lat,lon,fullpathout,inputinf)
  if ('prec95' in locals()) & ('prec99' in locals()):
    em.put_variable_thfile(prec95,'prec95',"%sthresholds.nc" %(fullpathout),patches[:,p])
    em.put_variable_thfile(prec99,'prec99',"%sthresholds.nc" %(fullpathout),patches[:,p])

  
  



