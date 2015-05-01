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
import os

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
single_file=False
calc_Pext=True
calc_Text=True
missing_vals=inputinf['missing_vals']

###############################################
###############################################

if not os.path.exists(inputinf['outpath']):
   os.makedirs(inputinf['outpath'])

###############################################
###############################################

# Calculating precipitation extremes ##
if  calc_Pext==True:
  
  if single_file:
    print "Loading precipitation file..."
    filepr = nc.Dataset("%s/%s"%(inputinf['inpath_prec'],inputinf['file_prec']),"r")
  else:
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
    print ifiles_pr

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
  prec_mask_y,prec_mask_m=em.calc_qualitymask(prec,years,inputinf)
 



  R10mm,R20mm,Rnnmm,SDII=em.calc_R10mm(prec,years,float(inputinf['Rnnmm']))

  Rx1day,Rx5day=em.calc_Rx5day(prec,years,months)

  R95p,R99p,PRCPtot,prec95,prec99=em.calc_R95p(prec,years,inputinf)

  vnamesP={'Rx1day' :Rx1day  ,
          'Rx5day' :Rx5day  , 
          'R10mm'  :R10mm   ,
          'R20mm'  :R20mm   ,
          'Rnnmm'  :Rnnmm   ,
          'SDII'   :SDII    ,
          'R95p'   :R95p    ,
          'R99p'   :R99p    ,
          'PRCPtot':PRCPtot}

  for varext in vnamesP.keys():
    
  
    if varext in ['Rx1day','Rx5day']:

      otime=otime_m
      vnamesP[varext][prec_mask_m==True]=const.missingval
    else:

      otime=otime_y
      vnamesP[varext][prec_mask_y==True]=const.missingval
    
    em.create_fileout(varext,otime,"%s%s.nc"%(fullpathout,varext),lat,lon,inputinf)
    em.put_variable(vnamesP[varext],"%s" %(varext),"%s%s.nc"%(fullpathout, varext))



  filepr.close()
##############################################
##############################################
if  calc_Text==True:
  ## Calculating Temp extremes ##

  if single_file:
    print "Loading temp files..."
    filetx = nc.Dataset("%s/%s"%(inputinf['inpath_temp'],inputinf['file_tmax']),"r")
    filetn = nc.Dataset("%s/%s"%(inputinf['inpath_temp'],inputinf['file_tmin']),"r")
    
  else:  
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
    print ifiles_tx
    filetx=nc.MFDataset(ifiles_tx)

    print "Loading tmin files between %s and %s" %(syear,eyear)
    print ifiles_tn
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
  
  if int(inputinf['is_thresfile'])==0:
    if (lat.shape[0]>300) & (lat.shape[1]>300) & (len(time)>3000):
      nsplit=20
    else:
      nsplit=1
  else:
    nsplit=1
  patches=em.roughly_split(range(lat.shape[0]),nsplit)








  for p in range(patches.shape[1]):
  
  
    print "###############################################"  
    print "LOADING PATCH %s : from %s to %s" %(p,patches[0,p],patches[1,p])
    print "###############################################"
  
    tmax = filetx.variables[inputinf['tmax_name']][:,patches[0,p]:patches[1,p],:]

    tmin = filetn.variables[inputinf['tmin_name']][:,patches[0,p]:patches[1,p],:]


    print "Tmax dimensions are: ", tmax.shape[:]
    print "Tmin dimensions are: ", tmin.shape[:]
    ### Check dimensions ###

    dates_input = nc.num2date(filetx.variables[inputinf['time_username']][:],units=filetx.variables['time'].units,calendar='standard')
    years_input = np.asarray([dates_input[i].year for i in xrange(len(dates_input))])

    #Reduce variables to analysis period
    tmax=tmax[(years_input>=syear) & (years_input<=eyear),:,:]
    tmin=tmin[(years_input>=syear) & (years_input<=eyear),:,:]
    
    # Convert to Celsius
    if hasattr(filetx.variables[inputinf['tmax_name']],'units'):
      if filetx.variables[inputinf['tmax_name']].units=='K':
        tmax=tmax-const.tkelvin
    if hasattr(filetn.variables[inputinf['tmin_name']],'units'):
      if filetn.variables[inputinf['tmin_name']]=='K':
        tmin=tmin-const.tkelvin

    ## Calculate qualitymask for tmax, tmin

    tmax_mask_y,tmax_mask_m=em.calc_qualitymask(tmax,years,inputinf)
    tmin_mask_y,tmin_mask_m=em.calc_qualitymask(tmin,years,inputinf)


    FD,ID,SU,TR = em.calc_FD(tmax,tmin,years)
  
    TXx,TXn,TNn,TNx,DTR = em.calc_TXx(tmax,tmin,years,months)

  
    ###############################################
    ###############################################
    ###PARALLEL ####
    if int(inputinf['is_thresfile'])==0:
      njobs=16
    else:
      njobs=1
    p_parallel=em.roughly_split(range(lat.shape[1]),njobs)
  
  
    if njobs>1:
  
      var_in=[(tmax[:,:,p_parallel[0,0]:p_parallel[1,0]],tmin[:,:,p_parallel[0,0]:p_parallel[1,0]],dates,inputinf)]
    
    
      for proc in xrange(1,njobs):
        var_in.append((tmax[:,:,p_parallel[0,proc]:p_parallel[1,proc]],tmin[:,:,p_parallel[0,proc]:p_parallel[1,proc]],dates,inputinf))
  
      all_var=Parallel(n_jobs=njobs)(delayed(em.calc_TX10p)(*var_in[j]) for j in xrange(njobs))
  
      TNp,TXp,tminp,tmaxp,tminpbs,tmaxpbs=zip(*all_var)
      TNp=np.concatenate(TNp,axis=-1)
      TXp=np.concatenate(TXp,axis=-1)
      tminp=np.concatenate(tminp,axis=-1)
      tmaxp=np.concatenate(tmaxp,axis=-1)
      tminpbs=np.concatenate(tminpbs,axis=-1)
      tmaxpbs=np.concatenate(tmaxpbs,axis=-1)
  
    else:
      print "NOT PARALLEL"
      TNp,TXp,tminp,tmaxp,tminpbs,tmaxpbs=em.calc_TX10p(tmax,tmin,dates,inputinf)
  
    ###############################################
    ###############################################
    vnamesT={"TN10p":TNp[0,:] ,
             "TN50p":TNp[1,:] , 
             "TN90p":TNp[2,:] ,
             "TX10p":TXp[0,:] ,
             "TX50p":TXp[1,:] ,
             "TX90p":TXp[2,:] ,
             "FD"   :FD       ,
             "ID"   :ID       ,
             "SU"   :SU       ,
             "TR"   :TR       ,
             "TXx"  :TXx      ,
             "TXn"  :TXn      ,
             "TNx"  :TNx      ,
             "TNn"  :TNn      ,
             "DTR"  :DTR      ,
       }
    if p==0:
      for varext in vnamesT.keys():
        if varext in ['FD','ID','SU','TR']:
          otime=otime_y
        else:
          otime=otime_m
        em.create_fileout(varext,otime,"%s%s.nc"%(fullpathout,varext),lat,lon,inputinf)
    
    for varext in vnamesT.keys():
      if varext in ['TN10p','TN50p','TN90p','TNx','TNn']:
        vnamesT[varext][tmin_mask_m==True]=const.missingval
      
      elif varext in  ['FD','TR']:
         vnamesT[varext][tmin_mask_y==True]=const.missingval
      
      elif varext in ['TX10p','TX50p','TX90p','TXx','TXn']:
        vnamesT[varext][tmax_mask_m==True]=const.missingval
      
      elif varext in ['ID','SU']:
        vnamesT[varext][tmax_mask_y==True]=const.missingval
        
      elif varext in ['DTR']:
        vnamesT[varext][((tmin_mask_m==True) | (tmax_mask_m==True))]=const.missingval

      em.put_variable(vnamesT[varext],"%s" %(varext),"%s%s.nc"%(fullpathout, varext),patches[:,p])
  
    if int(inputinf['save_thresholds'])==1:
      if p==0:
        em.create_thresfile(lat,lon,fullpathout,inputinf)
      em.put_variable_thfile(tminp,'tminp',"%sthresholds.nc" %(fullpathout),patches[:,p])
      em.put_variable_thfile(tminpbs,'tminpbs',"%sthresholds.nc" %(fullpathout),patches[:,p])
      em.put_variable_thfile(tmaxp,'tmaxp',"%sthresholds.nc" %(fullpathout),patches[:,p])
      em.put_variable_thfile(tmaxpbs,'tmaxpbs',"%sthresholds.nc" %(fullpathout),patches[:,p])
      if ('prec95' in locals()) & ('prec99' in locals()) & (p==0):
        em.put_variable_thfile(prec95,'prec95',"%sthresholds.nc" %(fullpathout))
        em.put_variable_thfile(prec99,'prec99',"%sthresholds.nc" %(fullpathout))
  
    # print "###############################################"  
    # print "FINISHED PATCH %s : from %s to %s" %(i,patches[0,p],patches[1,p])
    # print "###############################################" 


  filetx.close()
  filetn.close()


  
  



