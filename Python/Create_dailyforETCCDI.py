#!/usr/bin/env python

""" Create_dailyforETCCDI.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Fri May 30 11:13:59 EST 2014

"""

import netCDF4 as nc
import numpy as np
from optparse import OptionParser
from datetime import date, datetime
from scipy.stats import gamma
from netCDF4 import num2date, date2num, date2index
import sys
import os
import glob
import ccrc_utils as cu
from cdo import *   # python version
cdo = Cdo()

varw="tasmin"
pathin="/home/z3393020/Analyses/NARCliM/Bias_corrected/"
GCM_names=['MIROC3.2','CCCMA3.1','ECHAM5','CSIRO-MK30']
RCM_names=['R1','R2','R3']
Period_names=['1990-2010','2020-2040','2060-2080']
Domain_names=['d01','d02']
Period_covers={'1990-2010':'1990-1999',
               '2020-2040':'2020-2029',
               '2060-2080':'2060-2069'}
pathout='/srv/ccrc/data13/z3393020/Analyses/NARCliM/ForETCCDI/'



for gind,gname in enumerate(GCM_names):
  for rind,rname in enumerate(RCM_names):
    for pind,pname in enumerate(Period_names):
      for dind,dname in enumerate(Domain_names):
        print "Creating file for: ", gname,rname,pname,dname
        print "variable: ",varw
        fullpathin=cu.get_postproc_location(gname,rname,pname)[0]
        fullpathout='%s/%s/%s/%s/%s/' %(pathout,gname,rname,pname,dname)
        if not os.path.exists(fullpathout):
          os.makedirs(fullpathout)
        
        filesin=sorted(glob.glob('%s/%s/CCRC_NARCliM_DAY_*_%s.nc' %(fullpathin,dname,varw)))  
        ofile="%s/CCRC_NARCliM_DAY_%s_%s.nc" %(fullpathout,pname,varw)
        cdo.cat(input=filesin,output=ofile)
