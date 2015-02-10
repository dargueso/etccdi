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

varw="tasmax"
GCM_names=['NNRP']
RCM_names=['R1','R2','R3']
Period_names=['1950-2010']
Domain_names=['d01','d02']
pathout='/srv/ccrc/data13/z3393020/Analyses/NARCliM/ForETCCDI/'



for gind,gname in enumerate(GCM_names):
  for rind,rname in enumerate(RCM_names):
    for pind,pname in enumerate(Period_names):
      for dind,dname in enumerate(Domain_names):
        print "Creating file for: ", gname,rname,pname,dname
        print "variable: ",varw
        if varw == 'pracc_fl':
          fullpathin='/srv/ccrc/data30/z3393020/NARCliM/filtered/%s/%s/%s/'%(gname,rname,pname)
        else:  
          fullpathin=cu.get_postproc_location(gname,rname,pname)[0]
        
        
        print '%s/%s/CCRC_NARCliM_DAY_*_%s.nc' %(fullpathin,dname,varw)
        fullpathout='%s/%s/%s/%s/%s/' %(pathout,gname,rname,pname,dname)
        if not os.path.exists(fullpathout):
          os.makedirs(fullpathout)
        
        filesin=sorted(glob.glob('%s/%s/CCRC_NARCliM_DAY_*_%s.nc' %(fullpathin,dname,varw)))  
        if os.path.exists("%s/CCRC_NARCliM_DAY_%s_%s.nc" %(fullpathout,pname,varw)):
          os.remove("%s/CCRC_NARCliM_DAY_%s_%s.nc" %(fullpathout,pname,varw))
          
        ofile="%s/CCRC_NARCliM_DAY_%s_%s.nc" %(fullpathout,pname,varw)
        cdo.cat(input=filesin,output=ofile)
