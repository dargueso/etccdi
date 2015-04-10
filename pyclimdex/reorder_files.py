#!/usr/bin/env python

""" Script_name.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Fri Apr 10 09:27:58 AEST 2015

"""
import subprocess as subprocess
import os.path
import shutil
import os
import ccrc_utils as cu
import pdb


GCM_names=['MIROC3.2','CCCMA3.1','ECHAM5']#,'CSIRO-MK3.0']
RCM_names=['R1','R2','R3']
Period_names=['1990-2010','2020-2040','2060-2080']
#Period_names=['2020-2040']

Domain_names=['d01','d02']

outpath_old="/srv/ccrc/data14/z3393020/NARCliM/ETCCDI/Python/"
inpattern="CCRC_NARCliM_DAY_"
indeck="etccdi_multifile.nml.deck"

for gind,gname in enumerate(GCM_names):
  for rind,rname in enumerate(RCM_names):
    for pind,pname in enumerate(Period_names):
      for dind,dname in enumerate(Domain_names):
        
        outpath_generic="/srv/ccrc/data14/z3393020/NARCliM/ETCCDI/Python/"
        outpath_old=str.join("/",[outpath_old,gname,rname,pname,dname,"/"])
        outpath_new=str.join("/",[outpath_generic,gname,rname,pname,dname,"/"])
        if not os.path.exists(outpath_new):
          os.makedirs(outpath_new)
        
        print outpath_old,outpath_new

        os.system("mv %s/*nc %s/" %(outpath_old,outpath_new))  
        