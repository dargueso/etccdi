
#!/usr/bin/env python

""" Run_allmodels_etccdi.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Fri May 30 14:49:46 EST 2014

"""
import subprocess as subprocess
import os.path
import shutil
import os
import ccrc_utils as cu
import pdb


# GCM_names=['NNRP']
# RCM_names=['R1','R2','R3']
# Period_names=['1950-2010']#,'2020-2040','2060-2080']

GCM_names=['MIROC3.2','CCCMA3.1','ECHAM5','CSIRO-MK3.0']
RCM_names=['R1','R2','R3']
Period_names=['1990-2010','2020-2040','2060-2080']

Domain_names=['d01','d02']
#data_type='bc'



inpattern="CCRC_NARCliM_DAY_"
indeck="etccdi_multifile.nml.deck"
for data_type in ['pp','bc']:
  if data_type=='pp':
    outpath_generic="/srv/ccrc/data14/z3393020/NARCliM/ETCCDI/Raw/"
    txvarname='tasmax'
    tnvarname='tasmin'
    prvarname='pracc_fl'
    data_type_pr='flt'
  
  elif data_type=='bc':
    outpath_generic="/srv/ccrc/data14/z3393020/NARCliM/ETCCDI/Bias-corrected/"
    txvarname='tasmax_bc'
    tnvarname='tasmin_bc'
    prvarname='pracc_bc'
    data_type_pr=data_type

  thres_version='bootstrap'

  for gind,gname in enumerate(GCM_names):
    for rind,rname in enumerate(RCM_names):
      for pind,pname in enumerate(Period_names):
        for dind,dname in enumerate(Domain_names):
          fin = open (indeck,"r")
          fout = open ("etccdi_NARCliM.nml","w")
        
          print "Calculating ETCCDI for: ",gname,rname,pname,dname
        
        
          ### LINKING FILES ###
        
        
          filename_tx=str.join("",[inpattern,"*",txvarname,".nc"])
          filename_tn=str.join("",[inpattern,"*",tnvarname,".nc"])
          filename_pr=str.join("",[inpattern,"*",prvarname,".nc"])
        
        
          outpath=str.join("/",[outpath_generic,gname,rname,pname,dname,"/"])
       
          if not os.path.exists(outpath):
            os.makedirs(outpath)
        
          outname="CCRC_NARCliM"
          log_file=str.join("_",["ETCCDI",gname,rname,pname,dname,".log"])
        
          if pname=='1990-2010':
            is_thresfile=0
            thres_filename=""
          else:
            thres_path=str.join("/",[outpath_generic,gname,rname,"1990-2010",dname])
            filename=str.join("_",[outname,"1990-2009","thresholds.nc"])
            thres_filename= "%s/%s" %(thres_path,filename)
            is_thresfile=1

        
          print cu.get_location(gname,rname,pname,data_type)[0]
        
          namelist_dic={'%inpath_temp%': str.join("/",[cu.get_location(gname,rname,pname,data_type)[0],dname]),
                        '%inpath_prec%': str.join("/",[cu.get_location(gname,rname,pname,data_type_pr)[0],dname]),
                        '%outpath%': outpath,
                        '%inpattern%':  inpattern,
                        '%data_source%':       "NARCLIM",
                        '%file_prec%': filename_pr,
                        '%file_tmax%': filename_tx,
                        '%file_tmin%': filename_tn,
                        '%outname%':       outname,
                        '%save_thres%':    str(1),
                        '%log_file%':      log_file,
                        '%tmax_name%':     txvarname,
                        '%tmin_name%':     tnvarname,
                        '%prec_name%':     prvarname,
                        '%lon_username%':     "lon",
                        '%lat_username%':     "lat",
                        '%time_username%':     "time",
                        '%is_thresfile%':      str(is_thresfile),
                        '%thres_filename%': thres_filename,
                        '%thres_version%': thres_version,
                        '%is_rcm%':       str(1),
                        '%syear%':        str(int(pname[0:4])),
                        '%eyear%':        str(int(pname[5:9])-1),
                        '%basesyear%':    str(1990),
                        '%baseeyear%':    str(2009)}
        
          for line in fin.readlines():
            for linerep in namelist_dic.keys():
              line=line.replace(linerep,namelist_dic[linerep])
            fout.write(line)
          fin.close()
          fout.close()
        
          os.system("python ./pyclimdex_multifile.py -i etccdi_NARCliM.nml")
          #shutil. copy ( "etccdi_NARCliM.nml" , "./etccdi_%s_%s_%s_%s.nml" %(gname,rname,pname,dname) )
                                          
















         



	
	
