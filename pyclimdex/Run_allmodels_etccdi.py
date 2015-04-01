
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



GCM_names=['MIROC3.2','CCCMA3.1','ECHAM5','CSIRO-MK30']
RCM_names=['R1','R2','R3']
Period_names=['1990-2010','2020-2040','2060-2080']
#Period_names=['1990-2010']

Domain_names=['d01','d02']

outpath="/srv/ccrc/data14/z3393020/NARCliM/ETCCDI/Python/"
pathin="/srv/ccrc/data13/z3393020/Analyses/NARCliM/ForETCCDI/"
outpatt="CCRC_NARCliM_DAY_"
indeck="etccdi.nml.deck"

thres_version='bootstrap'

for gind,gname in enumerate(GCM_names):
  for rind,rname in enumerate(RCM_names):
    for pind,pname in enumerate(Period_names):
      for dind,dname in enumerate(Domain_names):
        fin = open (indeck,"r")
        fout = open ("etccdi_NARCliM.nml","w")
        
        print "Calculating ETCCDI for: ",gname,rname,pname,dname
        
        
        ### LINKING FILES ###
        fullpathin=str.join("/",[pathin,gname,rname,pname,dname])
        
        filename_tx=str.join("",[outpatt,pname,"_tasmax",".nc"])
        filename_tn=str.join("",[outpatt,pname,"_tasmin",".nc"])
        filename_pr=str.join("",[outpatt,pname,"_pracc_fl",".nc"])
        
        
        inpath=fullpathin
        outpath=str.join("/",[outpath,gname,rname,pname,dname,"/"])
       
        if not os.path.exists(outpath):
          os.makedirs(outpath)
        
        outname="CCRC_NARCliM"
        log_file=str.join("_",["ETCCDI",gname,rname,pname,dname,".log"])
        
        if pname=='1990-2010':
          is_thresfile=0
          thres_filename=""
        else:
          thres_path=str.join("/",[outpath,gname,rname,"1990-2010",dname])
          filename=str.join("_",[outpatt,"1990-2009","thresholds.nc"])
          thres_filename= "%s/%s" %(thres_path,filename)
          is_thresfile=1
        
        
        
        namelist_dic={'%inpath%': inpath,
                      '%outpath%': outpath,
                      '%file_prec%': filename_pr,
                      '%file_tmax%': filename_tx,
                      '%file_tmin%': filename_tn,
                      '%outname%':       outname,
                      '%save_thres%':    str(1),
                      '%log_file%':      log_file,
                      '%tmax_name%':     "tasmax",
                      '%tmin_name%':     "tasmin",
                      '%prec_name%':     "pracc_fl",
                      '%lon_username%':     "lon",
                      '%lat_username%':     "lon",
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
        
        os.system("python ./fclimdex.py -i etccdi_NARCliM.nml")
        shutil. copy ( "etccdi_NARCliM.nml" , "./etccdi_%s_%s_%s_%s.nml" %(gname,rname,pname,dname) )
                                          
















         



	
	
