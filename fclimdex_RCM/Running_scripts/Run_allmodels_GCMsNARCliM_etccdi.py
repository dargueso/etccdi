
#!/usr/bin/env python

""" Run_allmodels_etccdi.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Fri May 30 14:49:46 EST 2014

"""
import subprocess as subprocess
import os.path
import shutil
import pdb



GCM_names=['MIROC3.2','CCCMA3.1','ECHAM5','CSIRO-MK3.0']
#Period_names=['2020-2040','2060-2080']
Period_names=['1850-2100']



outpath="/srv/ccrc/data14/z3393020/NARCliM/ETCCDI/GCMs/"
pathin="/srv/ccrc/data13/z3393020/Analyses/NARCliM/ForETCCDI/GCMs/"
patt="CCRC_NARCliM_DAY_"
indeck="input.nml.GCMsNARCliM.deck"

for gind,gname in enumerate(GCM_names):
    for pind,pname in enumerate(Period_names):
        fin = open (indeck,"r")
        fout = open ("input.nml","w")
        
        print "Calculating ETCCDI for: ",gname,pname
        
        
        ### LINKING FILES ###
        #fullpathin=str.join("/",[pathin,gname])
        
        fullpathin=pathin
        
        filename=str.join("",[gname,"_tas.A2.1850-2100.AUS.nc"])
        pdb.set_trace()
        subprocess.call("ln -sf %s/%s ./tmax.nc" %(fullpathin,filename),shell=True)
        
        subprocess.call("ln -sf %s/%s ./tmin.nc" %(fullpathin,filename),shell=True)
        
        filename=str.join("",[gname,"_pr.A2.1850-2100.AUS.nc"])
        subprocess.call("ln -sf %s/%s ./pre.nc" %(fullpathin,filename),shell=True)
        
        ipt_dir="./"
        opt_dir=str.join("/",[outpath,gname,pname,"/"])
        if not os.path.exists(opt_dir):
          os.makedirs(opt_dir)
        
        outpatt="CCRC_NARCliM_GCMs"
        log_file=str.join("_",["ETCCDI",gname,pname,".log"])
        
        
        is_thres="F"
        
        
        
        namelist_dic={'%ipt_dir%': ipt_dir,
                      '%opt_dir%': opt_dir,
                      '%tmax_filename%': "./tmax.nc",
                      '%tmin_filename%': "./tmin.nc",
                      '%prcp_filename%': "./pre.nc",
                      '%outname%':       outpatt,
                      '%save_thres%':    "T",
                      '%log_file%':      log_file,
                      '%tmax_name%':     "tas",
                      '%tmin_name%':     "tas",
                      '%prec_name%':     "pr",
                      '%is_thres%':      is_thres}
        
        for line in fin.readlines():
          for linerep in namelist_dic.keys():
            line=line.replace(linerep,namelist_dic[linerep])
          fout.write(line)
        fin.close()
        fout.close()
        
        subprocess.call("./fclimdex.exe")
        shutil. copy ( "./input.nml" , "./input_%s_%s.nml" %(gname,pname) )
                                          
















         



	
	
