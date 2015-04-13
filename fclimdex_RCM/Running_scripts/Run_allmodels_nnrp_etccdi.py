
#!/usr/bin/env python

""" Run_allmodels_etccdi.py

Author: Daniel Argueso @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Fri May 30 14:49:46 EST 2014

"""
import subprocess as subprocess
import os.path
import shutil



GCM_names=['NNRP']
RCM_names=['R1','R2','R3']
#Period_names=['2020-2040','2060-2080']
Period_names=['1950-2010']

Domain_names=['d01','d02']

outpath="/srv/ccrc/data14/z3393020/NARCliM/ETCCDI/"
pathin="/srv/ccrc/data13/z3393020/Analyses/NARCliM/ForETCCDI/"
patt="CCRC_NARCliM_DAY_"
indeck="input.nml.nnrp.deck"

for gind,gname in enumerate(GCM_names):
  for rind,rname in enumerate(RCM_names):
    for pind,pname in enumerate(Period_names):
      for dind,dname in enumerate(Domain_names):
        fin = open (indeck,"r")
        fout = open ("input.nml","w")
        
        print "Calculating ETCCDI for: ",gname,rname,pname,dname
        
        
        ### LINKING FILES ###
        fullpathin=str.join("/",[pathin,gname,rname,pname,dname])
        
        filename=str.join("",[patt,pname,"_tasmax",".nc"])
        subprocess.call("ln -sf %s/%s ./tmax.nc" %(fullpathin,filename),shell=True)
        
        filename=str.join("",[patt,pname,"_tasmin",".nc"])
        subprocess.call("ln -sf %s/%s ./tmin.nc" %(fullpathin,filename),shell=True)
        
        filename=str.join("",[patt,pname,"_pracc_fl",".nc"])
        subprocess.call("ln -sf %s/%s ./pre.nc" %(fullpathin,filename),shell=True)
        
        ipt_dir="./"
        opt_dir=str.join("/",[outpath,gname,rname,pname,dname,"/"])
        if not os.path.exists(opt_dir):
          os.makedirs(opt_dir)
        
        outpatt="CCRC_NARCliM"
        log_file=str.join("_",["ETCCDI",gname,rname,pname,dname,".log"])
        
        if pname=='1950-2010':
          is_thres="F"
        else:
          thres_path=str.join("/",[outpath,gname,rname,"1950-2010",dname])
          filename=str.join("_",[outpatt,"1950-2009","thresholds.nc"])
          subprocess.call("ln -sf %s/%s ./thresholds.nc" %(thres_path,filename),shell=True)
          is_thres="T"
        
        
        
        namelist_dic={'%ipt_dir%': ipt_dir,
                      '%opt_dir%': opt_dir,
                      '%tmax_filename%': "./tmax.nc",
                      '%tmin_filename%': "./tmin.nc",
                      '%prcp_filename%': "./pre.nc",
                      '%outname%':       outpatt,
                      '%save_thres%':    "T",
                      '%log_file%':      log_file,
                      '%tmax_name%':     "tasmax",
                      '%tmin_name%':     "tasmin",
                      '%prec_name%':     "pracc",
                      '%is_thres%':      is_thres}
        
        for line in fin.readlines():
          for linerep in namelist_dic.keys():
            line=line.replace(linerep,namelist_dic[linerep])
          fout.write(line)
        fin.close()
        fout.close()
        
        subprocess.call("./fclimdex.exe")
        shutil. copy ( "./input.nml" , "./input_%s_%s_%s_%s.nml" %(gname,rname,pname,dname) )
                                          
















         



	
	
