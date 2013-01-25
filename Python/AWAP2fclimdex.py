#!/usr/bin/env python
## 
# python script to adapt AWAP files to fclimdex
#
# D. Argueso, CCRC-UNSW, August 2012
#
## g.e. # ./AWAP2fclimdex.py -f  filename
# 
from Scientific.IO.NetCDF import NetCDFFile
import numpy
from optparse import OptionParser

### Options

parser = OptionParser()

parser.add_option("-f", "--file", dest="file",
                  help="file to be processed", metavar="FILENAME")

(opts, args) = parser.parse_args()

#######
path = opts.file


##path="/srv/ccrc/data18/z3393242/studies/NARCliMGreatSydney/postprocess/nnrp/Std_LU/1990-2010/postPROJ/"
##patt="CCRC_NARCliM_Sydney_"

f=NetCDFFile(path,'w' )
lat = f.variables['lat']
lon = f.varaibles['lon']



newlat = f.createVariable('lat', 'f', ('y','x'))
newlat = tile(lat, (1, lon.shape(1)))
newlat = newlat.T

newlon = f.createVariable('lon', 'f', ('y','x'))
newlon = tile(lon, (1, lat.shape(1)))
newlon = newlon.T


