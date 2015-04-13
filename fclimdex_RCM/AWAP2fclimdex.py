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
import os

### Options

parser = OptionParser()

parser.add_option("-i", "--input", dest="infile",
                  help="input filename", metavar="INFILENAME")

parser.add_option("-o", "--output", dest="outfile",
				  help="output filename", metavar="OUTFILENAME")
				
parser.add_option("-v", "--var", dest="var",
								  help="variable", metavar="VARNAME")

(opts, args) = parser.parse_args()

#######


f=NetCDFFile(opts.infile,'r')
var = f.variables[opts.var]
lat = f.variables['lat']
lon = f.variables['lon']
latID = f.dimensions['lat'] 
lonID = f.dimensions['lon'] 
allDimNames = f.dimensions.keys() 
print allDimNames[0]
print var.shape[2]



 
fo=NetCDFFile(opts.outfile,'a')
fo.createDimension('lat',lat.shape[0])
fo.createDimension('lon',lon.shape[0])
fo.createDimension('time',None)
 
 # 
 # newlat = f.createVariable('newlat', 'f', ('nlat','nlon'))
 # newlat = tile(lat, (1, lon.shape(1)))
 # newlat = newlat.T
 # 
 # newlon = f.createVariable('newlon', 'f', (allDimNames[0],allDimNames[1]))
 # newlon = tile(lon, (1, lat.shape(1)))
 # newlon = newlon.T
fo.close()

