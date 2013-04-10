;Name
; prepare_wrf_for_fclimdex.ncl
;
;Purpose
;
;Prepare data outputs from WRF to the appropriate NETCDF format for fclimdex
;
;Calling sequence 
;
; $ncl prepare_wrf_for_fclimdex.ncl
;
;Inputs:
;
; Daily data from hrly_to_dly.d.ncl
; And aggregated precipitation from make_prec.d.ncl
; they must be concatenated using nrcat
;
;Outputs:
;
; Three files: PRECIP, TMAX and TMIN.
;
;Author:
; Daniel Argüeso @UNSW
;d.argueso@unsw.edu.au
;
;Creation:
;	18/04/2012



load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"

begin

	print("_____________ START PROGRAM _________________")
	print("_______ PREPARE DATA FOR FCLIMDEX ___________")
	
	;define directories
	
	diri="/scratch/z3393020/Sydney2km_test/"
	diro="/scratch/z3393020/Sydney2km_test/"
	
;PERIOD
syear=1995
eyear=2009
	
;	READING DATA
;PRECIP


fc=addfile(diri+"PRE."+syear+"-"+eyear+".day.nc","r")
PRE    = fc->PRE
fout=addfile(diro+"WRF_PR_"+syear+"-"+eyear+".rcm.nc","c")

dsizes = dimsizes(PRE)
nlon = dsizes(2) 
nlat = dsizes(1)
ntim =  -1

dim_names = (/ "lat",  "lon", "time" /)
dim_sizes = (/ nlat ,  nlon ,  ntim  /)
dimUnlim = (/ False,  False ,  True  /)
filedimdef( fout, dim_names, dim_sizes, dimUnlim )

time = fc->time
time@units = "days since 1989-11-01 00:00:00"
time@calendar = "standard"
delete(time@long_name)
delete(time@missing_value)


; ;-GENERATING ONE-DIMENSION COORDINATES--

 xlong = fc->XLONG
 xlat = fc->XLAT
xlong!1="lon"
xlong!0="lat"
xlat!1="lon"
xlat!0="lat"
xlong@standard_name = "longitude"
xlong@long_name = "longitude"
xlong@axis = "X"
delete(xlong@stagger)
delete(xlong@MemoryOrder)
delete(xlong@description)
delete(xlong@FieldType)


xlat@standard_name = "latitude"
xlat@long_name = "latitude"
xlat@axis = "Y"
delete(xlat@stagger)
delete(xlat@MemoryOrder)
delete(xlat@description)
delete(xlat@FieldType)

; ;-------------

 PRE!0="time"
 PRE!1="lat"
 PRE!2="lon"

fout->lon = xlong
fout->lat = xlat
fout->time = time
fout->PRE=PRE

;####### TMIN

fc=addfile(diri+"T2MIN."+syear+"-"+eyear+".day.nc","r")
T2MIN   = fc->T2MIN-273.15
T2MIN@long_name = "daily min 2-m temperature"
T2MIN@units = "degC"
T2MIN@_FillValue = -999.99

fout=addfile(diro+"WRF_TN_"+syear+"-"+eyear+".rcm.nc","c")
dsizes = dimsizes(T2MIN)
nlon = dsizes(2) 
nlat = dsizes(1)
ntim =  -1

dim_names = (/ "lon",  "lat", "time" /)
dim_sizes = (/ nlon ,  nlat ,  ntim  /)
dimUnlim = (/ False,  False ,  True  /)
filedimdef( fout, dim_names, dim_sizes, dimUnlim )


T2MIN!0 = "time"
T2MIN!1 = "lat"
T2MIN!2 = "lon"
fout->lon = xlong
fout->lat = xlat
fout->time = time
fout->TN=T2MIN

;######## TMAX

fc=addfile(diri+"T2MAX."+syear+"-"+eyear+".day.nc","r")
T2MAX   = fc->T2MAX-273.15
T2MAX@long_name = "daily max 2-m temperature"
T2MAX@units = "degC"
T2MAX@_FillValue = -999.99

fout=addfile(diro+"WRF_TX_"+syear+"-"+eyear+".rcm.nc","c")

dsizes = dimsizes(T2MAX)
nlon = dsizes(2) 
nlat = dsizes(1)
ntim =  -1

dim_names = (/ "lon",  "lat", "time" /)
dim_sizes = (/ nlon ,  nlat ,  ntim  /)
dimUnlim = (/ False,  False ,  True  /)
filedimdef( fout, dim_names, dim_sizes, dimUnlim )


T2MAX!0 = "time"
T2MAX!1 = "lat"
T2MAX!2 = "lon"
fout->lon = xlong
fout->lat = xlat
fout->time = time
fout->TX=T2MAX




	print("_____________ END PROGRAM _________________")


end






