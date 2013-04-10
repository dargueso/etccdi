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
; Daily data from Lluis' Postprocessing
; 
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
	
	diri="/scratch/z3393020/ETCCDI/WRF_rcm/"
	diro="/scratch/z3393020/ETCCDI/WRF_rcm/"
	
;PERIOD
syear=1990
eyear=2009
	
;	READING DATA
;PRECIP


fc=addfile(diri+"PRE."+syear+"-"+eyear+".day.nc","r")
PRE    = fc->pracc
fout=addfile(diro+"WRF_PR_"+syear+"-"+eyear+".nc","c")

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

 xlong = fc->lon
 xlat = fc->lat
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

xlong_v=xlong(0,:)
xlat_v=xlat(:,0)
xlat_v!0="lat"
xlong_v!0="lon"


; dsizescoord=dimsizes(xlong)
; do i=0,dsizescoord(1)-1
; 	xlong(i,:)=xlong_v
; end do
; do i=0,dsizescoord(0)-1
; 	xlat(:,i)=xlat_v
; end do

; ;-------------

 PRE!0="time"
 PRE!1="lat"
 PRE!2="lon"

fout->lon = xlong_v
fout->lat = xlat_v
fout->time = time
fout->PRE=PRE

;####### TMIN
delete(dsizes)
fc=addfile(diri+"T2MIN."+syear+"-"+eyear+".day.nc","r")
T2MIN   = fc->tasmintstep-273.15
T2MIN@long_name = "daily min 2-m temperature"
T2MIN@units = "degC"
T2MIN@_FillValue = -999.99

fout=addfile(diro+"WRF_TN_"+syear+"-"+eyear+".nc","c")
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
fout->lon = xlong_v
fout->lat = xlat_v
fout->time = time
fout->TN=T2MIN

;######## TMAX
delete(dsizes)
fc=addfile(diri+"T2MAX."+syear+"-"+eyear+".day.nc","r")
T2MAX   = fc->tasmaxtstep-273.15
T2MAX@long_name = "daily max 2-m temperature"
T2MAX@units = "degC"
T2MAX@_FillValue = -999.99

fout=addfile(diro+"WRF_TX_"+syear+"-"+eyear+".nc","c")

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
fout->lon = xlong_v
fout->lat = xlat_v
fout->time = time
fout->TX=T2MAX


; 
; 
; fc=addfile(diri+"PRE."+syear+"-"+eyear+".day.nc","r")
; 
; dims = getfilevardims(fc,"PRE")
; PRE    = fc->PRE
; time = fc->time
; time@units = "days since 1995-01-01 00:00:00"
; time@calendar = "standard"
; delete(time@long_name)
; delete(time@missing_value)
; ;PREnew = PRE( $dims(2)$ | :,  $dims(1)$ | :, $dims(0)$ | :)
; ;xlong = fc->XLONG($dims(2)$ | :, $dims(1)$ | :)
; ;xlat = fc->XLAT($dims(2)$ | :, $dims(1)$ | :)
; PREnew = PRE
; xlong = fc->XLONG
; xlat = fc->XLAT
; 
; xlong!0="lon"
; xlong!1="lat"
; xlat!0="lon"
; xlat!1="lat"
; xlong@standard_name = "longitude"
; xlong@long_name = "longitude"
; xlong@axis = "X"
; delete(xlong@stagger)
; delete(xlong@MemoryOrder)
; delete(xlong@description)
; delete(xlong@FieldType)
; 
; 
; xlat@standard_name = "latitude"
; xlat@long_name = "latitude"
; xlat@axis = "Y"
; delete(xlat@stagger)
; delete(xlat@MemoryOrder)
; delete(xlat@description)
; delete(xlat@FieldType)
; 
; xlong_v=xlong(:,0)
; xlat_v=xlat(0,:)
; xlat_v!0="lat"
; xlong_v!0="lon"
; 
; 
; dsizes=dimsizes(xlong)
; do i=0,dsizes(1)-1
; 	xlong(:,i)=xlong_v
; end do
; do i=0,dsizes(0)-1
; 	xlat(i,:)=xlat_v
; end do
; 
; 
; 
; 
; 
; print(dimsizes(xlat))
; PREnew!0="time"
; PREnew!1="lat"
; PREnew!2="lon"
; 
; PREnew@_FillValue = -999.99
; 
; 
; 
; 
; 
; fout=addfile(diro+"WRF_PR_"+syear+"-"+eyear+".nc","c")
; 
; 
; fout->lon = xlong_v
; fout->lat = xlat_v
; fout->time = time
; fout->PRE=PREnew
; 
; 
; ;TMIN
; fc=addfile(diri+"T2MIN."+syear+"-"+eyear+".day.nc","r")
; 
; dims = getfilevardims(fc,"T2MIN")
; T2MIN   = fc->T2MIN
; time = fc->time
; time@units = "days since 1995-01-01 00:00:00"
; time@calendar = "standard"
; delete(time@long_name)
; delete(time@missing_value)
; ;T2MINnew= T2MIN( $dims(2)$ | :,  $dims(1)$ | :, $dims(0)$ | :)
; ;xlong = fc->XLONG
; ;xlat = fc->XLAT
; T2MINnew= T2MIN
; 
; T2MINnew!0="time"
; T2MINnew!1="lat"
; T2MINnew!2="lon"
; 
; T2MINnew@_FillValue = -999.99
; 
; fout=addfile(diro+"WRF_TN_"+syear+"-"+eyear+".nc","c")
; 
; 
; fout->lon = xlong_v
; fout->lat = xlat_v
; fout->time = time
; fout->TN=T2MINnew
; 
; ;TMAX
; fc=addfile(diri+"T2MAX."+syear+"-"+eyear+".day.nc","r")
; 
; dims = getfilevardims(fc,"T2MAX")
; T2MAX   = fc->T2MAX
; time = fc->time
; time@units = "days since 1995-01-01 00:00:00"
; time@calendar = "standard"
; delete(time@long_name)
; delete(time@missing_value)
; ;T2MAXnew= T2MAX( $dims(2)$ | :,  $dims(1)$ | :, $dims(0)$ | :)
; ;xlong = fc->XLONG
; ;xlat = fc->XLAT
; T2MAXnew= T2MAX
; 
; T2MAXnew!0="time"
; T2MAXnew!1="lat"
; T2MAXnew!2="lon"
; 
; T2MAXnew@_FillValue = -999.99
; 
; fout=addfile(diro+"WRF_TX_"+syear+"-"+eyear+".nc","c")
; 
; 
; fout->lon = xlong_v
; fout->lat = xlat_v
; fout->time = time
; fout->TX=T2MAXnew




	print("_____________ END PROGRAM _________________")


end






