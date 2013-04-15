;+==============================
; Daniel Argueso
; University of New South Wales
; 14 Jan 2013
; Last version: 14 Jan 2013
; email: d.argueso@unsw.edu.au
;==============================



; Script based on Markus': timeseries+trend_AUSDEX_BAK.ncl
; Purpose: Extracts information from the fclimdex_RCM output files and plot area-averaged timeseries.


;_

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;************************************************
;************************************************
undef("calc_trend")
function calc_trend (data,index)
local x,tval,nxy,numyrs,b
begin
    x = ispan(0,dimsizes(data)-1,1)*1.
    ;count the number of non-missing (annual) values at each location
    numyrs = dim_num_n(.not.ismissing(data),0)
    tval = new(1,float,1e20)
    nxy = new(1,integer)
    rc = regcoef (x,data,tval,nxy)
    rc@long_name = "linear Trend "+index
    rc=where(numyrs.lt.30,rc@_FillValue,rc)
    
    b = tval
    b = 0.5
    prob  = betainc(nxy/(nxy+tval^2),nxy/2.0,b)
    prob@_FillValue=1e20
    ; set prob to missing where rc is missing (result otherwise 0)
    prob=where(ismissing(rc),rc,prob) 
    return(/rc,prob/)
end


begin
indices=(/"Rx5day","R99p","R95p","SDII","TN10p","TX10p","TN90p","TX90p"/)
;indices=(/"CDD","CWD","GSL","PRCPTOT","R20mm","R99p","Rx1day","SDII","TNn","TR","TXn","WSDI","CSDI","DTR","FD","ID","R10mm","R95p","Rx5day","SU","TN10p","TN90p","TNx","TX10p","TX90p","TXx"/)
;indices=(/"TNn","TXn","TXx","Rx1day"/) ;;TNx
;indices=("TN10p") 

geof=addfile("/srv/ccrc/data13/z3393020/Analyses/Sydney2km/geo_em_files/geo_em.d01.narclim.AUS_ONLY.nc","r")
landmask=geof->LANDMASK(0,:,:)
lat2d=geof->XLAT_M(0,:,:)
lon2d=geof->XLONG_M(0,:,:)
nlat  = dimsizes(lat2d(:,0))
nlon  = dimsizes(lon2d(0,:))

do ii=0,dimsizes(indices)-1
  index=indices(ii)

  if (ii.gt.0) then
;   delete(data_had)
;   delete(data_ghcn)
   delete(data1)
   delete(data2)
   delete(data1_FieldAvg)
   delete(data1_FieldAvg_wgt)
   delete(data1_FieldAvg_wgt_11yr)
  end if

  
  DATADIR1="/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"
  INDAT1=DATADIR1+"CCRC_NARCliM_NNRP_R1_1950-2009_"+index+".nc"
  DATADIR2="/scratch/z3393020/ETCCDI/fclimdex_RCM_R2/index/"
  INDAT2=DATADIR2+"CCRC_NARCliM_NNRP_R2_1950-2009_"+index+".nc"
  DATADIR3="/scratch/z3393020/ETCCDI/fclimdex_RCM_R3/index/"
  INDAT3=DATADIR3+"CCRC_NARCliM_NNRP_R3_1950-2009_"+index+".nc"
  ;AWAP
  DATADIR4="/srv/ccrc/data07/z3356123/AWAP_out/indices/fclimdex/0.5deg/"
  INDAT4=DATADIR4+"AWAP_1900-2011_"+index+".nc"
  MASK4=DATADIR4+"AWAP_1900-2011_R99p.nc"
  filemask4=addfile(MASK4,"r")
  latAWAP=filemask4->lat
  lonAWAP=filemask4->lon
  maskWA=filemask4->Annual(30,:,:)
  maskWA=where(.not.ismissing(maskWA),1,maskWA)
  maskWA=where(ismissing(maskWA),0,maskWA)
   
   maskWA_WRF=rgrid2rcm(latAWAP(::-1),lonAWAP,maskWA(::-1,:),lat2d,lon2d,1)

  ;HadEX2
  DATADIR5="/srv/ccrc/data03/z3356123/HadEX2/HadEX2_3.75x2.5deg/"
  INDAT5=DATADIR5+"H2_"+index+"_1901-2010_RegularGrid_global_3.75x2.5deg_LSmask.nc"

  ; ;GHCNDEX
  ;  DATADIR5="/srv/ccrc/data07/z3356123/GHCNDEX_2012/ghcndex/gridout/GHCNDEX_2.5x2.5deg/"
  ;  INDAT5=DATADIR5+"GHCND_"+index+"_1951-2013_RegularGrid_global_2.5x2.5deg_LSmask.nc"
  
  ;NCEP 1
  DATADIR6="/srv/ccrc/data03/z3356123/NCEP1/fclimdex/indices_Jan2012/"
  INDAT6=DATADIR6+"NCEP_1948-2011_"+index+"_LSmask.nc"

   
;;;;;;;
;; fill up with e.g. diferent resolutions/parameters/...


syear=1951
eyear=2009
nyrs = eyear-syear+1
yrs=ispan(syear,eyear,1)
x = ispan(0,nyrs-1,1)

sy_NARCliM=syear-1950
sy_AWAP=syear-1900
sy_NNRP=syear-1948
;sy_GHCN=syear-1951
sy_GHCN=syear-1901 ;HadEX2

file1  = addfile(INDAT1,"r")
 file2  = addfile(INDAT2,"r")
 file3  = addfile(INDAT3,"r")
 file4  = addfile(INDAT4,"r") ;AWAP
 file5  = addfile(INDAT5,"r") ;GHCNDEX or HadEX2
 file6  = addfile(INDAT6,"r");NCEP1
 

 

 data1  = file1->Annual(sy_NARCliM:sy_NARCliM+nyrs-1,:,:)
 data2  = file2->Annual(sy_NARCliM:sy_NARCliM+nyrs-1,:,:)
 data3  = file3->Annual(sy_NARCliM:sy_NARCliM+nyrs-1,:,:)
 data4  = file4->Annual(sy_AWAP:sy_AWAP+nyrs-1,:,:)
 
 data5_aux  = file5->Ann(sy_GHCN:sy_GHCN+nyrs-1,:,:)
 data5 = data5_aux(:,{-10.4:-39},{110:156})
 lat5_red = file5->lat({-10.4:-39})
   
 data6_aux  = file6->Annual(sy_NNRP:sy_NNRP+nyrs-1,:,:)
 data6 = data6_aux(:,{-10.4:-39},{110:156})
 lat6_red = file6->lat({-10.4:-39})

;; Masking land

; printVarSummary(landmask)
; printVarSummary(data1)


data1 = mask (data1, conform(data1, landmask*maskWA_WRF, (/1,2/)), 1)  
data2 = mask (data2, conform(data2, landmask*maskWA_WRF, (/1,2/)), 1)  
data3 = mask (data3, conform(data3, landmask*maskWA_WRF, (/1,2/)), 1)
data4 = mask (data4, conform(data4, maskWA  , (/1,2/)), 1)

;;;;;;;
;; fill up for additional infiles


  ; 1 array for all timeseries to be plotted
;;  yAnn = new ((/4,nyrs/) , float, data_had@_FillValue) 
  yAnn = new ((/6,nyrs/) , float, data1@_FillValue) 
  yAnn_11yr = yAnn

; *********************************************

; set data to missing on gridboxes where <40 years of data
   ; nomiss_dat1=dim_num_n(.not.ismissing(data1),0)  ; count missing values at each location
   ;  data1 = mask( data1, (nomiss_dat1.ge.40), True)  ; mask
   ;  nomiss_dat2=dim_num_n(.not.ismissing(data2),0)  ; count missing values at each location
   ;  data2 = mask( data2, (nomiss_dat2.ge.40), True)  ; mask
   ;  nomiss_dat3=dim_num_n(.not.ismissing(data3),0)  ; count missing values at each location
   ;  data3 = mask( data3, (nomiss_dat3.ge.40), True)  ; mask
   nomiss_dat4=dim_num_n(.not.ismissing(data4),0)  ; count missing values at each location
   data4 = mask( data4, (nomiss_dat4.ge.40), True)  ; mask
   
   nomiss_dat5=dim_num_n(.not.ismissing(data5),0)  ; count missing values at each location
   data5 = mask( data5, (nomiss_dat5.ge.40), True)  ; mask

   nomiss_dat6=dim_num_n(.not.ismissing(data6),0)  ; count missing values at each location
   data6 = mask( data6, (nomiss_dat6.ge.40), True)  ; mask
   
;;;;;;;
;; fill up for additional infiles


; Area Mean
  data1_FieldAvg = dim_avg_n_Wrap(data1, (/1,2/))
  data2_FieldAvg = dim_avg_n_Wrap(data2, (/1,2/))
  data3_FieldAvg = dim_avg_n_Wrap(data3, (/1,2/))
  data4_FieldAvg = dim_avg_n_Wrap(data4, (/1,2/))
  data5_FieldAvg = dim_avg_n_Wrap(data5, (/1,2/))
  data6_FieldAvg = dim_avg_n_Wrap(data6, (/1,2/))

; weighted Area Mean (use the cosine of the latitudes for weight)
  lat1   = lat2d
  rad    = 4.0*atan(1.0)/180.0
  clat1   = cos(lat1*rad)
  data1_FieldAvg_wgt = wgt_areaave2(data1, clat1, 0) ;arg4=0->calc for non-miss values

  lat2   = lat2d
  rad    = 4.0*atan(1.0)/180.0
  clat2   = cos(lat2*rad)
  data2_FieldAvg_wgt = wgt_areaave2(data2, clat2, 0)
  
  lat3   = lat2d
  rad    = 4.0*atan(1.0)/180.0
  clat3   = cos(lat3*rad)
  data3_FieldAvg_wgt = wgt_areaave2(data3, clat3, 0)

  lat4_oned = file4->lat
  lat4 = conform_dims((/dimsizes(data4(0,:,:))/),lat4_oned,0)
  rad = 4.0*atan(1.0)/180.0
  clat4   = cos(lat4*rad)
  data4_FieldAvg_wgt = wgt_areaave2(data4, clat4, 0)
  
  lat5_oned = lat5_red
  lat5 = conform_dims((/dimsizes(data5(0,:,:))/),lat5_oned,0)
  rad = 4.0*atan(1.0)/180.0
  clat5   = cos(lat5*rad)
  data5_FieldAvg_wgt = wgt_areaave2(data5, clat5, 0)
  
  lat6_oned = lat6_red
  lat6 = conform_dims((/dimsizes(data6(0,:,:))/),lat6_oned,0)
  rad = 4.0*atan(1.0)/180.0
  clat6   = cos(lat6*rad)
  data6_FieldAvg_wgt = wgt_areaave2(data6, clat6, 0)

;;;;;;;
;; fill up for additional infiles


; remove mean 1961-1990 (index 10:39) to center around 0
;;  data1_FieldAvg_wgt_norm = data1_FieldAvg_wgt - avg(data1_FieldAvg_wgt(10:39))
;;  data2_FieldAvg_wgt_norm = data2_FieldAvg_wgt - avg(data2_FieldAvg_wgt(10:39))


;; Gaussian smoothing (11yr)
  w   = filwgts_normal (11, 1.0, 0) ;; calculate Gaussian weights
  data1_FieldAvg_wgt_11yr = wgt_runave_Wrap (data1_FieldAvg_wgt,w,0)
  data2_FieldAvg_wgt_11yr = wgt_runave_Wrap (data2_FieldAvg_wgt,w,0)
  data3_FieldAvg_wgt_11yr = wgt_runave_Wrap (data3_FieldAvg_wgt,w,0)
  data4_FieldAvg_wgt_11yr = wgt_runave_Wrap (data4_FieldAvg_wgt,w,0)
  data5_FieldAvg_wgt_11yr = wgt_runave_Wrap (data5_FieldAvg_wgt,w,0)
   data6_FieldAvg_wgt_11yr = wgt_runave_Wrap (data6_FieldAvg_wgt,w,0)

  printVarSummary(data1_FieldAvg_wgt)
  printVarSummary(yAnn)

  yAnn(0,:) = data1_FieldAvg_wgt
  yAnn(1,:) = data2_FieldAvg_wgt
  yAnn(2,:) = data3_FieldAvg_wgt
  yAnn(3,:) = data4_FieldAvg_wgt
  yAnn(4,:) = data5_FieldAvg_wgt
  yAnn(5,:) = data6_FieldAvg_wgt
;;  yAnn(2,:) = GHCN_FieldAvg_wgt_norm1
;;  yAnn(3,:) = GHCN_FieldAvg_wgt_norm2
  yAnn_11yr(0,:) = data1_FieldAvg_wgt_11yr
  yAnn_11yr(1,:) = data2_FieldAvg_wgt_11yr
  yAnn_11yr(2,:) = data3_FieldAvg_wgt_11yr
  yAnn_11yr(3,:) = data4_FieldAvg_wgt_11yr
  yAnn_11yr(4,:) = data5_FieldAvg_wgt_11yr
   yAnn_11yr(5,:) = data6_FieldAvg_wgt_11yr
;;  yAnn_11yr(2,:) = GHCN_FieldAvg_wgt_11yr1
;;  yAnn_11yr(3,:) = GHCN_FieldAvg_wgt_11yr2

;; output table with annual values
  outmat=new ((/nyrs,7/) , float, data1@_FillValue)
  outmat(:,0)=yrs
  outmat(:,1)= data1_FieldAvg_wgt
  outmat(:,2)= data2_FieldAvg_wgt
  outmat(:,3)= data3_FieldAvg_wgt
  outmat(:,4)= data4_FieldAvg_wgt
  outmat(:,5)= data5_FieldAvg_wgt
  outmat(:,6)= data6_FieldAvg_wgt
;;  outmat(:,3)= GHCN_FieldAvg_wgt1
  fmtf   = "f5.0,1x,3f10.2"
  opt = True
;;  opt@title  = "Year      HadEX     HadGHCN   GHCNDEX"
  ; opt@title  = "Year      R1        R2      R3      AWAP    GHCN    NCEP"
  opt@title  = "Year      R1        R2      R3      AWAP    HadEX2    NCEP"
  opt@fout = "/srv/ccrc/data13/z3393020/Analyses/NARCliM/ETCCDI/NARCliM_"+index+".txt"
  write_matrix (outmat, fmtf, opt)


; *****************************************************
; *** Trend calculation (of *FieldAvg_wgt)

rc_prob=calc_trend(data1_FieldAvg_wgt,index) 
print ("Trend (decadal) "+index+" Data1= "+rc_prob(0)*10+" (p="+rc_prob(1)+")")

rc_prob=calc_trend(data2_FieldAvg_wgt,index) 
print ("Trend (decadal) "+index+" Data2= "+rc_prob(0)*10+" (p="+rc_prob(1)+")")

rc_prob=calc_trend(data3_FieldAvg_wgt,index) 
print ("Trend (decadal) "+index+" Data3= "+rc_prob(0)*10+" (p="+rc_prob(1)+")")

rc_prob=calc_trend(data4_FieldAvg_wgt,index) 
print ("Trend (decadal) "+index+" Data4= "+rc_prob(0)*10+" (p="+rc_prob(1)+")")

rc_prob=calc_trend(data5_FieldAvg_wgt,index) 
print ("Trend (decadal) "+index+" Data5= "+rc_prob(0)*10+" (p="+rc_prob(1)+")")
   
rc_prob=calc_trend(data6_FieldAvg_wgt,index) 
print ("Trend (decadal) "+index+" Data6= "+rc_prob(0)*10+" (p="+rc_prob(1)+")")   

; ********************************************
; Plot Annual Values and smoothed line as Overlay...

  plotfile="/srv/ccrc/data13/z3393020/Analyses/NARCliM/ETCCDI/images/NARCliM_"+index+".png"
  wks_type = "png"
  wks_type@wkWidth = 1600
  wks_type@wkHeight = 1200
  wks = gsn_open_wks(wks_type,plotfile)

; **********************************************
; PLOT 1: Annual Values

 res                   = True
 res@gsnDraw           = False         ; Don't draw plots (for overlay)
 res@gsnFrame          = False         ; don't advance frame yet

 res@tiMainString      = "Australia Area Average "+index+" "   ; add title
 
 res@gsnMaximize   =True
 res@vpHeightF= 0.6                    ; change aspect ratio of plot
 res@vpWidthF = 0.9
 res@vpXF     = 0.05                   ; start plot at x ndc coord

;; res@gsnYRefLine           = 0.0             ; create a reference line  

 res@xyLineThicknesses = 0.5
 res@xyDashPattern = 1               ; dashed/solid(=0)
 res@trXMinF = 1950       ; set minimum X-axis value
 res@trXMaxF = 2009       ; set maximum X-axis value


; res@xyLineColors      = (/"black","blue","red","green","purple","orange"/)
 res@xyLineColors      = (/"blue","red","green","purple","orange","black"/)
; res@xyLineColors      = (/"blue","red","green"/)

 plot1  = gsn_csm_xy (wks,yrs,yAnn,res) ; create plot

; **********************************************
; PLOT 2: smoothed curves

  sres = True                                 ; set up a second resource list
  sres@gsnDraw = False                        ; do not draw the plot
  sres@gsnFrame = False                       ; do not advance the frame
  
  sres@xyLineThicknesses = (/3.0,3.0,3.0,3.0/)
  sres@xyDashPattern = 0               ; dashed/solid(=0)

; Legend with thick lines from smoothed curves
  sres@pmLegendDisplayMode    = "Always"   ; Display a legend.
;  LegendText= grids+" deg, m="+slopes
;  sres@xyExplicitLegendLabels = (/"R1","R2","R3","AWAP","GHCN","NCEP1"/)
  sres@xyExplicitLegendLabels = (/"R1","R2","R3","AWAP","HadEX2","NCEP1"/)
;  sres@xyExplicitLegendLabels = LegendText
  ;sres@lgItemOrder  = (/ 3,2,1,0 /)   ; Reorder the legends
  sres@lgItemOrder  = (/ 1,0 /)   ; Reorder the legends


  sres@pmLegendSide           = "Top"               ; Change location of
  sres@pmLegendParallelPosF   = .13 ;.2                  ; move units right
  sres@pmLegendOrthogonalPosF = -0.42 ;-0.45                ; move units down
  sres@pmLegendWidthF         = 0.08                ; Change width and
  sres@pmLegendHeightF        = 0.15 ;0.13                ; height of legend.
  sres@lgLabelFontHeightF     = .012                   ; change font height
  sres@lgPerimOn              = False                 ; no box around


;  sres@xyLineColors      = (/"black","blue","red","green","purple","orange"/)
   sres@xyLineColors      = (/"blue","red","green","purple","orange","black"/)
;  sres@xyLineColors      = (/"blue","red","green"/)

  plot2 = gsn_csm_xy (wks,yrs,yAnn_11yr,sres) ; create plot

; ***************************************************
; overlay...
  overlay(plot1,plot2)
  draw(plot1)
  frame(wks)

end do
end