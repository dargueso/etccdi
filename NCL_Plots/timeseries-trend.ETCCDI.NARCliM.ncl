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
begin

indices=(/"CDD","CWD","GSL","PRCPTOT","R20mm","R99p","Rx1day","SDII","TNn","TR","TXn","WSDI","CSDI","DTR","FD","ID","R10mm","R95p","Rx5day","SU","TN10p","TN90p","TNx","TX10p","TX90p","TXx"/)
;indices=(/"TNn","TXn","TXx","Rx1day"/) ;;TNx


geof_file=addfile("/srv/ccrc/data13/z3393020/Analyses/Sydney2km/geo_em_files/geo_em.d01.narclim.nc","r")
landmask=geof_file->LANDMASK(0,:,:)
lat2d=geof_file->XLAT_M(0,:,:)

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
   
   ;NCEP 1
   DATADIR5="/srv/ccrc/data03/z3356123/NCEP1/fclimdex/indices_Jan2012/"
   INDAT5=DATADIR5+"NCEP_1948-2011_"+index+"_LSmask.nc"
   
   ;GHCNDEX

   DATADIR6="/srv/ccrc/data07/z3356123/GHCNDEX_2012/ghcndex/gridout/GHCNDEX_2.5x2.5deg/"
   INDAT6=DATADIR6+"GHCND_"+index+"_1951-2013_RegularGrid_global_2.5x2.5deg_LSmask.nc"
   
;;;;;;;
;; fill up with e.g. diferent resolutions/parameters/...


  nyrs = 2009-1951+1
  yrs=ispan(1951,2009,1)

  file1  = addfile(INDAT1,"r")
  file2  = addfile(INDAT2,"r")
  file3  = addfile(INDAT3,"r")
  file4  = addfile(INDAT4,"r") ;AWAP
  file5  = addfile(INDAT5,"r") ;NCEP1
  file6  = addfile(INDAT6,"r")    ;GHCNDEX

  data1  = file1->Annual(1:,:,:)
  data2  = file2->Annual(1:,:,:)
  data3  = file3->Annual(1:,:,:)
  data4  = file4->Annual(51:109,:,:)
  
  data5_aux  = file5->Annual(2:61,:,:)
  data5 = data5_aux(:,{-10.4:-39},{110:156})
  lat5_red = file5->lat({-10.4:-39})
  
  data6_aux  = file6->Ann(0:59,:,:)
  data6 = data6_aux(:,{-10.4:-39},{110:156})
  lat6_red = file6->lat({-10.4:-39})
;; Masking land

printVarSummary(landmask)
printVarSummary(data1)

data1 = mask (data1, conform(data1, landmask, (/1,2/)), 1)  
data2 = mask (data2, conform(data2, landmask, (/1,2/)), 1)  
data3 = mask (data3, conform(data3, landmask, (/1,2/)), 1)   
  
;;;;;;;
;; fill up for additional infiles


  ; 1 array for all timeseries to be plotted
;;  yAnn = new ((/4,nyrs/) , float, data_had@_FillValue) 
  yAnn = new ((/5,nyrs/) , float, data1@_FillValue) 
  yAnn_11yr = yAnn

; *********************************************

; set data to missing on gridboxes where <40 years of data
   nomiss_dat1=dim_num_n(.not.ismissing(data1),0)  ; count missing values at each location
   data1 = mask( data1, (nomiss_dat1.ge.40), True)  ; mask
   nomiss_dat2=dim_num_n(.not.ismissing(data2),0)  ; count missing values at each location
   data2 = mask( data2, (nomiss_dat2.ge.40), True)  ; mask
   nomiss_dat3=dim_num_n(.not.ismissing(data3),0)  ; count missing values at each location
   data3 = mask( data3, (nomiss_dat3.ge.40), True)  ; mask
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
  outmat=new ((/nyrs,6/) , float, data1@_FillValue)
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
  opt@title  = "Year      R1        R2      R3      AWAP        NCEP    GHCN"
  opt@fout = "/srv/ccrc/data13/z3393020/Analyses/NARCliM/ETCCDI/NARCliM_"+index+".txt"
  write_matrix (outmat, fmtf, opt)


; *****************************************************
; *** Trend calculation (of *FieldAvg_wgt)

  y = data1_FieldAvg_wgt
  x = ispan(0,dimsizes(y)-1,1)*1.
  rc   = regline (x,y)
  df   = rc@nptxy-2
  prob = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.5) )
  print ("Trend (decadal) "+index+" Data1= "+rc*10+" (p="+prob+")")
  delete(x)
  delete(y)

  y = data2_FieldAvg_wgt
  x = ispan(0,dimsizes(y)-1,1)*1.
  rc   = regline (x,y)
  df   = rc@nptxy-2
  prob = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.5) )
  print ("Trend (decadal) "+index+" Data2= "+rc*10+" (p="+prob+")")
  delete(x)
  delete(y)
  
  y = data3_FieldAvg_wgt
  x = ispan(0,dimsizes(y)-1,1)*1.
  rc   = regline (x,y)
  df   = rc@nptxy-2
  prob = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.5) )
  print ("Trend (decadal) "+index+" Data3= "+rc*10+" (p="+prob+")")
  delete(x)
  delete(y)
  
  y = data4_FieldAvg_wgt
  x = ispan(0,dimsizes(y)-1,1)*1.
  rc   = regline (x,y)
  df   = rc@nptxy-2
  prob = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.5) )
  print ("Trend (decadal) "+index+" Data4= "+rc*10+" (p="+prob+")")
  delete(x)
  delete(y)
  
  y = data5_FieldAvg_wgt
  x = ispan(0,dimsizes(y)-1,1)*1.
  rc   = regline (x,y)
  df   = rc@nptxy-2
  prob = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.5) )
  print ("Trend (decadal) "+index+" Data5= "+rc*10+" (p="+prob+")")
  delete(x)
  delete(y)

  y = data6_FieldAvg_wgt
  x = ispan(0,dimsizes(y)-1,1)*1.
  rc   = regline (x,y)
  df   = rc@nptxy-2
  prob = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.6) )
  print ("Trend (decadal) "+index+" Data6= "+rc*10+" (p="+prob+")")
  delete(x)
  delete(y)
  
; ********************************************
; Plot Annual Values and smoothed line as Overlay...

  plotfile="/srv/ccrc/data13/z3393020/Analyses/NARCliM/ETCCDI/images/NARCliM_"+index+".png"
  wks_type = "png"
  wks_type@wkWidth = 800
  wks_type@wkHeight = 600
  wks = gsn_open_wks(wks_type,plotfile)

; **********************************************
; PLOT 1: Annual Values

 res                   = True
 res@gsnDraw           = False         ; Don't draw plots (for overlay)
 res@gsnFrame          = False         ; don't advance frame yet

 res@tiMainString      = "Area Average AUS"+index+" "   ; add title
 
 res@gsnMaximize   =True
 res@vpHeightF= 0.4                    ; change aspect ratio of plot
 res@vpWidthF = 0.8
 res@vpXF     = 0.1                   ; start plot at x ndc coord

;; res@gsnYRefLine           = 0.0             ; create a reference line  

 res@xyLineThicknesses = 0.5
 res@xyDashPattern = 1               ; dashed/solid(=0)
 res@trXMinF = 1950		; set minimum X-axis value
 res@trXMaxF = 2009		; set maximum X-axis value


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
  sres@xyExplicitLegendLabels = (/"R1","R2","R3","AWAP","NCEP1","GHCN"/)
;  sres@xyExplicitLegendLabels = LegendText
  ;sres@lgItemOrder  = (/ 3,2,1,0 /)   ; Reorder the legends
  sres@lgItemOrder  = (/ 1,0 /)   ; Reorder the legends


  sres@pmLegendSide           = "Top"               ; Change location of
  sres@pmLegendParallelPosF   = .13 ;.2                  ; move units right
  sres@pmLegendOrthogonalPosF = -0.42 ;-0.45                ; move units down
  sres@pmLegendWidthF         = 0.08                ; Change width and
  sres@pmLegendHeightF        = 0.11 ;0.13                ; height of legend.
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