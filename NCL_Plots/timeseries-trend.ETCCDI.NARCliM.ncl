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
   INDAT1=DATADIR1+"CCRC_NARCliM_R1_1950-2009_"+index+".nc"
   DATADIR2="/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"
   INDAT2=DATADIR2+"CCRC_NARCliM_R2_1950-2009_"+index+".nc"
   DATADIR3="/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"
   INDAT3=DATADIR3+"CCRC_NARCliM_R3_1950-2009_"+index+".nc"

;;;;;;;
;; fill up with e.g. diferent resolutions/parameters/...


  nyrs = 2009-1950+1
  yrs=ispan(1950,2009,1)

  file1  = addfile(INDAT1,"r")
  file2  = addfile(INDAT2,"r")
  file3  = addfile(INDAT3,"r")

  data1  = file1->Annual
  data2  = file2->Annual
  data3  = file3->Annual
  
  
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
  yAnn = new ((/3,nyrs/) , float, data1@_FillValue) 
  yAnn_11yr = yAnn

; *********************************************

; set data to missing on gridboxes where <40 years of data
   nomiss_dat1=dim_num_n(.not.ismissing(data1),0)  ; count missing values at each location
   data1 = mask( data1, (nomiss_dat1.ge.40), True)  ; mask
   nomiss_dat2=dim_num_n(.not.ismissing(data2),0)  ; count missing values at each location
   data2 = mask( data2, (nomiss_dat2.ge.40), True)  ; mask
   nomiss_dat3=dim_num_n(.not.ismissing(data3),0)  ; count missing values at each location
   data3 = mask( data3, (nomiss_dat3.ge.40), True)  ; mask
   

;;;;;;;
;; fill up for additional infiles


; Area Mean
  data1_FieldAvg = dim_avg_n_Wrap(data1, (/1,2/))
  data2_FieldAvg = dim_avg_n_Wrap(data2, (/1,2/))
  data3_FieldAvg = dim_avg_n_Wrap(data3, (/1,2/))


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

  printVarSummary(data1_FieldAvg_wgt)
  printVarSummary(yAnn)

  yAnn(0,:) = data1_FieldAvg_wgt
  yAnn(1,:) = data2_FieldAvg_wgt
  yAnn(2,:) = data3_FieldAvg_wgt
;;  yAnn(2,:) = GHCN_FieldAvg_wgt_norm1
;;  yAnn(3,:) = GHCN_FieldAvg_wgt_norm2
  yAnn_11yr(0,:) = data1_FieldAvg_wgt_11yr
  yAnn_11yr(1,:) = data2_FieldAvg_wgt_11yr
  yAnn_11yr(2,:) = data3_FieldAvg_wgt_11yr
;;  yAnn_11yr(2,:) = GHCN_FieldAvg_wgt_11yr1
;;  yAnn_11yr(3,:) = GHCN_FieldAvg_wgt_11yr2

;; output table with annual values
  outmat=new ((/nyrs,4/) , float, data1@_FillValue)
  outmat(:,0)=yrs
  outmat(:,1)= data1_FieldAvg_wgt
  outmat(:,2)= data2_FieldAvg_wgt
  outmat(:,3)= data3_FieldAvg_wgt
;;  outmat(:,3)= GHCN_FieldAvg_wgt1
  fmtf   = "f5.0,1x,3f10.2"
  opt = True
;;  opt@title  = "Year      HadEX     HadGHCN   GHCNDEX"
  opt@title  = "Year      R1     R2    R3"
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
 res@xyLineColors      = (/"blue","red","green","purple"/)
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
  sres@xyExplicitLegendLabels = (/"Data1","Data2","Data3"/)
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
   sres@xyLineColors      = (/"blue","red","green","purple"/)
;  sres@xyLineColors      = (/"blue","red","green"/)

  plot2 = gsn_csm_xy (wks,yrs,yAnn_11yr,sres) ; create plot

; ***************************************************
; overlay...
  overlay(plot1,plot2)
  draw(plot1)
  frame(wks)

end do
end