;; ncl script to calculate and plot 
;; the climate mean 1961-1990 and trends since 1950 of all indices
;; 
;; Markus Donat, 18 August 2011
;; modified 11 Jan 2012
;; m.donat@unsw.edu.au

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;************************************************
begin

gridsizes=(/"0.5","1.0","2.0","3.0","5.0"/)
;;mparams=
;;distances=

;;indices=(/"CDD","CWD","GSL","PRCPTOT","R20mm","R99p","Rx1day","SDII","TN50p","TNn","TR","TX50p","TXn","WSDI","CSDI","DTR","FD","ID","R10mm","R95p","R95pT","Rnnmm","Rx5day","SU","TN10p","TN90p","TNx","TX10p","TX90p","TXx"/)
indices=(/"Rx1day","TNn","TXx","TNx","TXn"/)


do igrid=0,dimsizes(gridsizes)-1
gridsize=gridsizes(igrid)
do ii=0,dimsizes(indices)-1
  index=indices(ii)
   print ("plotting index "+index+"...")
   DATADIR="../Gridding/gridout/AUSDEX_HQ_"+gridsize+"deg_2000km_m1/"

   INDAT=DATADIR+"AUSDEX_HQ_"+index+"_1951-2012_RegularGrid_AUS_"+gridsize+"deg_m1_LSmask.nc"

   PLOTDIR=DATADIR+"PLOTS"
   system("mkdir -p "+PLOTDIR)
   plotfile1=PLOTDIR+"/AUSDEX_"+index+"_trend_1951-2011_"+gridsize+"deg"
   plotfile2=PLOTDIR+"/AUSDEX_"+index+"_mean_1961-1990_"+gridsize+"deg"

;;   if (igrid.gt.0) then
;;     delete(data1)
;;
;;  end if

   if (ii.gt.0.or.igrid.gt.0) then
     delete(file1)
     delete(data1)
     delete(y)
     delete(x)
     delete(tval)
     delete(nxy)
     delete(numyrs)
     delete(rc)
     delete(b)
     delete(prob)
     delete(TimAvg)
     delete(res)
     delete(sres)
     delete(res3)
   end if

   file1  = addfile(INDAT,"r")
   data1  = file1->Ann

printVarSummary(data1)

; *******************************************************
; trend calculation (where at least 10 data points available)

   y=data1(:,:,:) ; years 1951-2012
   x    = ispan(0,dimsizes(y&time)-1,1)*1.

    dimy = dimsizes(y)    ; determine the dimensions of "y"
    ntim = dimy(0)
    nlat = dimy(1)
    nlon = dimy(2)

    ;count the number of non-missing (annual) values at each location
    numyrs = dim_num_n(.not.ismissing(y),0)

    tval = new ((/nlat,nlon/) , float,-99.99)
;    tval@_FillValue=-99.99
    nxy  = new ((/nlat,nlon/) , integer)

    rc = y(0,:,:)        ;copy coordinates + attributes to rc
    rc   = regcoef (x,y(lat|:,lon|:,time|:), tval,nxy)
    ;;rc10  = where(numyrs.lt.10,rc@_FillValue,regcoef (x,y(lats|:,lons|:,time|:), tval,nxy))
    rc@long_name = "linear Trend "+index

    ; set rc to missing where less than 30 data points contribute to regression
    rc=where(numyrs.lt.30,rc@_FillValue,rc)

    b  = tval            ; b must be same size as tval (and df)
    b = 0.5
    prob  = betainc(nxy/(nxy+tval^2),nxy/2.0,b)
    copy_VarCoords(y(0,:,:), prob)        ;copy coordinates to prob
    prob@_FillValue=-99.99

    ; set prob to missing where rc is missing (result otherwise 0)
    prob=where(ismissing(rc),rc,prob) 



; *******************************************************
; plot Overlay trend + probability

; wks = gsn_open_wks("X11" ,plotfile)
; wks = gsn_open_wks("eps" ,plotfile)
  wks_type = "png"
  wks_type@wkWidth = 800
  wks_type@wkHeight = 600
  wks = gsn_open_wks(wks_type,plotfile1)

;  gsn_define_colormap(wks,"hotcolr_19lev")

  greymap=(/"grey"/) ; colormap with grey
  gsn_merge_colormaps(wks,"hotcolr_19lev",greymap) ; merge two color maps 
;  gsn_draw_colormap(wks) ; (to see all available colors)

; **********************************************
; PLOT 1.1: Trend

  res                   = True
  res@gsnMaximize          = True            ; Maximize plot in frame
  res@gsnDraw           = False         ; Don't draw plots
  res@gsnFrame          = False         ; Don't advance frame
  res@cnFillOn          = True          ; Turn on contour fill
  res@cnFillMode       = "CellFill"     ; Raster Mode
  res@cnLinesOn         = False         ; Turn off contour lines

  res@lbLabelAutoStride = True          ; Control labelbar spacing

;; Robinson
;;  res@mpProjection          = "Robinson"
;;  ;res@mpOutlineOn            = True        ; turn on map outlines
;;  res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
;;  res@mpGridLineDashPattern = 2
;;  res@mpPerimOn         = False             ; turn off box around plot

  res@gsnAddCyclic          = False    ;for plotting subregions
  res@mpMinLatF            = -46       ; range to zoom in on
  res@mpMaxLatF            = -8.
  res@mpMinLonF            = 108.
  res@mpMaxLonF            = 158.

  res@gsnCenterString     = ""
  res@gsnLeftString     = ""
  res@gsnRightString     = ""
  res@cnInfoLabelOn     = False         ; turn off info label
  res@gsnTickMarksOn   = False          ; turn off all tick marks

  res@gsnSpreadColors       = True     ; use full range of color map
  res@gsnSpreadColorStart  = 2              ; Start at color index 1
  res@gsnSpreadColorEnd    = 21              ; Stop at 2nd color from end

  ;;; index specific settings.....
  if(index.eq."TXx".or.index.eq."TNn".or.index.eq."TNx".or.index.eq."TXn".or.index.eq."TX90p".or.index.eq."TN90p") then
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.1               ; set min contour level
   res@cnMaxLevelValF       =  0.1               ; set max contour level
   res@cnLevelSpacingF      =  0.02               ; set contour spacing
  end if
  if(index.eq."TX10p".or.index.eq."TN10p") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res@gsnSpreadColorStart  = 3   
   res@gsnSpreadColorEnd    = 22   
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;   res@cnMinLevelValF       = -0.15               ; set min contour level
;   res@cnMaxLevelValF       =  0.15               ; set max contour level
;   res@cnLevelSpacingF      =  0.03               ; set contour spacing
   res@cnMinLevelValF       = -0.4               ; set min contour level
   res@cnMaxLevelValF       =  0.4               ; set max contour level
   res@cnLevelSpacingF      =  0.1               ; set contour spacing
  end if
  if(index.eq."TX90p".or.index.eq."TN90p") then
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.4               ; set min contour level
   res@cnMaxLevelValF       =  0.4               ; set max contour level
   res@cnLevelSpacingF      =  0.1               ; set contour spacing
  end if
  if(index.eq."TX50p".or.index.eq."TN50p") then
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.4               ; set min contour level
   res@cnMaxLevelValF       =  0.4               ; set max contour level
   res@cnLevelSpacingF      =  0.1               ; set contour spacing
  end if
  if(index.eq."CDD".or.index.eq."GSL") then
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -1.               ; set min contour level
   res@cnMaxLevelValF       =  1.               ; set max contour level
   res@cnLevelSpacingF      =  0.2               ; set contour spacing
  end if
  if(index.eq."CWD") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res@gsnSpreadColorStart  = 3   
   res@gsnSpreadColorEnd    = 22              
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.2               ; set min contour level
   res@cnMaxLevelValF       =  0.2               ; set max contour level
   res@cnLevelSpacingF      =  0.05               ; set contour spacing
  end if
  if(index.eq."WSDI") then
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.3               ; set min contour level
   res@cnMaxLevelValF       =  0.3               ; set max contour level
   res@cnLevelSpacingF      =  0.1               ; set contour spacing
  end if
  if(index.eq."CSDI") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res@gsnSpreadColorStart  = 3   
   res@gsnSpreadColorEnd    = 22 
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.3               ; set min contour level
   res@cnMaxLevelValF       =  0.3               ; set max contour level
   res@cnLevelSpacingF      =  0.1               ; set contour spacing
  end if
  if(index.eq."DTR") then
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.06               ; set min contour level
   res@cnMaxLevelValF       =  0.06               ; set max contour level
   res@cnLevelSpacingF      =  0.01               ; set contour spacing
  end if
  if(index.eq."FD".or.index.eq."ID") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res@gsnSpreadColorStart  = 3   
   res@gsnSpreadColorEnd    = 22              
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.6               ; set min contour level
   res@cnMaxLevelValF       =  0.6               ; set max contour level
   res@cnLevelSpacingF      =  0.1               ; set contour spacing
  end if
  if(index.eq."TR".or.index.eq."SU") then
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.6               ; set min contour level
   res@cnMaxLevelValF       =  0.6               ; set max contour level
   res@cnLevelSpacingF      =  0.1               ; set contour spacing
  end if
  if(index.eq."R20mm".or.index.eq."R10mm".or.index.eq."Rnnmm".or.index.eq."SDII") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res@gsnSpreadColorStart  = 3   
   res@gsnSpreadColorEnd    = 22              
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -0.3               ; set min contour level
   res@cnMaxLevelValF       =  0.3               ; set max contour level
   res@cnLevelSpacingF      =  0.1               ; set contour spacing
  end if
  if(index.eq."Rx1day".or.index.eq."Rx5day") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res@gsnSpreadColorStart  = 3   
   res@gsnSpreadColorEnd    = 22              
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -2.               ; set min contour level
   res@cnMaxLevelValF       =  2.               ; set max contour level
   res@cnLevelSpacingF      =  0.5               ; set contour spacing
  end if
  if(index.eq."PRCPTOT") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res@gsnSpreadColorStart  = 3   
   res@gsnSpreadColorEnd    = 22              
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -10.               ; set min contour level
   res@cnMaxLevelValF       =  10.               ; set max contour level
   res@cnLevelSpacingF      =  2.               ; set contour spacing
  end if
  if(index.eq."R99p".or.index.eq."R95p") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res@gsnSpreadColorStart  = 3   
   res@gsnSpreadColorEnd    = 22              
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -3.               ; set min contour level
   res@cnMaxLevelValF       =  3.               ; set max contour level
   res@cnLevelSpacingF      =  1.               ; set contour spacing
  end if
  if(index.eq."R95pT") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res@gsnSpreadColorStart  = 3   
   res@gsnSpreadColorEnd    = 22              
   res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res@cnMinLevelValF       = -.2               ; set min contour level
   res@cnMaxLevelValF       =  .2               ; set max contour level
   res@cnLevelSpacingF      =  .05               ; set contour spacing
  end if

;;  res@cnLevelSelectionMode = "ExplicitLevels"    ; set explicit contour levels
;;  res@cnLevels    = (/ -0.1,-0.05,-0.02,-0.01,-0.005,0,0.005,0.01,0.02,0.05,0.1/)  ; set levels

   res@tiMainString        = "AUSDEX "+index+" Trend 1951-2011"

  plot1 = gsn_csm_contour_map_ce(wks,rc,res)  ; Create filled contours

; **********************************************
; PLOT 1.2: Significance

  sres = True                                 ; set up a second resource list
  sres@gsnDraw = False                        ; do not draw the plot
  sres@gsnFrame = False                       ; do not advance the frame

;  sres@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;  sres@cnMinLevelValF       = 0               ; set min contour level
;  sres@cnMaxLevelValF       = 0.1               ; set max contour level
  sres@cnInfoLabelOn     = False         ; turn off info label
  sres@cnLinesOn         = False         ; Turn off contour lines
  sres@cnLineLabelsOn      = False        ; no contour line labels



 plot2 = gsn_csm_contour(wks,prob,sres)  ; Create data2 contour plot

; ***********************************************
; 2b: SHADING
; use pattern fill #2 to fill all areas less than the first contour
; less than or equal to 0.1, and use pattern fill #17 to fill all areas greater
; than the first contour greater than or equal to 1.
;
  opt = True
  opt@gsnShadeFillType = "pattern"      ; pattern fill
  opt@gsnShadeHigh = 2                  ; use pattern #2
  opt@gsnShadeLow = 17               ; use pattern #6 (crossed lines)#17 (dots)
  plot2 = gsn_contour_shade(plot2,0.1,1.,opt)

; ***************************************************
; overlay...
  overlay(plot1,plot2)
  draw(plot1)
  frame(wks)


; ***************************************************
; Mean 1961-1990
  ;;TimAvg = dim_avg_n_Wrap(data1(12:41,:,:),0)
  TimAvg = dim_avg_n_Wrap(data1(11:40,:,:),0)

; ***************************************************
; Plot Mean

  wks2_type = "png"
  wks2_type@wkWidth = 800
  wks2_type@wkHeight = 600
  wks2 = gsn_open_wks(wks2_type,plotfile2)
  gsn_define_colormap(wks2,"rainbow+white+gray")

  res3                   = True
  res3@gsnMaximize          = True            ; Maximize plot in frame
  res3@cnFillOn          = True          ; Turn on contour fill
  res3@cnFillMode       = "CellFill"     ; Raster Mode
  res3@cnLinesOn         = False         ; Turn off contour lines

  res3@lbLabelAutoStride = True          ; Control labelbar spacing

;; Robinson
;;  res3@mpProjection          = "Robinson"
;;  res3@mpOutlineOn            = True        ; turn on map outlines
;;  res3@mpGridAndLimbOn   = True              ; turn on lat/lon lines
;;  res3@mpGridLineDashPattern = 2
;;  res3@mpPerimOn         = False             ; turn off box around plot

  res3@gsnAddCyclic          = False    ;for plotting subregions
  res3@mpMinLatF            = -46       ; range to zoom in on
  res3@mpMaxLatF            = -8.
  res3@mpMinLonF            = 108.
  res3@mpMaxLonF            = 158.

  res3@gsnCenterString     = ""
  res3@gsnLeftString     = ""
  res3@gsnRightString     = ""
  res3@cnInfoLabelOn     = False         ; turn off info label
  res3@gsnTickMarksOn   = False          ; turn off all tick marks

  res3@gsnSpreadColors       = True     ; use full range of color map
  res3@gsnSpreadColorStart  = 40            
  res3@gsnSpreadColorEnd    = 220             


  ;;; index specific settings (if not default is wanted)
  if(index.eq."TX90p".or.index.eq."TN90p") then
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  8               ; set min contour level
   res3@cnMaxLevelValF       =  12               ; set max contour level
   res3@cnLevelSpacingF      =  0.5               ; set contour spacing
  end if
  if(index.eq."TX10p".or.index.eq."TN10p") then
   ;;?!gsn_reverse_colormap(wks2)           ; Reverse the color map.
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  8               ; set min contour level
   res3@cnMaxLevelValF       =  12               ; set max contour level
   res3@cnLevelSpacingF      =  0.5               ; set contour spacing
  end if
  if(index.eq."TX50p".or.index.eq."TN50p") then
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  48               ; set min contour level
   res3@cnMaxLevelValF       =  52               ; set max contour level
   res3@cnLevelSpacingF      =  0.5               ; set contour spacing
  end if
  if(index.eq."Rx1day") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  100               ; set max contour level
   res3@cnLevelSpacingF      =  10               ; set contour spacing
  end if
  if(index.eq."Rx5day") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  200               ; set max contour level
   res3@cnLevelSpacingF      =  20               ; set contour spacing
  end if
  if(index.eq."CDD") then
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  200               ; set max contour level
   res3@cnLevelSpacingF      =  20               ; set contour spacing
  end if
  if(index.eq."CWD") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  10               ; set max contour level
   res3@cnLevelSpacingF      =  2               ; set contour spacing
  end if
  if(index.eq."CSDI".or.index.eq."WSDI") then
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  16               ; set max contour level
   res3@cnLevelSpacingF      =  2               ; set contour spacing
  end if
  if(index.eq."TR") then
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  300               ; set max contour level
   res3@cnLevelSpacingF      =  30               ; set contour spacing
  end if
  if(index.eq."PRCPTOT") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  1500               ; set max contour level
   res3@cnLevelSpacingF      =  150               ; set contour spacing
  end if
  if(index.eq."R10mm") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  80               ; set max contour level
   res3@cnLevelSpacingF      =  10               ; set contour spacing
  end if
  if(index.eq."R20mm") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  50               ; set max contour level
   res3@cnLevelSpacingF      =  5               ; set contour spacing
  end if
  if(index.eq."Rnnmm") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  36               ; set max contour level
   res3@cnLevelSpacingF      =  4               ; set contour spacing
  end if
  if(index.eq."R95p") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  600               ; set max contour level
   res3@cnLevelSpacingF      =  50               ; set contour spacing
  end if 
  if(index.eq."R95pT") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  16               ; set min contour level
   res3@cnMaxLevelValF       =  30               ; set max contour level
   res3@cnLevelSpacingF      =  2               ; set contour spacing
  end if 
  if(index.eq."R99p") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  300               ; set max contour level
   res3@cnLevelSpacingF      =  30               ; set contour spacing
  end if
  if(index.eq."SDII") then
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   res3@gsnSpreadColorStart  = 200            
   res3@gsnSpreadColorEnd    = 20 
   res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   res3@cnMinLevelValF       =  0               ; set min contour level
   res3@cnMaxLevelValF       =  25               ; set max contour level
   res3@cnLevelSpacingF      =  2.5               ; set contour spacing
  end if


   res3@tiMainString        = "AUSDEX "+index+" mean 1961-1990"

  plot3 = gsn_csm_contour_map_ce(wks2,TimAvg,res3)  ; Create filled contours

end do ;;end loop indices
end do ;;end loop gridsizes

end

; montage:

;files=`ls *TXx*par1_*trend*`
;
