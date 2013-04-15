
;+==============================
; Daniel Argueso
; University of New South Wales
; 09 April 2013
; Last version: 14 Jan 2013
; email: d.argueso@unsw.edu.au
;==============================

; Script based on Markus': trend+mean_AUSDEX_BAK.ncl
; Purpose: Extracts information from the fclimdex_RCM output files and plot maps of index trends.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;************************************************
undef("calc_trend")
function calc_trend (data,nyrs,index)
local x,tval,nxym,numyrs,b
begin
    x = ispan(0,nyrs-1,1)
    ;count the number of non-missing (annual) values at each location
    numyrs = dim_num_n(.not.ismissing(data),0)
    tval = new(dimsizes(data(0,:,:)),float,1e20)
    nxy = new(dimsizes(data(0,:,:)),integer)
    rc = regcoef (x,data(lat|:,lon|:,time|:),tval,nxy)
    rc@long_name = "linear Trend "+index
    rc=where(numyrs.lt.30,rc@_FillValue,rc)
    
    b = tval
    b = 0.5
    prob  = betainc(nxy/(nxy+tval^2),nxy/2.0,b)
    copy_VarCoords(data(0,:,:),prob)
    prob@_FillValue=1e20
    ; set prob to missing where rc is missing (result otherwise 0)
    prob=where(ismissing(rc),rc,prob) 
    return(/rc,prob/)
end

begin

;indices=("TN10p") 
;indices=(/"CDD","CWD","GSL","PRCPTOT","R20mm","R99p","Rx1day","SDII","TN50p","TNn","TR","TX50p","TXn","WSDI","CSDI","DTR","FD","ID","R10mm","R95p","R95pT","Rnnmm","Rx5day","SU","TN10p","TN90p","TNx","TX10p","TX90p","TXx"/)
;indices=(/"CDD","CWD","GSL","PRCPTOT","R20mm","R99p","Rx1day","SDII","TNn","TR","TXn","WSDI","CSDI","DTR","FD","ID","R10mm","R95p","Rx5day","SU","TN10p","TN90p","TNx","TX10p","TX90p","TXx"/)
indices=(/"Rx5day","R99p","R95p","SDII","TN10p","TX10p","TN90p","TX90p"/)
geof=addfile("/srv/ccrc/data13/z3393020/Analyses/Sydney2km/geo_em_files/geo_em.d01.narclim.AUS_ONLY.nc","r")
landmask=geof->LANDMASK(0,:,:)
lat2d=geof->XLAT_M(0,:,:)
lon2d=geof->XLONG_M(0,:,:)
nlat  = dimsizes(lat2d(:,0))
nlon  = dimsizes(lon2d(0,:))

;adding AUS map
;fshape="/srv/ccrc/data13/z3393020/Analyses/aus_bounds/AUS_adm0.shp"


do ii=0,dimsizes(indices)-1
  index=indices(ii)
  
  plotfile="/srv/ccrc/data13/z3393020/Analyses/NARCliM/ETCCDI/images/NARCliM_"+index+"map"
  wks_type = "png"
  wks_type@wkWidth = 1600
  wks_type@wkHeight = 1200
  wks = gsn_open_wks(wks_type,plotfile)
  
  
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
  

  
  data1 = mask (data1, conform(data1, landmask, (/1,2/)), 1)  
  data2 = mask (data2, conform(data2, landmask, (/1,2/)), 1)  
  data3 = mask (data3, conform(data3, landmask, (/1,2/)), 1)
  data4 = mask (data4, conform(data4, maskWA  , (/1,2/)), 1)

  ; *********************************************
   
     ; set data to missing on gridboxes where <40 years of data
        ; nomiss_dat1=dim_num_n(.not.ismissing(data1),0)  ; count missing values at each location
        ; data1 = mask( data1, (nomiss_dat1.ge.40), True)  ; mask
        ; nomiss_dat2=dim_num_n(.not.ismissing(data2),0)  ; count missing values at each location
        ; data2 = mask( data2, (nomiss_dat2.ge.40), True)  ; mask
        ; nomiss_dat3=dim_num_n(.not.ismissing(data3),0)  ; count missing values at each location
        ; data3 = mask( data3, (nomiss_dat3.ge.40), True)  ; mask
        ; nomiss_dat4=dim_num_n(.not.ismissing(data4),0)  ; count missing values at each location
        ; data4 = mask( data4, (nomiss_dat4.ge.40), True)  ; mask
        ; nomiss_dat5=dim_num_n(.not.ismissing(data5),0)  ; count missing values at each location
        ; data5 = mask( data5, (nomiss_dat5.ge.40), True)  ; mask
        ; nomiss_dat6=dim_num_n(.not.ismissing(data6),0)  ; count missing values at each location
        ; data6 = mask( data6, (nomiss_dat6.ge.40), True)  ; mask
   
   
   
       print ("plotting index "+index+"...")
       greymap=(/"gray"/) ; colormap with grey
       gsn_merge_colormaps(wks,"hotcolr_19lev",greymap) ;
       
       rc_prob1=calc_trend(data1,nyrs,index)
       rc_prob1@lon2d=lon2d
       rc_prob1@lat2d=lat2d
       rc_prob1@_FillValue=1e20
       
       rc_prob2=calc_trend(data2,nyrs,index)
       rc_prob2@lon2d=lon2d
       rc_prob2@lat2d=lat2d
       rc_prob2@_FillValue=1e20
       
       rc_prob3=calc_trend(data3,nyrs,index)
       rc_prob3@lon2d=lon2d
       rc_prob3@lat2d=lat2d
       rc_prob3@_FillValue=1e20
       
       rc_prob4=calc_trend(data4,nyrs,index)
       copy_VarCoords(data4(0,:,:),rc_prob4(0,:,:))
       copy_VarCoords(data4(0,:,:),rc_prob4(1,:,:))
       rc_prob4@_FillValue=1e20
       
       rc_prob5=calc_trend(data5,nyrs,index)
       copy_VarCoords(data5(0,:,:),rc_prob5(0,:,:))
       copy_VarCoords(data5(0,:,:),rc_prob5(1,:,:))
       rc_prob5@_FillValue=1e20
   
       rc_prob6=calc_trend(data6,nyrs,index)
       copy_VarCoords(data6(0,:,:),rc_prob6(0,:,:))
       copy_VarCoords(data6(0,:,:),rc_prob6(1,:,:))
       rc_prob6@_FillValue=1e20
       
   
   
   
   
   
   
   
       res = True
       
       ;Defining the domain
       
       res@mpCenterLonF           = geof@CEN_LON
     res@mpCenterLatF           = geof@POLE_LAT-90
       res@mpLimitMode       = "Corners"
         res@mpLeftCornerLatF       = lat2d(25,125)
         res@mpLeftCornerLonF       = lon2d(25,125)
         res@mpRightCornerLatF      = lat2d(102,30)
         res@mpRightCornerLonF      = lon2d(102,30)+360
   
       ;     res@mpLandFillColor = -1
       ; res@mpOceanFillColor = -1
       ; res@mpFillDrawOrder = "Draw"
     
     ;res@mpOutlineBoundarySets  = "NoBoundaries"
   
       res@gsnMaximize          = True            ; Maximize plot in frame
       res@gsnDraw           = False         ; Don't draw plots
       res@gsnFrame          = False         ; Don't advance frame
       res@gsnAddCyclic          = False    ;for plotting subregions
       res@gsnCenterString     = ""
       res@gsnLeftString     = ""
       res@gsnRightString     = ""
       res@gsnTickMarksOn   = False          ; turn off all tick marks
       
       res@cnFillOn          = True          ; Turn on contour fill
       res@cnFillMode       = "CellFill"     ; Raster Mode
       res@cnLinesOn         = False         ; Turn off contour lines
       res@cnInfoLabelOn     = False         ; turn off info label
       res@cnLineLabelsOn = False
       res@lbLabelBarOn              = False
       ;res@lbLabelAutoStride = True          ; Control labelbar spacing
   
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
      res@cnMinLevelValF       = -0.3               ; set min contour level
      res@cnMaxLevelValF       =  0.3               ; set max contour level
      res@cnLevelSpacingF      =  0.05               ; set contour spacing
     end if
     if(index.eq."TX90p".or.index.eq."TN90p") then
      res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
      res@cnMinLevelValF       = -0.3               ; set min contour level
      res@cnMaxLevelValF       =  0.3               ; set max contour level
      res@cnLevelSpacingF      =  0.05               ; set contour spacing
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
      res@cnLevelSpacingF      =  0.05               ; set contour spacing
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
      res@cnLevelSpacingF      =  0.5               ; set contour spacing
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
   
   
   
      res@tiMainString        = "NARCliM "+index+" Trend 1951-2011"
   
      pres             = True
      pres@gsLineColor = "black"
      pres@gsLineThicknessF = 0.4
     
     plot1=new(6,graphic)
     
     res@tiMainString        = "NARCliM R1"
     plot1(0) = gsn_csm_contour_map(wks,rc_prob1(0,:,:),res)  ; Create filled contours
      res@tiMainString        = "NARCliM R2"
     plot1(1) = gsn_csm_contour_map(wks,rc_prob2(0,:,:),res)  ; Create filled contours
      res@tiMainString        = "NARCliM R3"
     plot1(2) = gsn_csm_contour_map(wks,rc_prob3(0,:,:),res)  ; Create filled contours
      res@tiMainString        = "AWAP"
     plot1(3) = gsn_csm_contour_map(wks,rc_prob4(0,:,:),res)  ; Create filled contours
      res@tiMainString        = "HadEX2"
     plot1(4) = gsn_csm_contour_map(wks,rc_prob5(0,:,:),res)  ; Create filled contours
      res@tiMainString        = "NNRP 1"
     plot1(5) = gsn_csm_contour_map(wks,rc_prob6(0,:,:),res)  ; Create filled contours
   ; **********************************************
   ; PLOT 1.2: Significance
   
     sres = True                                 ; set up a second resource list
     sres@gsnDraw = False                        ; do not draw the plot
     sres@gsnFrame = False                       ; do not advance the frame
   
   
     sres@cnInfoLabelOn     = False         ; turn off info label
     sres@cnLinesOn         = False         ; Turn off contour lines
     sres@cnLineLabelsOn      = False        ; no contour line labels
   
   
    plot2=new(6,graphic)
    plot2(0) = gsn_csm_contour(wks,rc_prob1(1,:,:),sres)  ; Create data2 contour plot
    plot2(1) = gsn_csm_contour(wks,rc_prob2(1,:,:),sres)  ; Create data2 contour plot
    plot2(2) = gsn_csm_contour(wks,rc_prob3(1,:,:),sres)  ; Create data2 contour plot
    plot2(3) = gsn_csm_contour(wks,rc_prob4(1,:,:),sres)  ; Create data2 contour plot
    plot2(4) = gsn_csm_contour(wks,rc_prob5(1,:,:),sres)  ; Create data2 contour plot
    plot2(5) = gsn_csm_contour(wks,rc_prob6(1,:,:),sres)  ; Create data2 contour plot
   ; ***********************************************
   ; 2b: SHADING
   ; use pattern fill #2 to fill all areas less than the first contour
   ; less than or equal to 0.1, and use pattern fill #17 to fill all areas greater
   ; than the first contour greater than or equal to 1.
   ;
     opt = True
     opt@gsnShadeFillType = "pattern"      ; pattern fill
     opt@gsnShadeHigh = 2                  ; use pattern #2
     opt@gsnShadeLow = 6              ; use pattern #6 (crossed lines)#17 (dots)
     
     
     do i=0,dimsizes(plot1)-1
     plot2(i) = gsn_contour_shade(plot2(i),0.1,1.,opt)
     end do
   
   ; ***************************************************
   ; overlay...
   do i=0,dimsizes(plot1)-1
     overlay(plot1(i),plot2(i))
     end do
   
     
     pres=True
     pres@txString = index+" trend 1950-2009"
     pres@gsnPanelLabelBar           = True
     pres@lbAutoManage           = True
     pres@lbTitleString          = "units/year"
     pres@lbTitlePosition            = "Left"
     pres@lbTitleDirection       = "Across"
     pres@lbTitleFontHeightF     = 0.012
     pres@lbLabelFontHeightF     = 0.012
   
     
     gsn_panel(wks,plot1,(/2,3/),pres)
     
    
   

end do ;;end loop indices

end

