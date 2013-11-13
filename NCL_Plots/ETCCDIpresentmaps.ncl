
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

begin

;indices=("TN10p") 
;indices=(/"CDD","CWD","GSL","PRCPTOT","R20mm","R99p","Rx1day","SDII","TN50p","TNn","TR","TX50p","TXn","WSDI","CSDI","DTR","FD","ID","R10mm","R95p","R95pT","Rnnmm","Rx5day","SU","TN10p","TN90p","TNx","TX10p","TX90p","TXx"/)
indices=(/"CDD","CWD","GSL","PRCPTOT","R20mm","R99p","Rx1day","SDII","TNn","TR","TXn","WSDI","CSDI","DTR","FD","ID","R10mm","R95p","Rx5day","SU","TN10p","TN90p","TNx","TX10p","TX90p","TXx"/)
;indices=(/"Rx5day","R99p","R95p","SDII","TN10p","TX10p","TN90p","TX90p"/)
geof=addfile("/srv/ccrc/data13/z3393020/Analyses/share/geo_em_files/geo_em.d01.narclim.AUS_ONLY.nc","r")
landmask=geof->LANDMASK(0,:,:)
lat2d=geof->XLAT_M(0,:,:)
lon2d=geof->XLONG_M(0,:,:)
nlat  = dimsizes(lat2d(:,0))
nlon  = dimsizes(lon2d(0,:))

;adding AUS map
;fshape="/srv/ccrc/data13/z3393020/Analyses/aus_bounds/AUS_adm0.shp"


do ii=0,dimsizes(indices)-1
  index=indices(ii)
  
  plotfile="/srv/ccrc/data13/z3393020/Analyses/NARCliM/ETCCDI/images/NARCliM_GHCNDEX_mean_"+index+"map"
  wks_type = "png"
  wks_type@wkWidth = 1600
  wks_type@wkHeight = 1200
  wks = gsn_open_wks(wks_type,plotfile)
  
  
  DATADIR1="/home/z3393020//Analyses/NARCliM/ETCCDI/indices/R1/NNRP/Bias_corrected/"
  INDAT1=DATADIR1+"CCRC_NARCliM_NNRP_R1_bc_1950-2009_"+index+".nc"
  DATADIR2="/home/z3393020//Analyses/NARCliM/ETCCDI/indices/R2/NNRP/Bias_corrected/"
  INDAT2=DATADIR2+"CCRC_NARCliM_NNRP_R2_bc_1950-2009_"+index+".nc"
  DATADIR3="/home/z3393020//Analyses/NARCliM/ETCCDI/indices/R3/NNRP/Bias_corrected/"
  INDAT3=DATADIR3+"CCRC_NARCliM_NNRP_R3_bc_1950-2009_"+index+".nc"
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
  ;DATADIR5="/srv/ccrc/data03/z3356123/HadEX2/HadEX2_3.75x2.5deg/"
  ;INDAT5=DATADIR5+"H2_"+index+"_1901-2010_RegularGrid_global_3.75x2.5deg_LSmask.nc"

  ; ;GHCNDEX
  DATADIR5="/srv/ccrc/data07/z3356123/GHCNDEX_2012/ghcndex_Apr2013/gridout/GHCNDEX_2.5x2.5deg/"
  
  ;  DATADIR5="/srv/ccrc/data07/z3356123/GHCNDEX_2012/ghcndex/gridout/GHCNDEX_2.5x2.5deg/"
    INDAT5=DATADIR5+"GHCND_"+index+"_1951-2013_RegularGrid_global_2.5x2.5deg_LSmask.nc"
  
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
  sy_GHCN=syear-1951
  ;sy_GHCN=syear-1901 ;HadEX2
  
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

      
  
       print ("plotting index "+index+"...")
       greymap=(/"gray"/) ; colormap with grey
       ;gsn_merge_colormaps(wks,"hotcolr_19lev",greymap) ;
       gsn_merge_colormaps(wks, "WhViBlGrYeOrRe",greymap)
       
       rc_prob1=dim_avg_n(data1,0)
       rc_prob1@lon2d=lon2d
       rc_prob1@lat2d=lat2d
       rc_prob1@_FillValue=1e20
       
       rc_prob2=dim_avg_n(data2,0)
       rc_prob2@lon2d=lon2d
       rc_prob2@lat2d=lat2d
       rc_prob2@_FillValue=1e20
       
       rc_prob3=dim_avg_n(data3,0)
       rc_prob3@lon2d=lon2d
       rc_prob3@lat2d=lat2d
       rc_prob3@_FillValue=1e20
       
       rc_prob4=dim_avg_n(data4,0)
       print (dimsizes(data4))
       print (dimsizes(rc_prob4))
       copy_VarCoords(data4(0,:,:),rc_prob4(:,:))
       rc_prob4@_FillValue=1e20
       
       rc_prob5=dim_avg_n(data5,0)
       copy_VarCoords(data5(0,:,:),rc_prob5(:,:))
       rc_prob5@_FillValue=1e20
   
       rc_prob6=dim_avg_n(data6,0)
       copy_VarCoords(data6(0,:,:),rc_prob6(:,:))
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
   
   
   
      ;---SETUP
    		limit_bot=0
    		stride=1
    	if (index.eq."TR") then
    	limit_up=35
    	spacing=5.
    	end if

    	if (index.eq."SU") then
    	limit_up=300
    	spacing=5.
    	end if

    	if (index.eq."FD") then
    	limit_up=50
    	spacing=5.
    	end if

    	if (index.eq."ID") then
    	limit_up=10
    	spacing=0.5
    	end if

    	if (index.eq."TN10p") then
    	limit_up=30
    	spacing=2.5
    	end if

    	if (index.eq."TN50p") then
    	limit_up=100
    	spacing=5.
    	end if

    	if (index.eq."TN90p") then
    	limit_up=30
    	spacing=2.5
    	end if

    	if (index.eq."TX10p") then
    	limit_up=30
    	spacing=2.5
    	end if

    	if (index.eq."TX50p") then
    	limit_up=100
    	spacing=5.
    	end if

    	if (index.eq."TX90p") then
    	limit_up=30
    	spacing=2.5
    	end if

    	if (index.eq."GSL") then
    	limit_bot=240	
    	limit_up=365
    	spacing=15.
    	end if

    	if (index.eq."TXx") then
    	limit_bot=25
    	limit_up=45
    	spacing=1.
    	end if

    	if (index.eq."TNn") then
    	limit_bot=-10
    	limit_up=8
    	spacing=1.
    	end if

    	if (index.eq."TXn") then
    	limit_bot=-2
    	limit_up=16
    	spacing=1.
    	end if

    	if (index.eq."TNx") then
    	limit_bot=10
    	limit_up=30
    	spacing=1.
    	end if

    	if (index.eq."DTR") then
    	limit_bot=6
    	limit_up=16
    	spacing=1.
    	end if

    	if (index.eq."WSDI") then
    	limit_bot=0
    	limit_up=4
    	spacing=0.25
    	stride=2
    	end if

    	if (index.eq."CSDI") then
    	limit_bot=0
    	limit_up=3
    	spacing=0.25
    	stride=2
    	end if

    	;PRECIP EXTREMES

    	if (index.eq."R10mm") then
    	limit_bot=0
    	limit_up=60
    	spacing=5.
    	end if

    	if (index.eq."R20mm") then
    	limit_bot=0
    	limit_up=30
    	spacing=2.5
    	stride=2
    	end if

    	if (index.eq."SDII") then
    	limit_bot=0
    	limit_up=15
    	spacing=1.
    	stride=2
    	end if

    	if (index.eq."Rx1day") then
    	limit_bot=0
    	limit_up=200
    	spacing=10.
    	stride=2
    	end if

    	if (index.eq."Rx5day") then
    	limit_bot=0
    	limit_up=400
    	spacing=20.
    	stride=2
    	end if

    	if (index.eq."CDD") then
    	limit_bot=0
    	limit_up=60
    	spacing=2.5
    	stride=2
    	end if

    	if (index.eq."CWD") then
    	limit_bot=0
    	limit_up=20
    	spacing=1.
    	stride=2
    	end if

    	if (index.eq."R95p") then
    	limit_bot=0
    	limit_up=600
    	spacing=25.
    	stride=2
    	end if

    	if (index.eq."R99p") then
    	limit_bot=0
    	limit_up=300
    	spacing=25.
    	stride=2
    	end if

    	if (index.eq."PRCPTOT") then
    	limit_bot=0
    	limit_up=2500
    	spacing=100.
    	stride=2
    	end if
    	
    	res@gsnSpreadColors = True
      res@gsnSpreadColorStart  = 0              ; Start at color index 1
      res@gsnSpreadColorEnd    = 102              ; Stop at 2nd color from end
      res@mpLandFillColor = 103
    	res@mpOceanFillColor = 0
    	res@mpFillDrawOrder = "PreDraw"
    	
      res@cnLevelSelectionMode	= "ManualLevels"
    	res@cnMinLevelValF			= limit_bot			;same manual levels
    	res@cnMaxLevelValF			= limit_up			;as other plots
    	res@cnLevelSpacingF			= spacing		;
   
   
   
      res@tiMainString        = "NARCliM "+index+" Trend 1951-2011"
   
      pres             = True
      pres@gsLineColor = "black"
      pres@gsLineThicknessF = 0.4
     
     plot1=new(6,graphic)
     
     res@tiMainString        = "NARCliM R1"
     plot1(0) = gsn_csm_contour_map(wks,rc_prob1(:,:),res)  ; Create filled contours
      res@tiMainString        = "NARCliM R2"
     plot1(1) = gsn_csm_contour_map(wks,rc_prob2(:,:),res)  ; Create filled contours
      res@tiMainString        = "NARCliM R3"
     plot1(2) = gsn_csm_contour_map(wks,rc_prob3(:,:),res)  ; Create filled contours
      res@tiMainString        = "AWAP"
     plot1(3) = gsn_csm_contour_map(wks,rc_prob4(:,:),res)  ; Create filled contours
      res@tiMainString        = "GHCNDEX"
     plot1(4) = gsn_csm_contour_map(wks,rc_prob5(:,:),res)  ; Create filled contours
      res@tiMainString        = "NNRP 1"
     plot1(5) = gsn_csm_contour_map(wks,rc_prob6(:,:),res)  ; Create filled contours
   ; **********************************************

   
     
     pres=True
     pres@txString = index+"mean annual values 1950-2009"
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

