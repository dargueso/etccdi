
;Name:
;	plot_ETCCDI_maps.ncl
;
;Purpose:
; This script reads files from fclimdex_RCM and plots the values of the indices for different months.
;
;Author:
;	Daniel Argüeso @ UNSW 
;
;Creation:
;	10/08/2012
;
;Last Modified:


begin

	print("____________START OF PROGRAM________________________")

	;---LIBRARIES
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"




	;---WRF DOMAINS

	index="TR"
	;pattern of the filenames
	patt_pres="CCRC_NARCliM_Sydney_1989-2009_"
	patt_fut="CCRC_NARCliM_Sydney_2039-2059_"

	patt_imgout = "wrf2k_Sydney_ETCCDI_"
	;---SETUP

	zoom=True

	syear_pres = 1990
	eyear_pres = 2009

	syear_fut = 2040
	eyear_fut = 2059
	;directories

	imgdiro = "~/Analyses/Sydney2km/ETCCDI/images/"
	diri_inv = "/home/z3393020/Analyses/Sydney2km/geo_em_files/" ; files to get the invariant fields.



	diri_pres = "/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"
	diri_fut = "/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"

	;---SETUP LANDUSE MASK
	;use landuse data to mask out ocean from plots
	luf = addfile(diri_inv+"geo_em.d03.nc","r")
	lud_wrf =luf->LANDMASK(0,:,:)
	lud_xlat=luf->XLAT_M(0,:,:)
	lud_xlon=luf->XLONG_M(0,:,:)
	lat2d = luf->XLAT_M(0,:,:)
	lon2d = luf->XLONG_M(0,:,:)







	;#######################
	;--READ ETCCDI PRESENT
	;#######################

	f = addfile(diri_pres+patt_pres+index+".nc","r")

	wrf_pres=f->Annual
	dsizes=dimsizes(wrf_pres)


	wrf_pres@lat2d = lat2d
	wrf_pres@lon2d = lon2d
	wrf_pres!1 = "lat"
	wrf_pres!2 = "lon"
	wrf_pres!0 = "time"
	wrf_pres@units=f@units
	print(wrf_pres@units)

	; 	;#######################
	; 	;--READ ETCCDI FUTURE
	; 	;#######################
	; 
	f2 = addfile(diri_pres+patt_pres+index+".nc","r")

	wrf_fut=f2->Annual
	dsizes_f=dimsizes(wrf_fut)


	wrf_fut@lat2d = lat2d
	wrf_fut@lon2d = lon2d
	wrf_fut!1 = "lat"
	wrf_fut!2 = "lon"
	wrf_fut!0 = "time"
	wrf_fut@units=f2@units
	print(wrf_fut@units)

	;##################################
	;#### END OF READING DATA #########
	;##################################

	print("calculate stats")
	wrfann_pres=dim_avg_n(wrf_pres(1:,:,:),0)
	wrfann_fut=dim_avg_n(wrf_fut(1:dsizes_f(0)-1,:,:),0)
	wrfann_pres@lat2d = lat2d
	wrfann_pres@lon2d = lon2d
	wrfann_fut@lat2d = lat2d
	wrfann_fut@lon2d = lon2d
	
	
	wrfann_pres=mask(wrfann_pres,lud_wrf,1)
	wrfann_fut=mask(wrfann_fut,lud_wrf,1)
	limit_up=max((/max(wrfann_pres),max(wrfann_fut)/))


	;calculate changes

	ind_changes=wrfann_fut-wrfann_pres
	ind_changes@lat2d = lat2d
	ind_changes@lon2d = lon2d


	; #################################################
	; ########### PLOTTING ############################
	; #################################################
	;from example panel_15.ncl
	print("move to plotting")
	wks_type="ps"
	wks_type@wkOrientation="landscape"
	wks = gsn_open_wks(wks_type,imgdiro+patt_imgout+index+"_changes")
	gsn_define_colormap(wks, "WhViBlGrYeOrRe")
	plot=new(2,graphic)

	res=True
	res@cnFillOn		= True			;; turn on color
	res@cnLinesOn		= False
	res@gsnSpreadColors = True
	res@lbLabelBarOn	= False
	res@gsnDraw			= False
	res@gsnFrame		= False
	res@gsnAddCyclic 	= False

	res@cnLevelSelectionMode	= "ManualLevels"
	res@cnMinLevelValF			= 0				;same manual levels
	res@cnMaxLevelValF			= limit_up			;as other plots
	res@cnLevelSpacingF			= 5.0			;

	;###### MY RESOURCES #####
	
	res@cnLineLabelsOn		= False
	res@cnInfoLabelOn		= False
	
	res@gsnSpreadColorStart	= 10
	res@gsnSpreadColorEnd	= -1
	res@cnFillMode			= "CellFill"		
	res@cnRasterSmoothingOn 	= False

	res@mpDataBaseVersion	= "HighRes"
	res@mpLimitMode		= "Corners"
	res@pmTickMarkDisplayMode 	= "Never"

	res@mpOutlineDrawOrder		= "PostDraw"  ; force map tp be drawn 1st 
	res@mpGridLineDashPattern		= 2           ; lat/lon lines as dashed
	res@mpPerimOn				= True
	res@mpPerimDrawOrder			= "PostDraw"
	res@mpOutlineOn				= True
	res@mpOutlineBoundarySets		= "National"
	res@mpGeophysicalLineThicknessF = 1
	res@mpDataBaseVersion			= "HighRes"
	res@mpLabelDrawOrder			= "PostDraw"

	res@mpLandFillColor = -1
	res@mpOceanFillColor = 0
	res@mpFillDrawOrder = "PostDraw"

	res@mpLeftCornerLatF	= min(lat2d)	;
	res@mpLeftCornerLonF	= min(lon2d)	;extend beyond strict
	res@mpRightCornerLatF	= max(lat2d)	;grt syd domain
	res@mpRightCornerLonF	= max(lon2d)	;

	if (zoom) then
		;zoom to Sydney urban region
		res@mpLeftCornerLatF     = -32.94     ;
		res@mpLeftCornerLonF     = 149.90          ;zoom
		res@mpRightCornerLatF     = -34.54          ;to urban
		res@mpRightCornerLonF     = 151.5          ;


	end if

	res@mpGridAndLimbOn			= True		;reference lat/
	res@mpGridSpacingF			= 5.			;lon grids

	res@gsnRightString			= ""			;suppress
	res@gsnCenterString			= ""			;labelling
	res@gsnLeftString			= "" 		;for each panel
	res@lbLabelBarOn				= False
	res@mpProjection				= "LambertConformal"
	res@mpLambertParallel1F		= -27.5
	res@mpLambertParallel2F		= -35.
	res@mpLambertMeridianF		= 146.1

	;###### END OF MY RESOURCES #####

	plot(0)=gsn_csm_contour_map(wks,wrfann_pres(:,:),res)
	plot(1)=gsn_csm_contour_map(wks,wrfann_fut(:,:),res)



	;Panel first two plots


	pres1                     = True
	pres1@gsnPanelLabelBar    = True       ; common label bar
	pres1@gsnFrame            = False      ; don't advance frame yet
	pres1@lbLabelAutoStride   = True 
	pres1@lbTitleString	= wrf_pres@units
	pres1@lbTitlePosition			= "Left"
	pres1@lbTitleDirection		= "Across"
	pres1@lbTitleFontHeightF		= 0.014
	
	pres1@vpWidthF=0.55
	pres1@gsnPanelRight =0.6
	pres1@gsnPanelFigureStrings			= (/"Present","Future"/)
	gsn_panel(wks,plot,(/1,2/),pres1)




	;***************************************
	; create third individual plots
	;***************************************

	res@cnMinLevelValF		= 0			;same manual levels
	res@cnMaxLevelValF		= 20.0			;as other plots
	res@cnLevelSpacingF	= 2			;
	plot2 = gsn_csm_contour_map(wks,ind_changes(:,:),res)


	pres2                     = True
	pres2@gsnPanelLabelBar    = True       ; common label bar
	pres2@gsnFrame            = False      ; don't advance frame yet
	pres2@lbLabelAutoStride   = True 
	pres2@lbTitleString	= wrf_pres@units
	pres2@lbTitlePosition			= "Left"
	pres2@lbTitleDirection		= "Across"
	pres2@lbTitleFontHeightF		= 0.014

	pres2@gsnPanelLeft  = 0.6
	pres2@gsnPanelRight = 0.9
	pres2@gsnPanelFigureStrings			= (/"Changes"/)
	
	gsn_panel(wks,plot2,(/1,1/),pres2)









	txres               = True                     ; text mods desired
	txres@txFontHeightF = 0.02                 ; font smaller. default big
	txres@txJust = "topCenter" ;; justification of text coordinates
	gsn_text_ndc(wks,"Changes in "+f@long_name+" (MK3.5)",0.45,0.8,txres)
	frame(wks)

	; 
	; 
	; 	


	print("_______________END OF PROGRAM_____________________")
end
