
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
	patt_nnrp="CCRC_NARCliM_Sydney_1989-2009_"
	patt_oeh="CCRC_NARCliM_Sydney_2039-2059_"

	patt_imgout = "wrf2k_Sydney_ETCCDI_"
	;---SETUP

	zoom=True

	syear_pres = 1990
	eyear_pres = 2009

	;directories

	imgdiro = "~/Analyses/Sydney2km/ETCCDI/images/"
	diri_inv = "/home/z3393020/Analyses/Sydney2km/geo_em_files/" ; files to get the invariant fields.



	diri_nnrp = "/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"
	diri_std = "/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"

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

	f = addfile(diri_nnrp+patt_nnrp+index+".nc","r")

	wrf_nnrp=f->Annual
	dsizes=dimsizes(wrf_nnrp)


	wrf_nnrp@lat2d = lat2d
	wrf_nnrp@lon2d = lon2d
	wrf_nnrp!1 = "lat"
	wrf_nnrp!2 = "lon"
	wrf_nnrp!0 = "time"
	wrf_nnrp@units=f@units
	print(wrf_nnrp@units)

	; 	;#######################
	; 	;--READ ETCCDI STD
	; 	;#######################
	; 
	f2 = addfile(diri_std+patt_std+index+".nc","r")

	wrf_std=f2->Annual
	dsizes_f=dimsizes(wrf_std)


	wrf_std@lat2d = lat2d
	wrf_std@lon2d = lon2d
	wrf_std!1 = "lat"
	wrf_std!2 = "lon"
	wrf_std!0 = "time"
	wrf_std@units=f2@units
	print(wrf_std@units)


	; 	;#######################
	; 	;--READ ETCCDI OEH
	; 	;#######################
	; 
	f3 = addfile(diri_oeh+patt_nnrp+index+".nc","r")

	wrf_oeh=f3->Annual
	dsizes_f=dimsizes(wrf_oeh)


	wrf_oeh@lat2d = lat2d
	wrf_oeh@lon2d = lon2d
	wrf_oeh!1 = "lat"
	wrf_oeh!2 = "lon"
	wrf_oeh!0 = "time"
	wrf_oeh@units=f2@units
	print(wrf_oeh@units)

	;##################################
	;#### END OF READING DATA #########
	;##################################

	print("calculate stats")
	wrfann_nnrp=dim_avg_n(wrf_nnrp(1:,:,:),0)
	wrfann_std=dim_avg_n(wrf_std(1:dsizes_f(0)-1,:,:),0)
	wrfann_nnrp@lat2d = lat2d
	wrfann_nnrp@lon2d = lon2d
	wrfann_std@lat2d = lat2d
	wrfann_std@lon2d = lon2d
	wrfann_oeh@lat2d = lat2d
	wrfann_oeh@lon2d = lon2d
	
	wrfann_nnrp=mask(wrfann_nnrp,lud_wrf,1)
	wrfann_std=mask(wrfann_std,lud_wrf,1)
	wrfann_oeh=mask(wrfann_oeh,lud_wrf,1)
	limit_up=max((/max(wrfann_nnrp),max(wrfann_std),max(wrfann_oeh)/))



	; #################################################
	; ########### PLOTTING ############################
	; #################################################
	;from example panel_15.ncl
	print("move to plotting")
	wks_type="ps"
	wks_type@wkOrientation="landscape"
	wks = gsn_open_wks(wks_type,imgdiro+patt_imgout+index)
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

	plot(0)=gsn_csm_contour_map(wks,wrfann_nnrp(:,:),res)
	plot(1)=gsn_csm_contour_map(wks,wrfann_std(:,:),res)
	plot(2)=gsn_csm_contour_map(wks,wrfann_oeh(:,:),res)



	;Panel first two plots


	pres1                     = True
	pres1@gsnPanelLabelBar    = True       ; common label bar
	pres1@gsnFrame            = False      ; don't advance frame yet
	pres1@lbLabelAutoStride   = True 
	pres1@lbTitleString	= wrf_pres@units
	pres1@lbTitlePosition			= "Left"
	pres1@lbTitleDirection		= "Across"
	pres1@lbTitleFontHeightF		= 0.01
	pres1@gsnPanelFigureStrings			= (/"NNRP","MK3.5_Std","MK3.5_OEH"/)
	gsn_panel(wks,plot,(/1,3/),pres1)














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
