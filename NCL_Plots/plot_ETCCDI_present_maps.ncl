
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

	index="TN90p"
	;pattern of the filenames
	patt_awap="CCRC_NARCliM_Sydney_AWAP_1990-2009_"
	patt_nnrp="CCRC_NARCliM_Sydney_NNRP_1990-2009_"
	patt_mk35="CCRC_NARCliM_Sydney_MK35_bc_1990-2009_"



	patt_imgout = "Sydney_ETCCDI_"
	;---SETUP
		limit_bot=0
		stride=1
	if (index.eq."TR") then
	limit_up=35
	spacing=5
	end if
	
	if (index.eq."SU") then
	limit_up=150
	spacing=5
	end if
	
	if (index.eq."FD") then
	limit_up=50
	spacing=5
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
	spacing=5
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
	spacing=5
	end if
	
	if (index.eq."TX90p") then
	limit_up=30
	spacing=2.5
	end if
	
	if (index.eq."GSL") then
	limit_bot=240	
	limit_up=365
	spacing=15
	end if
	
	if (index.eq."TXx") then
	limit_bot=25
	limit_up=45
	spacing=1
	end if
	
	if (index.eq."TNn") then
	limit_bot=-10
	limit_up=8
	spacing=1
	end if
	
	if (index.eq."TXn") then
	limit_bot=-2
	limit_up=16
	spacing=1
	end if
	
	if (index.eq."TNx") then
	limit_bot=10
	limit_up=30
	spacing=1
	end if
	
	if (index.eq."DTR") then
	limit_bot=6
	limit_up=16
	spacing=1
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
	spacing=5
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
	spacing=0.5
	stride=2
	end if
	
	if (index.eq."Rx1day") then
	limit_bot=0
	limit_up=200
	spacing=10
	stride=2
	end if
	
	if (index.eq."Rx5day") then
	limit_bot=0
	limit_up=400
	spacing=20
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
	spacing=1
	stride=2
	end if
	
	if (index.eq."R95p") then
	limit_bot=0
	limit_up=600
	spacing=25
	stride=2
	end if
	
	if (index.eq."R99p") then
	limit_bot=0
	limit_up=300
	spacing=25
	stride=2
	end if
	
	if (index.eq."PRCPTOT") then
	limit_bot=0
	limit_up=2500
	spacing=100
	stride=2
	end if
	
	
	zoom=True

	syear_pres = 1990
	eyear_pres = 2009
	nyears=eyear_pres-syear_pres+1
	;directories

	imgdiro = "~/Analyses/Sydney2km/ETCCDI/images/"
	diri_inv = "/home/z3393020/Analyses/Sydney2km/geo_em_files/" ; files to get the invariant fields.


	diri_awap = "/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"
	diri_nnrp = "/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"
	diri_mk35 = "/scratch/z3393020/ETCCDI/fclimdex_RCM/index/"

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

	;AWAP
	
	f = addfile(diri_awap+patt_awap+index+".nc","r")

	awap=f->Annual
	dsizes=dimsizes(awap)
	
	awap_aux=addfile("/scratch/z3393020/ETCCDI/fclimdex_RCM/tmin.SYD.1990-2009.nc","r")
	obslat=	awap_aux->lat
	obslon=	awap_aux->lon
	
	
	awap!1 = "lat"
	awap!2 = "lon"
	awap!0 = "time"
	awap&lat = obslat
	awap&lon = obslon
	awap@units=f@units
	print(awap@units)
	awap@_FillValue=awap@missing_value
	
		
	lud_awap = rcm2rgrid(lud_xlat,lud_xlon,lud_wrf,awap&lat,awap&lon,0)

	;WRF-NNRP
	
	f2 = addfile(diri_nnrp+patt_nnrp+index+".nc","r")

	wrf_nnrp=f2->Annual
	dsizes_f=dimsizes(wrf_nnrp)


	wrf_nnrp@lat2d = lat2d
	wrf_nnrp@lon2d = lon2d
	wrf_nnrp!1 = "lat"
	wrf_nnrp!2 = "lon"
	wrf_nnrp!0 = "time"
	wrf_nnrp@units=f2@units
	print(wrf_nnrp@units)

	if (index.eq."TN10p" .or. index.eq."TN50p" .or. index.eq."TN90p" .or. index.eq."TX10p" .or. index.eq."TX50p" .or. index.eq."TX90p") then
		wrf_nnrp@units="% of days"
	end if
		

	;WRF-MK35

	f3 = addfile(diri_mk35+patt_mk35+index+".nc","r")

	wrf_mk35=f3->Annual
	dsizes_f=dimsizes(wrf_mk35)


	wrf_mk35@lat2d = lat2d
	wrf_mk35@lon2d = lon2d
	wrf_mk35!1 = "lat"
	wrf_mk35!2 = "lon"
	wrf_mk35!0 = "time"
	wrf_mk35@units=f2@units
	print(wrf_mk35@units)

	if (index .eq. "GSL") then
		awap2=awap(0:nyears-2,:,:)
		wrf_nnrp2=wrf_nnrp(0:nyears-2,:,:)
		wrf_mk352=wrf_mk35(0:nyears-2,:,:)
		
		delete([/awap,wrf_nnrp,wrf_mk35/])
		
		awap=awap2
		wrf_nnrp=wrf_nnrp2
		wrf_mk35=wrf_mk352
	end if
		
	;##################################
	;#### END OF READING DATA #########
	;##################################


	printVarSummary(awap)


	print("calculate stats")
	awap_ann=dim_avg_n(awap,0)
	nnrp_ann=dim_avg_n(wrf_nnrp,0)
	mk35_ann=dim_avg_n(wrf_mk35,0)
	
	; if (index.eq."Rx1day" .or. index.eq."Rx5day")then
	; 	awap_ann=dim_max_n(awap,0)
	; 	nnrp_ann=dim_max_n(wrf_nnrp,0)
	; 	mk35_ann=dim_max_n(wrf_mk35,0)
	; end if
	
	
	
	awap_ann@_FillValue=awap@missing_value
	nnrp_ann@_FillValue=wrf_nnrp@missing_value
	mk35_ann@_FillValue=wrf_mk35@missing_value

	awap_ann!0 = "lat"
	awap_ann!1 = "lon"
	awap_ann&lat = awap&lat
	awap_ann&lon = awap&lon
	nnrp_ann@lat2d = lat2d
	nnrp_ann@lon2d = lon2d
	mk35_ann@lat2d = lat2d
	mk35_ann@lon2d = lon2d
	
		printVarSummary(mk35_ann)

	
	
	; awap_ann=mask(awap_ann,lud_awap,1)
	 	nnrp_ann=mask(nnrp_ann,lud_wrf,1)
	 	mk35_ann=mask(mk35_ann,lud_wrf,1)
	; limit_up=max((/max(awap_ann),max(nnrp_ann),max(mk35_ann)/))



	; #################################################
	; ########### PLOTTING ############################
	; #################################################
	;from example panel_15.ncl
	print("move to plotting")
	wks_type="pdf"
	wks_type@wkOrientation="landscape"
	wks = gsn_open_wks(wks_type,imgdiro+patt_imgout+index)
	gsn_define_colormap(wks, "WhViBlGrYeOrRe")
	plot=new(3,graphic)

	res=True
	res@cnFillOn		= True			;; turn on color
	res@cnLinesOn		= False
	res@gsnSpreadColors = True
	res@lbLabelBarOn	= False
	res@gsnDraw			= False
	res@gsnFrame		= False
	res@gsnAddCyclic 	= False

	res@cnLevelSelectionMode	= "ManualLevels"
	res@cnMinLevelValF			= limit_bot			;same manual levels
	res@cnMaxLevelValF			= limit_up			;as other plots
	res@cnLevelSpacingF			= spacing		;

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
	printVarSummary(awap_ann)



	plot(0)=gsn_csm_contour_map(wks,awap_ann,res)
	plot(1)=gsn_csm_contour_map(wks,nnrp_ann,res)
	plot(2)=gsn_csm_contour_map(wks,mk35_ann,res)



	;Panel first two plots


	pres1                     = True
	pres1@gsnPanelLabelBar    = True       ; common label bar
	pres1@gsnFrame            = False      ; don't advance frame yet
	pres1@lbLabelAutoStride   = False
	pres1@lbLabelStride = stride
	pres1@lbTitleString	= wrf_nnrp@units
	pres1@lbTitlePosition			= "Left"
	pres1@lbTitleDirection		= "Across"
	pres1@lbTitleFontHeightF		= 0.01
	pres1@gsnPanelFigureStrings			= (/"AWAP","NNRP","MK3.5_OEH"/)
	gsn_panel(wks,plot,(/1,3/),pres1)














	txres               = True                     ; text mods desired
	txres@txFontHeightF = 0.02                 ; font smaller. default big
	txres@txJust = "topCenter" ;; justification of text coordinates
	gsn_text_ndc(wks,f@long_name+" "+definition+" (1990-2009)",0.45,0.8,txres)
	frame(wks)

	; 
	; 
	; 	


	print("_______________END OF PROGRAM_____________________")
end
