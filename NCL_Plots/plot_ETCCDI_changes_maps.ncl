
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
; 31/10/2012 Fixed no ZOOM option


begin

	print("____________START OF PROGRAM________________________")

	;---LIBRARIES
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


	syear=1990
	eyear=2009
	nyears=eyear-syear+1

	;---WRF DOMAINS

	index="TX90p"
	
	
		limit_bot=0.
		stride=2.
		spacing_c=0.5
		stride_c=2.	
		spacing=5.
		definition=" "

	if (index.eq."TR") then
		limit_up=60.
		spacing=5
		stride=2
		spacing_c=5
		stride_c=2
		
	definition="(Tmin>20~S~o~N~C)"
	end if
	
	if (index.eq."SU") then
		limit_up=180.
		spacing=15.
		stride=2.
			spacing_c=5
			stride_c=2
	definition="(Tmax>25~S~o~N~C)"
	end if
	
	if (index.eq."FD") then
		limit_up=70.
		spacing=2.5
		stride=4.
	definition="(Tmin<0~S~o~N~C)"
	end if
	
	if (index.eq."ID") then

	spacing=0.5
	definition="(Tmax<0~S~o~N~C)"
	end if
	
	if (index.eq."TN10p") then
	limit_up=15
	spacing=1
	stride=2
	spacing_c=1
	
	end if
	
	if (index.eq."TN50p") then
		limit_up=100
		limit_bot=40
		spacing=2.5
		stride=2.
		spacing_c=2.5
		stride_c=2.
	end if
	
	if (index.eq."TN90p") then
			limit_bot=5.
			limit_up=90.
			spacing=5
			stride=2
			spacing_c=5
	
	end if

	if (index.eq."TX10p") then
		limit_up=15
		spacing=1
		stride=2
		spacing_c=1
	end if
	
	if (index.eq."TX50p") then
		limit_up=100
		limit_bot=40
		spacing=2.5
		stride=2.
		spacing_c=2.5
		stride_c=2.
	end if
	
	if (index.eq."TX90p") then
		limit_bot=5.
		limit_up=30
		spacing=1
		stride=2
		spacing_c=5
	end if
	
	if (index.eq."GSL") then
		limit_bot=240	
		limit_up=365
		spacing=15	
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
	limit_up=10
	spacing=0.5
	stride=2
	spacing_c=1
	stride_c=2
	end if
	
	if (index.eq."CSDI") then
	limit_bot=0
	limit_up=2
	spacing=0.2
	stride=2
	spacing_c=0.2
	stride_c=2
	end if
	
	;PRECIP EXTREMES
	
	if (index.eq."R10mm") then
		limit_up=60
		spacing=2.5
		spacing_c=1
		stride_c=2.
	definition=" "
	end if
	
	if (index.eq."R20mm") then
		limit_up=30
		spacing=2.5
		stride=2
	spacing_c=1
	stride_c=2
	end if
	
	if (index.eq."SDII") then

		limit_up=15
		spacing=0.5
		stride=2
	end if
	
	if (index.eq."Rx1day") then
			limit_up=200
	spacing=10
	stride=2
	end if
	
	if (index.eq."Rx5day") then
			limit_up=400
	spacing=20
	stride=2
	end if
	
	if (index.eq."CDD") then
		spacing_c=1
		stride_c=2.
		limit_up=60
		spacing=2.5
		stride=2
	end if
	
	if (index.eq."CWD") then
		limit_up=15
		spacing=1
		stride=2
	spacing=1
	stride=2
	end if
	
	if (index.eq."R95p") then

		limit_up=600
		spacing=25
		stride=2
	end if
	
	if (index.eq."R99p") then
		
		limit_up=300
		spacing=25
		stride=2
	end if
	
	if (index.eq."PRCPTOT") then
		limit_up=2500
		spacing=100
		stride=2
	end if
	
	
	
	
	;pattern of the filenames
	patt_pres="CCRC_NARCliM_Sydney_MK35_bc_1990-2009_"
	patt_fut="CCRC_NARCliM_Sydney_MK35_bc_2040-2059_"

	patt_imgout = "Sydney_ETCCDI_"
	;---SETUP

	zoom=False

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

	domsize=dimsizes(lat2d)





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

	if (index.eq."TN10p" .or. index.eq."TN50p" .or. index.eq."TN90p" .or. index.eq."TX10p" .or. index.eq."TX50p" .or. index.eq."TX90p") then
		wrf_pres@units="% of days"
	end if

	; 	;#######################
	; 	;--READ ETCCDI FUTURE
	; 	;#######################
	; 
	f2 = addfile(diri_fut+patt_fut+index+".nc","r")

	wrf_fut=f2->Annual
	dsizes_f=dimsizes(wrf_fut)


	wrf_fut@lat2d = lat2d
	wrf_fut@lon2d = lon2d
	wrf_fut!1 = "lat"
	wrf_fut!2 = "lon"
	wrf_fut!0 = "time"
	wrf_fut@units=f2@units
	print(wrf_fut@units)
	printVarSummary(wrf_fut)
	wrf_fut=mask(wrf_fut,wrf_fut.eq."nan",False)
	wrf_fut@_FillValue=wrf_fut@missing_value


	if (index .eq. "GSL") then
		
		wrf_pres2=wrf_pres(0:nyears-2,:,:)
		wrf_fut2=wrf_fut(0:nyears-2,:,:)
		
		delete([/wrf_pres,wrf_fut/])
		
	
		wrf_pres=wrf_pres2
		wrf_fut=wrf_fut2
	end if
	
	wrf_pres(:,0:5,:)=-999.9;
	wrf_pres(:,:,0:5)=-999.9;
	wrf_pres(:,:,domsize(1)-6:domsize(1)-1)=-999.9;
	wrf_pres(:,domsize(0)-6:domsize(0)-1,:)=-999.9;
	
	wrf_fut(:,0:5,:)=-999.9;
	wrf_fut(:,:,0:5)=-999.9;
	wrf_fut(:,:,domsize(1)-6:domsize(1)-1)=-999.9;
	wrf_fut(:,domsize(0)-6:domsize(0)-1,:)=-999.9;
	;##################################
	;#### END OF READING DATA #########
	;##################################

	print("calculate stats")
	wrfann_pres=dim_avg_n(wrf_pres,0)
	wrfann_fut=dim_avg_n(wrf_fut,0)
	wrfann_pres@lat2d = lat2d
	wrfann_pres@lon2d = lon2d
	wrfann_fut@lat2d = lat2d
	wrfann_fut@lon2d = lon2d
	
	; if (index.eq."Rx1day" .or. index.eq."Rx5day")then
	; 	wrfann_pre=dim_max_n(wrf_pres,0)
	; 	wrfann_fut=dim_max_n(wrf_fut,0)
	; end if
	
	
	wrfann_pres=mask(wrfann_pres,lud_wrf,1)
	wrfann_fut=mask(wrfann_fut,lud_wrf,1)
	print(ceil(max((/max(wrfann_pres),max(wrfann_fut)/))/5)*5)

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
	res@gsnMaximize		= True

	res@cnLevelSelectionMode	= "ManualLevels"
	res@cnMinLevelValF			= tofloat(limit_bot)				;same manual levels
	res@cnMaxLevelValF			= tofloat(limit_up)			;as other plots
	res@cnLevelSpacingF			= tofloat(spacing)			;

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


	if (zoom) then
		;zoom to Sydney urban region
		res@mpLeftCornerLatF     = -32.94     ;
		res@mpLeftCornerLonF     = 149.90          ;zoom
		res@mpRightCornerLatF     = -34.54          ;to urban
		res@mpRightCornerLonF     = 151.5          ;
		res@mpProjection				= "LambertConformal"
		res@mpLambertParallel1F		= -27.5
		res@mpLambertParallel2F		= -35.
		res@mpLambertMeridianF		= 146.1

	else
		
			res@mpLeftCornerLatF	= lat2d(0,0)	;
			res@mpLeftCornerLonF	= lon2d(0,0)	;extend beyond strict
			res@mpRightCornerLatF	= lat2d(domsize(0)-1,domsize(1)-1)	;grt syd domain
			res@mpRightCornerLonF	= lon2d(domsize(0)-1,domsize(1)-1)	;
			truelat1 = luf@TRUELAT1
			truelat2 = luf@TRUELAT2
			cen_lat = luf@CEN_LAT
			stand_lon = luf@STAND_LON

			res@mpProjection  = "LambertConformal"
			res@mpLambertParallel1F = truelat1
			res@mpLambertParallel2F = truelat2
			res@mpLambertMeridianF = stand_lon
	end if

	res@mpGridAndLimbOn			= True		;reference lat/
	res@mpGridSpacingF			= 5.			;lon grids

	res@gsnRightString			= ""			;suppress
	res@gsnCenterString			= ""			;labelling
	res@gsnLeftString			= "" 		;for each panel
	res@lbLabelBarOn				= False


	;###### END OF MY RESOURCES #####

	plot(0)=gsn_csm_contour_map(wks,wrfann_pres(:,:),res)
	plot(1)=gsn_csm_contour_map(wks,wrfann_fut(:,:),res)



	;Panel first two plots


	pres1                     = True
	pres1@gsnPanelLabelBar    = True       ; common label bar
	pres1@gsnFrame            = False      ; don't advance frame yet
	pres1@gsnMaximize		= True
	pres1@lbLabelStride   = stride
	pres1@lbTitleString	= wrf_pres@units
	
	pres1@lbTitlePosition			= "Left"
	pres1@lbTitleDirection		= "Across"
	pres1@lbTitleFontHeightF		= 0.014
	pres1@lbLabelFontHeightF = 0.012
	
	pres1@vpWidthF=0.44
	pres1@gsnPanelRight =0.6
	pres1@gsnPanelFigureStrings			= (/"Present","Future"/)
	pres1@gsnPanelFigureStringsFontHeightF= 0.012
	gsn_panel(wks,plot,(/1,2/),pres1)




	;***************************************
	; create third individual plots
	;***************************************
	printVarSummary(ind_changes)
	
	if(index.eq."FD".or.index.eq."ID".or.index.eq."TN10p".or.index.eq."TX10p".or.index.eq."CSDI")then
	gsn_reverse_colormap(wks)
	end if
	
	; res@cnMinLevelValF		= (/floor(min(min(ind_changes))/2)*2/)			;same manual levels
	; 		res@cnMaxLevelValF		= (/ceil(max(max(ind_changes))/2)*2/)			;as other plots
	 res@cnMinLevelValF		= 0
	 res@cnMaxLevelValF		= 90
	res@cnLevelSpacingF	= spacing_c			;
	plot2 = gsn_csm_contour_map(wks,ind_changes(:,:),res)


	pres2                     = True
	pres2@gsnPanelLabelBar    = True       ; common label bar
	pres2@gsnFrame            = False      ; don't advance frame yet
	pres2@gsnMaximize		= True
	pres2@lbLabelStride   = stride_c 
	;pres2@lbTitleString	= wrf_pres@units
	pres2@lbTitlePosition			= "Left"
	pres2@lbTitleDirection		= "Across"
	pres2@lbTitleFontHeightF		= 0.014
	pres2@lbLabelFontHeightF		= 0.012

	pres2@vpWidthF=0.25
	pres2@gsnPanelLeft  = 0.6
	pres2@gsnPanelRight = 0.9
	pres2@gsnPanelFigureStrings			= (/"Changes"/)
	pres2@gsnPanelFigureStringsFontHeightF= 0.012
	gsn_panel(wks,plot2,(/1,1/),pres2)









	txres               = True                     ; text mods desired
	txres@txFontHeightF = 0.02                 ; font smaller. default big
	txres@txJust = "topCenter" ;; justification of text coordinates
	gsn_text_ndc(wks,"Changes in "+index+" "+definition+" MK3.5",0.45,0.8,txres)
	frame(wks)

	; 
	; 
	; 	


	print("_______________END OF PROGRAM_____________________")
end
