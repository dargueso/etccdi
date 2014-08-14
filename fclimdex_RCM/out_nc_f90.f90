!...   This program is to create an NetCDF file.
!...   And the dimention is three dimentions (lon, lat, time)

    subroutine out_nc_f90(nindx,nlevel,data)
    use COMM
    use netcdf

!... define the variables used in this program !
    integer nindx, ncid,nlevel,Ititle
    integer LATID, LONID, TIMEID                ! dimention's ID
    integer varID_lon,varID_lat,varID_time      ! variable's ID
    integer dataDim(3),coordim(2)
    integer,allocatable :: ID_data(:)
    real data(Nlon, Nlat, YRS,nlevel)      ! variable array
    integer start(3), count(3)
    data start /1,1,1/
    character*150 O_file
    character*100 ctmp
    character*10, allocatable :: ctit(:)
    character*100 ctitle(3)
    data ctitle/', 13 levels -- 12 monthly and 1 Annual values.',', 1 level -- Annual value.',', 4 levels -- 4 seasons.'/

    if(nlevel/=13 .and. nlevel/=1 .and. nlevel/=4) stop 'Error: nlevel in out_nc is wrong !'
    allocate(ctit(nlevel),ID_data(nlevel))

    if(nlevel==13) then
        Ititle=1
        ctit=cmonth
    else if(nlevel==1 ) then
        Ititle=2
        ctit=ann
    else
        Ititle=3
        ctit=season
    endif
    count(1)=Nlon; count(2)=Nlat; count(3)=YRS
!        write(O_file,'(a,a,"/",a,"")') trim(opt_dir),folder(nindx),folder(nindx)
    write(ctmp,'("index ",a,a)') trim(folder(nindx)),trim(ctitle(Ititle))
    if(sub_folder) then
      O_file=trim(opt_dir)//trim(folder(nindx))//sub//trim(Oname)//trim(folder(nindx))//'.nc'
    else
      O_file=trim(opt_dir)//trim(Oname)//trim(folder(nindx))//'.nc'
    endif
!        print*,O_file,ctmp

!... begin to creat a new NetCDF file 
    call err_handle(NF90_create(O_file, NF90_clobber, ncid),'create nc file')
!... begin to define your dimentions
    call err_handle(NF90_def_dim(ncid, 'lon', Nlon, LONID),'define lon dimension')
    call err_handle(NF90_def_dim(ncid, 'lat', Nlat, LATID),'define lat dimension')
    call err_handle(NF90_def_dim(ncid, 'time', YRS, TIMEID),'define time dimension')

if (is_rcm) then

	coordim(1)=lonid
	coordim(2)=latid
    call err_handle(NF90_def_var(ncid, 'lon', NF90_float, coordim,varID_lon),'define var lon')
    call err_handle(NF90_def_var(ncid, 'lat', NF90_float, coordim,varID_lat),'define var lat')
    call err_handle(NF90_def_var(ncid, 'time', NF90_int, timeID,varID_time),'define var time')


else
    call err_handle(NF90_def_var(ncid, 'lon', NF90_float, lonID,varID_lon),'define var lon')
    call err_handle(NF90_def_var(ncid, 'lat', NF90_float, latID,varID_lat),'define var lat')
    call err_handle(NF90_def_var(ncid, 'time', NF90_int, timeID,varID_time),'define var time')
end if


!... begin to define varibles
    datadim(1)=lonid
    datadim(2)=latid
    datadim(3)=timeid
    do i=1,nlevel
       call err_handle(NF90_def_var(ncid, trim(ctit(i)), NF90_float, datadim, ID_data(i)),'define data')
       call err_handle(NF90_put_att(ncid, ID_data(i), 'missing_value',MISSING),'put att for data')
       call err_handle(NF90_put_att(ncid, ID_data(i), '_FillValue', MISSING),'put att for data')
    enddo

!... Put attribute of variables
    call err_handle(NF90_put_att(ncid, NF90_global, 'Title', trim(ctmp)),'define global title')
    call err_handle(NF90_put_att(ncid, NF90_global, 'author',  'Hongang Yang - hongang.yang@unsw.edu.au ') ,'define global author')
    call err_handle(NF90_put_att(ncid, NF90_global, 'history', 'Created from Fclimdex version 3.1.1 '),'define global history')
    call err_handle(NF90_put_att(ncid, NF90_global, 'units',  trim(units(nindx))),'define global unit')
    call err_handle(NF90_put_att(ncid, NF90_global, 'long_name', trim(long_names(nindx))),'define global long_name')
!    call err_handle(NF90_put_att(ncid, NF90_global, 'missing_value', ,MISSING),'put att for data')
!    call err_handle(NF90_put_att(ncid, NF90_global, '_FillValue', MISSING),'put att for data')

    call err_handle(NF90_put_att(ncid, varID_lon, 'long_name', 'Longitude'),'put att for lon')
    call err_handle(NF90_put_att(ncid, varID_lon, 'units', 'degrees_east'),'put att for lon')
    call err_handle(NF90_put_att(ncid, varID_lon, 'axis', 'X'),'put att for lon')

    call err_handle(NF90_put_att(ncid, varID_lat, 'long_name', 'Latitude'),'put att for lat')
    call err_handle(NF90_put_att(ncid, varID_lat, 'units',  'degrees_north'),'put att for lat')
    call err_handle(NF90_put_att(ncid, varID_lat, 'axis', 'Y'),'put att for lat')

!	call err_handle(NF90_put_att(ncid, varID_time, 'units','hours since 1-1-1 0:0:0')
    call err_handle(NF90_put_att(ncid, varID_time, 'units', 'day as %Y%m%d.%f'),'put att for time')
    call err_handle(NF90_put_att(ncid, varID_time, 'calendar', 'proleptic_gregorian'),'put att for time')


!... end of define mode
    call err_handle(NF90_enddef(ncid),'end define')

!... To get data !

if (is_rcm) then
	!print*, lon2d(1:2,:)
	!print*, lat2d(:,1:2)
    call err_handle(NF90_put_var(ncid, varID_lon, lon2d),'save lon var')
    call err_handle(NF90_put_var(ncid, varID_lat, lat2d),'save lat var')
    call err_handle(NF90_put_var(ncid, varID_time, time),'save time var')
else
    call err_handle(NF90_put_var(ncid, varID_lon, lon),'save lon var')
    call err_handle(NF90_put_var(ncid, varID_lat, lat),'save lat var')
    call err_handle(NF90_put_var(ncid, varID_time, time),'save time var')
end if
!... begin to put data in to this file!
     do i=1,nlevel
        call err_handle(NF90_put_var(ncid, ID_data(i), data(:,:,:,i)),'save data')
     enddo

!... close this file..
    call err_handle(NF90_close(ncid),'close nc file')
    deallocate(ctit,ID_data)

     return
     end subroutine out_nc_f90



   subroutine out_nc2_f90(nindx,data)
    use COMM
    use netcdf

!... define the variables used in this program !
    integer nindx, ncid,ID_data(6)
    integer LATID, LONID, TIMEID                ! dimention's ID
    integer varID_lon,varID_lat,varID_time      ! variable's ID
    integer dataDim(3),coordim(2)
    real data(Nlon, Nlat, DoY,6)      ! variable array
    integer start(3), count(3)
    character*150 :: O_file
    character*100 :: ctmp, ctit(6)

    data start /1,1,1/
    data ctit/'thresan10','thresan50','thresan90','thresax10','thresax50','thresax90'/
    
    print*,'saving threshold data ...'
    count(1)=Nlon; count(2)=Nlat; count(3)=DoY
    write(ctmp,'("index ",a,a)') trim(folder(nindx)),' 6 levels: 3 for Tmin, 3 for Tmax.'

    if(sub_folder) then
      O_file=trim(opt_dir)//trim(folder(nindx))//sub//trim(Oname)//trim(folder(nindx))//'.nc'
    else
      O_file=trim(opt_dir)//trim(Oname)//trim(folder(nindx))//'.nc'
    endif

!... begin to creat a new NetCDF file 
    call err_handle(NF90_create(O_file, NF90_clobber, ncid),'create nc file')

!... begin to define your dimentions
    call err_handle(NF90_def_dim(ncid, 'lon', Nlon, LONID),'define lon dimension')
    call err_handle(NF90_def_dim(ncid, 'lat', Nlat, LATID),'define lat dimension')
    call err_handle(NF90_def_dim(ncid, 'doy', DoY, TIMEID),'define time dimension')
if (is_rcm) then
	coordim(1)=lonid
	coordim(2)=latid
    call err_handle(NF90_def_var(ncid, 'lon', NF90_float, coordim,varID_lon),'define var lon')
    call err_handle(NF90_def_var(ncid, 'lat', NF90_float, coordim,varID_lat),'define var lat')
   ! call err_handle(NF90_def_var(ncid, 'lon', NF90_float, lonID,varID_lon),'define var lon')
   ! call err_handle(NF90_def_var(ncid, 'lat', NF90_float, latID,varID_lat),'define var lat')
    call err_handle(NF90_def_var(ncid, 'doy', NF90_int, timeID,varID_time),'define var time')
else

    call err_handle(NF90_def_var(ncid, 'lon', NF90_float, lonID,varID_lon),'define var lon')
    call err_handle(NF90_def_var(ncid, 'lat', NF90_float, latID,varID_lat),'define var lat')
    call err_handle(NF90_def_var(ncid, 'doy', NF90_int, timeID,varID_time),'define var time')
end if
!... begin to define varibles
    datadim(1)=lonid
    datadim(2)=latid
    datadim(3)=timeid
    do i=1,6
       call err_handle(NF90_def_var(ncid, trim(ctit(i)), NF90_float, datadim, ID_data(i)),'define data')
       call err_handle(NF90_put_att(ncid, ID_data(i), 'missing_value', MISSING),'put att for data')
       call err_handle(NF90_put_att(ncid, ID_data(i), '_FillValue', MISSING),'put att for data')
    enddo

!... Put attribute of variables
    call err_handle(NF90_put_att(ncid, NF90_global, 'Title', trim(ctmp)),'define global title')
    call err_handle(NF90_put_att(ncid, NF90_global, 'author', 'Hongang Yang - hongang.yang@unsw.edu.au ') ,'define global author')
    call err_handle(NF90_put_att(ncid, NF90_global, 'history', 'Created from Fclimdex version 3.1.1 '),'define global history')
    call err_handle(NF90_put_att(ncid, NF90_global, 'units',trim(units(nindx))),'define global unit')
    call err_handle(NF90_put_att(ncid, NF90_global, 'long_name',trim(long_names(nindx))), 'define global long_name')

    call err_handle(NF90_put_att(ncid, varID_lon, 'long_name', 'Longitude'),'put att for lon')
    call err_handle(NF90_put_att(ncid, varID_lon, 'units', 'degrees_east'),'put att for lon')
    call err_handle(NF90_put_att(ncid, varID_lon, 'axis',  'X'),'put att for lon')

    call err_handle(NF90_put_att(ncid, varID_lat, 'long_name', 'Latitude'),'put att for lat')
    call err_handle(NF90_put_att(ncid, varID_lat, 'units', 'degrees_north'),'put att for lat')
    call err_handle(NF90_put_att(ncid, varID_lat, 'axis', 'Y'),'put att for lat')

    call err_handle(NF90_put_att(ncid, varID_time, 'units', 'day since 1-1-0 00:00:00'),'put att for time')

!... end of define mode
    call err_handle(NF90_enddef(ncid),'end define')

!... To get data !
if (is_rcm) then
    call err_handle(NF90_put_var(ncid, varID_lon, lon2d),'save lon var')
    call err_handle(NF90_put_var(ncid, varID_lat, lat2d),'save lat var')
    call err_handle(NF90_put_var(ncid, varID_time, (/(i,i=1,DoY)/)),'save time var')
else

    call err_handle(NF90_put_var(ncid, varID_lon, lon),'save lon var')
    call err_handle(NF90_put_var(ncid, varID_lat, lat),'save lat var')
    call err_handle(NF90_put_var(ncid, varID_time, (/(i,i=1,DoY)/)),'save time var')
end if
!... begin to put data in to this file!
     do i=1,6
        call err_handle(NF90_put_var(ncid, ID_data(i), data(:,:,:,i)),'save data')
     enddo
!... close this file..
    call err_handle(NF90_close(ncid),'close nc file')

     return
     end subroutine out_nc2_f90
