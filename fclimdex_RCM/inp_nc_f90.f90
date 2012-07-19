! We need to read data from netcdf files ...
! 
! different versions of NetCDF has different function name and augments, too confusing...
!
! question:
! 1. the type of lon,lat,time is 6 (double)
!  the type of tmax is 5 (real)
!  However, I read lon,lat,tmax as real, time as integer. Why?
!
subroutine inp_nc_f90
use COMM
use functions
use netcdf
implicit none
integer,dimension(3) :: ncid,ID_time,ID_lon,ID_lat,ID_var,len_time
integer :: i,j,k,it,ky,month,kth,ii,jj,kk,now,nchar,aux,len_lat,len_lon
integer :: ndim(3),indx(1),ndims_file(3)
integer :: len_spatial_dim(3,2),dimids_var(3,3)
real :: miss
character*(120)  :: file(3),title,cname,attribute,data_name(3),name(2)
real,allocatable :: value_coord(:,:,:)
integer,dimension(:),allocatable  :: t1,t2,t3,dimids_coord,dimids_time
real,allocatable,dimension(:,:,:) :: prc,tmx,tmn
character*(10)::lon_name,lat_name,time_name
logical :: exists
data_name(1)=trim(tmax_name)
data_name(2)=trim(tmin_name)
data_name(3)=trim(prec_name)
data lon_name/'lon'/, lat_name/'lat'/,time_name/'time'/


! define the name of the files
! No matter the order in the infilename_temp.txt
! But here, Tmax is the first, Tmin the second and Precip the third.

file(1)=trim(tmax_dir)!//trim(cTmax)
file(2)=trim(tmin_dir)!//trim(cTmin)
file(3)=trim(prcp_dir)!//trim(cPrcp)

!inquire(file=file(2),iostat=status,exist=lexist)
! print*,status,lexist
!print*,'netCDF version ',NF90_INQ_LIBVERS()

!Loop over the three files
do i=1,3

nchar=len_trim(file(i)) ! number of characters
! do they have ".nc" extension?
! if not, add it
if(file(i)(nchar-2:nchar) .ne. '.nc')   file(i)=trim(file(i))//'.nc'

! does the file exists?

inquire(file=file(i),Exist=exists)
if(.Not.exists) then
	print*,'the following input file doesn"t exist, code stops:'
	print*,trim(file(i))
	stop
	endif


!Open the file
!Get the number of dimensions in the file and their ids.
	call err_handle(nf90_open(file(i),NF90_nowrite,NCID(i)), 'open file')
	call err_handle(NF90_INQUIRE(NCID(i),ndims_file(i)),'inqure dim/var numbers')

enddo
!print*,'The files have',ndims_file,'dimensions'

!################### READ THE SPATIAL COORDINATES #########################################
!get the coordinates (dimensions) and check that the are the same for all three data....
!LONGITUDE


! Read longitude variable in all three files
do i=1,3 !Loop over the files
	call err_handle(NF90_INQ_VARID (NCID(i), lon_name, ID_lon(i)),'inquire longitude var ID')
	call err_handle(NF90_INQUIRE_VARIABLE(NCID(i),ID_lon(i),ndims=ndim(i)),'inquire longitude var dimensions')
end do


!Longitude in each of the files must have same number of dimensions (Generally 1 or 2)

if (ndim(1)/=ndim(2) .or. ndim(1)/=ndim(3)) stop 'Longitude variable has different dimensions in each file'


allocate(dimids_coord(ndim(1)))




do i=1,3
	call err_handle(NF90_INQUIRE_VARIABLE(NCID(i),ID_lon(i),ndims=ndim(i),dimIDs=dimids_coord),'inquire var dimensionsm IDs') 
	!IDs are always provided from the last to the first in the variable and the number refers their order in the FILE.
	do j=1,ndim(1)
		call err_handle(NF90_INQUIRE_DIMENSION (NCID(i),dimids_coord(j),name(j),len=len_spatial_dim(i,j)),'inq dimensions')
		enddo
		!print*,'grid sizes are ID:',dimids_coord,'size',len_spatial_dim(i,:)
		enddo
		deallocate (dimids_coord)


!Check the size of longitude
! y and x are always the rightmost dimensions (time might be included in the dimensions of the coordinate variables
!Time is omitted in this checking
if (len_spatial_dim(1,1)/=len_spatial_dim(2,1) .or. len_spatial_dim(1,1)/=len_spatial_dim(3,1)) stop 	'Error: different dimension sizes for lon ...'
if (ndim(1) .lt. 1) then
	if (len_spatial_dim(1,2)/=len_spatial_dim(2,2) .or. len_spatial_dim(1,2)/=len_spatial_dim(3,2)) stop 	'Error: different dimension sizes for lon ...'
end if


!Check that the longitude grid cover the same range.

if  (ndim(1)==1) then!Unidimensional coordinantes
	allocate(value_coord(len_spatial_dim(1,1),1,3))

do i=1,3
call err_handle(NF90_get_var(ncid(i),ID_lon(i),value_coord(:,1,i)),'get Var real')
	enddo


if(value_coord(1,1,1)/=value_coord(1,1,2).or.value_coord(1,1,1)/=value_coord(1,1,3)) stop 'Error: different lon range ...'


if(value_coord(len_spatial_dim(1,1),1,1)/=value_coord(len_spatial_dim(1,1),1,2).or.value_coord(len_spatial_dim(1,1),1,1)/=value_coord(len_spatial_dim(1,1),1,3)) stop 'Error: different lon range ...'
allocate(lon(len_spatial_dim(1,1)))
lon=value_coord(:,1,1)
else!Bidimensional coordinates

allocate(value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),3))

do i=1,3
call err_handle(NF90_get_var(ncid(i),ID_lon(i),value_coord(:,:,i)),'get Var real')
	enddo


if(value_coord(1,1,1)/=value_coord(1,1,2).or.value_coord(1,1,1)/=value_coord(1,1,3)) stop 'Error: different lon range ...'
if(value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),1)/=value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),2).or.value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),1)/=value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),3)) stop 'Error: different lon range ...'

allocate(lon(len_spatial_dim(1,1)))
lon=value_coord(:,len_spatial_dim(1,2)/2+1,1) !Gettin the central line approximately. Useful for RCMs non-regular grids

end if




! If the grid is not regular:

if (is_RCM) then

print*,'warning: you have set "is_rcm" to', is_rcm
print*,'that means your grid is not regular'
print*,'warning: this might cause wrong GSL values'


allocate(lon2d(len_spatial_dim(1,1),len_spatial_dim(1,2)))
lon2d=value_coord(:,:,1)
end if
deallocate(value_coord)

!LATITUDE
do i=1,3 !Loop over the files
	call err_handle(NF90_INQ_VARID (NCID(i), lat_name, ID_lat(i)),'inq var ID')

	call err_handle(NF90_INQUIRE_VARIABLE(NCID(i),ID_lat(i),ndims=ndim(i)),'inquire var dimensions')
enddo

if (ndim(1)/=ndim(2) .or. ndim(1)/=ndim(3)) stop 'Latitude variable has different dimensions in each file'

allocate(dimids_coord(ndim(1)))

do i=1,3
	call err_handle(NF90_INQUIRE_VARIABLE(NCID(i),ID_lat(i),ndims=ndim(i),dimIDs=dimids_coord),'inquire var dimensionsm IDs')
	do j=1,ndim(1)
	call err_handle(NF90_INQUIRE_DIMENSION (NCID(i),dimids_coord(j),len=len_spatial_dim(i,j)),'inq dimensions')
		enddo
		enddo
		deallocate (dimids_coord)

!Checking size of latitude
if (len_spatial_dim(1,1)/=len_spatial_dim(2,1) .or. len_spatial_dim(1,1)/=len_spatial_dim(3,1)) stop 	'Error: different dimension sizes for lat ...'
if (ndim(1) .lt. 1) then
	if (len_spatial_dim(1,2)/=len_spatial_dim(2,2) .or. len_spatial_dim(1,2)/=len_spatial_dim(3,2)) stop 	'Error: different dimension sizes for lat ...'
end if

!Check that the latitude grid cover the same range.
if (ndim(1)==1) then !Unidimensional coordinates
	allocate(value_coord(len_spatial_dim(1,1),1,3))
	do i=1,3
call err_handle(NF90_get_var(ncid(i),ID_lat(i),value_coord(:,1,i)),'get Var real')
		enddo
		if(value_coord(1,1,1)/=value_coord(1,1,2).or.value_coord(1,1,1)/=value_coord(1,1,3)) stop 'Error: different lat range ...'
		if(value_coord(len_spatial_dim(1,1),1,1)/=value_coord(len_spatial_dim(1,1),1,2).or.value_coord(len_spatial_dim(1,1),1,1)/=value_coord(len_spatial_dim(1,1),1,3)) stop 'Error: different lat range ...'
		allocate(lat(len_spatial_dim(1,1)))
		lat=value_coord(:,1,1)
		else !Bidimensional coordinates

allocate(value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),3))
do i=1,3
call err_handle(NF90_get_var(ncid(i),ID_lat(i),value_coord(:,:,i)),'get Var real')
	enddo
	if(value_coord(1,1,1)/=value_coord(1,1,2).or.value_coord(1,1,1)/=value_coord(1,1,3)) stop 'Error: different lat range ...'
	if(value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),1)/=value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),2).or.value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),1)/=value_coord(len_spatial_dim(1,1),len_spatial_dim(1,2),3)) stop 'Error: different lat range ...'
	allocate(lat(len_spatial_dim(1,2)))
	!Gettin the central line approximately. Useful for RCMs non-regular grids
	lat=value_coord(len_spatial_dim(1,1)/2+1,:,1)
end if


if (is_RCM) then

allocate(lat2d(len_spatial_dim(1,1),len_spatial_dim(1,2)))
lat2d=value_coord(:,:,1)
end if
deallocate(value_coord)

!definitive size of the arrays after checks
len_lat=size(lat)
len_lon=size(lon)

!################ END OF READING SPATIAL COORDINATES ############################3


!################ READING TIME DIMENSION ################################

! get dimension variable - time
do i=1,3
		call err_handle(NF90_INQ_VARID (NCID(i), time_name, ID_time(i)),'inq var ID')
	call err_handle(NF90_INQUIRE_VARIABLE(NCID(i),ID_time(i),ndims=ndim(i)),'inq number of dimensions')
	allocate(dimids_time(ndim(i)))
	call err_handle(NF90_INQUIRE_VARIABLE(NCID(i),ID_time(i),ndims=ndim(i),dimIDs=dimids_time),'inq var dimension IDs')
	call err_handle(NF90_INQUIRE_DIMENSION (NCID(i),dimids_time(1),len=len_time(i)),'inq dimensions for length') 
	deallocate(dimids_time)     
	enddo
	if (len_time(1) /= len_time(2)) stop 'Maximum and mininum temperatures do not have the same time coverage'

attribute=''
allocate(t1(len_time(1)),t2(len_time(2)),t3(len_time(3))) !time dimension lenght might be different


!Convert time units into the format YYYYMMDD

	call err_handle(NF90_get_var(ncid(1),ID_time(1),t1),'get Var int -- t1')
If(NF90_GET_ATT(ncid(1),ID_time(1),'units',attribute) == NF90_noerr) call time_convert(len_time(1),t1,attribute)


attribute=''
	call err_handle(NF90_get_var(ncid(2),ID_time(2),t2),'get Var int -- t2')
If(NF90_GET_ATT(ncid(2),ID_time(2),'units',attribute) == NF90_noerr) call time_convert(len_time(2),t2,attribute)

attribute=''
	call err_handle(NF90_get_var(ncid(3),ID_time(3),t3),'get Var int -- t3')
If(NF90_GET_ATT(ncid(3),ID_time(3),'units',attribute) == NF90_noerr) call time_convert(len_time(3),t3,attribute)

!The internal order of the files are TMAX, TMIN and PRECIP.
!It does not depend on the order we provided them in the infilename_temp.txt file
! get data Var...

do i=1,3 !Loop over the files
	call err_handle(NF90_INQ_VARID (NCID(i),data_name(i),ID_var(i)),'inq var IDs: Variable not loaded. Check namelist')

	call err_handle(NF90_INQUIRE_VARIABLE(NCID(i),ID_var(i),ndims=ndim(i),dimIDs=dimids_var(i,:)),'inq var ID for dimensions')
!print*,'variable dimension IDs', dimids_var(i,:)
end do

if (ndim(1)/= 3 .or. ndim(2)/=3 .or. ndim(3)/=3) stop 'The variables do not have three dimensions'

! GET MAXIMUM TEMPERATURE
allocate(tmx(len_lon,len_lat,len_time(1))) ! First dimension is time, before-last is latitude and last is longitude according to CF standars.


call err_handle(NF90_get_var(ncid(1),ID_var(1),tmx),'get Var real -- Tmax')
if(NF90_GET_ATT(NCID(1), ID_var(1),'missing_value',miss) == NF90_noerr) then
	where(tmx == miss) tmx=MISSING
		endif
		if(NF90_GET_ATT(NCID(1), ID_var(1),'_FillValue',miss) == NF90_noerr) then
			where(tmx == miss) tmx=MISSING
				endif

!GET MINIMUM TEMPERATURE
allocate(tmn(len_lon,len_lat,len_time(2)))

call err_handle(NF90_get_var(ncid(2),ID_var(2),tmn),'get Var real -- Tmin')
if(NF90_GET_ATT(NCID(2), ID_var(2),'missing_value',miss) == NF90_noerr) then
	where(tmn == miss) tmn=MISSING
		endif
		if(NF90_GET_ATT(NCID(2), ID_var(2),'_FillValue',miss) == NF90_noerr) then
			where(tmn == miss) tmn=MISSING
				endif


!GET PRECIPITATION
allocate(prc(len_lon,len_lat,len_time(3)))

call err_handle(NF90_get_var(ncid(3),ID_var(3),prc),'get Var real -- Prec')
if(NF90_GET_ATT(NCID(2), ID_var(3),'missing_value',miss) == NF90_noerr) then
	where(prc == miss) prc=MISSING
		endif
		if(NF90_GET_ATT(NCID(3), ID_var(3),'_FillValue',miss) == NF90_noerr) then
			where(prc == miss) prc=MISSING
				endif


!CLOSE THE FILES!
do i=1,3
	call err_handle(NF90_CLOSE(NCID(i)),'close file')
	enddo

Syear=min(t1(1),t2(1),t3(1))/10000
Eyear=max(t1(len_time(1)),t2(len_time(2)),t3(len_time(3)))/10000
YRS=Eyear-Syear+1





allocate(time(YRS))
tot=0
do i=Syear,Eyear
	tot=tot+365
	if(leapyear(i)==1) tot=tot+1
	time(i+1-Syear)=i*10000.+101.
	enddo
	print*,'Total data #:',tot

allocate(YMD(tot,3),data_tmax(len_lon,len_lat,tot),data_tmin(len_lon,len_lat,tot),data_prcp(len_lon,len_lat,tot))

it=0
do i=SYEAR,EYEAR
	Ky=leapyear(i)+1
	do month=1,12
		Kth=Mon(month,Ky)
		do k=1,kth
			it=it+1
			!      if(abs(it-790)<3) print*,it,ky,kth,k
			YMD(it,1)=i
			YMD(it,2)=month
			YMD(it,3)=k
			!        if(i==1992 .and. month==9) print*,k,it,tmx(1,1,it),tmn(1,1,it),prc(1,1,it)
			enddo
			enddo
			enddo
			if(it/=tot) stop 'Error: different TOT and IT in inp_nc !'
			!     print*,tot
			!stop

jj=1 ! This was part of the original QC subroutine -- get full data and insert MISSING values..
ii=1
kk=1
do i=1, tot
	now=YMD(i,1)*10000+YMD(i,2)*100+YMD(i,3)


if(ii>len_time(1)) then
	!         
	data_tmax(:,:,i)=MISSING
	else
		if(now.eq.t1(ii)) then
			data_tmax(:,:,i)=tmx(:,:,ii)
			ii=ii+1
			else
				data_tmax(:,:,i)=MISSING
				endif
				endif

if(jj>len_time(2)) then
	data_tmin(:,:,i)=MISSING
	else
		if(now.eq.t2(jj)) then
			data_tmin(:,:,i)=tmn(:,:,jj)
			jj=jj+1
			else
				data_tmin(:,:,i)=MISSING
				endif
				endif

if(kk>len_time(3)) then
	data_prcp(:,:,i)=MISSING
	else
		if(now.eq.t3(kk)) then
			data_prcp(:,:,i)=prc(:,:,kk)
			kk=kk+1
			else
				data_prcp(:,:,i)=MISSING
				endif
				endif

!     if(YMD(i,1)==1992 .and. YMD(i,2) ==9 .and. YMD(i,3)==1) then
!     endif

enddo

Nlon=len_lon;  Nlat=len_lat;  Ntime=tot !ndim(3)
deallocate(t1,t2,t3,prc,tmx,tmn)

print*,'finish reading Tmax,Tmin and PRCP data ...'
print*, size(data_tmax,1)
print*, size(data_tmax,2)
print*, size(data_tmax,3)

! stop

return
end subroutine inp_nc_f90