!        last modified 2008-05-06
!  add TMAXmean and TMINmean output in q!  function
!  in TN10p subroutine, add an 1e-5 term on all
!  thresholds, to eliminate computational error. ( 3.5 may store like
!  3.5000001 or 3.49999999 in thresholds )
!  changed TN10p subroutine, set missing value for monthly output the
!  same level as R, eg. >10 days missing in a month, then set this month
!  missing.
!   last modified 2008-06-15
!  changed percentile funcion to calculate multi-level percentiles in a single
!  routine, also changed threshold.
!
!  Based on the version provided by Lisa Alexander @ 2011.4.11
!  To simplify, modify and improve the code by H.Yang from 2011.4.12
!
! Version 3: read and output Netcdf files (using Fortran77 interface)
! also output thresholds for DoY
! change the input to be smarter, and output nc files with missingvalues ...
!   3.1:
!  input file is smarter.
!  check if output are in one folder or different sub folders.


! **********     Main program starts from here	 ********************* !
program Fclimdex
use COMM
integer :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,IDcpu,ID,status,i
integer(4)            :: stnnum
!      real(DP),dimension(2) :: time
double precision     :: tfrom,tto
real,allocatable,dimension(:,:,:,:) :: FDout,Rnmout,CDDout,R95pout,wcsdi,thresout
real,allocatable,dimension(:,:,:,:,:) :: QCout,TXXout,RXout,TX10pout
real,allocatable,dimension(:,:,:) :: GSLout,thresoutpr
!      real,allocatable :: indx_mon(:,:,:,:),indx_year(:,:,:)

namelist /input/ ipt_dir,opt_dir,para_file,tmax_dir,tmin_dir,prcp_dir,inf_file,log_file,  &
OS,save_thresholds,sub_folder,NoMissingThreshold,is_rcm,is_thresfile,&
prec_name,tmax_name,tmin_name, STDSPAN, BASESYEAR, BASEEYEAR, PRCPNN, Outname, &
lon_username,lat_username,time_username

!$OMP parallel private(IDcpu);  IDcpu=omp_get_thread_num();   print*,'OpenMP ID #',IDcpu,'of ',omp_get_num_threads()
!$OMP end parallel

SS=int(WINSIZE/2)
NoMissingThreshold=0.85  ! default value.        ! suggested by Li
!      Ofile_base='./data/climdex'  ! for Linux
!      Ofile_base='.\data\climdex' ! for Windows
save_thresholds=.false.
sub_folder=.false.



call print_head

call getID(ID)
open (ID, file="input.nml",status='OLD')
read(ID,nml=input)
close(ID)
if(OS/='windows'.and.OS/='linux') stop 'Error in input.nml: OS wrong !'

call getID(ID_log)
open (ID_log, file=trim(opt_dir)//trim(log_file),IOSTAT=status)
if(status.ne.0) then
close(ID_log)
stop 'ERROR in input.nml: the opt_dir does NOT exist !'
endif

call get_time(0,ID_log,tfrom)
if(sub_folder) then
call check_folders  ! check if sub_dir exists, if not, create it.
else
print*,'all data will be put in one folder...'
endif

stnnum=1

call inp_nc_f90
print('(a,i,2x,i,2x,i)'),'dimensions (Nlon, Nlat, YRS): ',Nlon,Nlat,YRS

write(Oname,'(a,"_",i4.4,"-",i4.4,"_")') trim(Outname),Syear,Eyear      ! name for output
	print*,'THE OUTPUT NAME IS: ',Oname
!else
!stop 'Error: the input data file names have different header ...'
!endif

BYRS=BASEEYEAR-BASESYEAR+1

allocate(Tmax(tot),Tmin(tot),PRCP(tot))
allocate(MNASTAT(YRS,12,3),YNASTAT(YRS,3),thresout(Nlon,Nlat,DoY,6),thresoutpr(Nlon,Nlat,2))
allocate(QCout(Nlon,Nlat,YRS,13,3),FDout(Nlon,Nlat,YRS,4),GSLout(Nlon,Nlat,YRS),TXXout(Nlon,Nlat,YRS,13,5), &
Rnmout(Nlon,Nlat,YRS,4),RXout(Nlon,Nlat,YRS,13,2),CDDout(Nlon,Nlat,YRS,2),R95pout(Nlon,Nlat,YRS,3), &
TX10pout(Nlon,Nlat,YRS,13,6),wcsdi(Nlon,Nlat,YRS,2))

thresout=MISSING; thresoutpr=MISSING; QCout=MISSING; FDout=MISSING; GSLout=MISSING; TXXout=MISSING; Rnmout=MISSING; RXout=MISSING
CDDout=MISSING; R95pout=MISSING; TX10pout=MISSING; wcsdi=MISSING

do lat_j=1,Nlat
LATITUDE=lat(lat_j)
do lon_i=1,Nlon
TMAX(:)=data_tmax(lon_i,lat_j,:)
TMIN(:)=data_tmin(lon_i,lat_j,:)
PRCP(:)=data_prcp(lon_i,lat_j,:)

Tmax_miss=.false.; Tmin_miss=.false.; Prcp_miss=.false.

if(count(TMAX(1:TOT) .eq. MISSING) .ge. TOT) Tmax_miss=.true.
if(count(TMIN(1:TOT) .eq. MISSING) .ge. TOT) Tmin_miss=.true.
if(count(PRCP(1:TOT) .eq. MISSING) .ge. TOT) Prcp_miss=.true.
if(Tmax_miss .and. Tmin_miss .and. Prcp_miss) goto 78
!             print*,'All PRCP are missing, so we will NOT calculate/output Prcp-related indices !'
!             goto,78
!          endif

call qc(QCout(lon_i,lat_j,:,:,:))   ! PRCPmiss,TMAXmiss,TMINmiss (13 months) in NASTAT 
call FD(FDout(lon_i,lat_j,:,:))     ! FD, SU, ID, TR (Annual)
call GSL(GSLout(lon_i,lat_j,:))    ! GSL (Annual)
call TXX(TXXout(lon_i,lat_j,:,:,:))    ! TXx, TXn, TNx, TNn, DTR (13 months)
call Rnnmm(Rnmout(lon_i,lat_j,:,:))  ! R10mm, R20mm, Rnnmm, SDII (Annual)
call RX5day(RXout(lon_i,lat_j,:,:,:)) ! Rx1day, Rx5day (13 months)
call CDD(CDDout(lon_i,lat_j,:,:))    ! CDD, CWD (Annual)
call R95p(R95pout(lon_i,lat_j,:,:),thresoutpr(lon_i,lat_j,:),lon_i,lat_j) ! R95p, R99p, PRCPtot (Annual)
call TX10p(TX10Pout(lon_i,lat_j,:,:,:),wcsdi(lon_i,lat_j,:,:),thresout(lon_i,lat_j,:,:),lon_i,lat_j)   ! tn10p,tn50p,tn90p,tx10p,tx50p,tx90p (13 months); wsdi,csdi (Annual)
78      continue
enddo
enddo
deallocate(Tmax,Tmin,PRCP)
print*,'finish calculating indices ...'

do i=1,3
call out_nc_f90(27+i,13,QCout(:,:,:,:,i))
enddo

do i=1,4
call out_nc_f90(i,1,FDout(:,:,:,i))
enddo
!   print*,'GSL'
call out_nc_f90(5,1,GSLout(:,:,:))
!   print*,'TXX'

do i=1,5
call out_nc_f90(5+i,13,TXXout(:,:,:,:,i))
enddo

do i=1,4
call out_nc_f90(10+i,1,Rnmout(:,:,:,i))
enddo

do i=1,2
call out_nc_f90(14+i,13,RXout(:,:,:,:,i))
enddo

do i=1,2
call out_nc_f90(16+i,1,CDDout(:,:,:,i))
enddo

do i=1,3
call out_nc_f90(18+i,1,R95pout(:,:,:,i))
enddo

do i=1,6
call out_nc_f90(21+i,13,TX10pout(:,:,:,:,i))
enddo

do i=1,2
call out_nc_f90(30+i,1,wcsdi(:,:,:,i))
enddo

!  print*,'thresholds:',thresout(1,1,1:5,1:2),thresout(1,1,360:365,1:2)
!  print*,maxval(thresout(:,:,1:5,:)),minval(thresout(:,:,1:5,:))

if(save_thresholds) call out_nc2_f90(33,thresout,thresoutpr)
call get_time(2,ID_log,tto)
print*,'finish saving netcdf data ...'

deallocate(MNASTAT,YNASTAT,QCout,FDout,GSLout,TXXout, lon,lat,time, &
Rnmout,RXout,CDDout,R95pout,TX10pout,wcsdi,YMD,data_tmax,data_tmin,data_prcp,thresout,thresoutpr)
stnnum=stnnum+1

write(ID_log,*) "Total ",stnnum," dataset(s) were calculated !"
call get_time(1,ID_log,tto)
close(ID_log)

print*,'Total time usage is about ',real(tto-tfrom),' seconds..'

stop
end program Fclimdex
