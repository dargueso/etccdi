      MODULE COMM		 ! common module used by all parts.
      IMPLICIT NONE
      integer,parameter::DP=kind(1.0)	!  single precision. H.Yang
      integer,parameter:: DoY=365,N_index=33  ! DoY=Day_of_Year

      character*80 :: StnID,ofile_base,ifile,odir
      integer(4)   :: STDSPAN, BASESYEAR, BASEEYEAR, PRCPNN,  &
              SYEAR, EYEAR, TOT, YRS, BYRS, WINSIZE, SS
      real(DP),dimension(:),allocatable,save :: PRCP,Tmax,Tmin
      integer,dimension(12,2)        :: Mon   ! we use new Mon in palce of old Mon and Monleap
      real(DP)   :: LATITUDE, MISSING=-999.9,NoMissingThreshold
      integer(4),allocatable,save :: YMD(:,:), MNASTAT(:,:,:),YNASTAT(:,:)
      integer    :: ID_log,ID_ifile,IDsave,iend,ID_inf
      character*80 :: folder(N_index),OS,Outname,Oname,BPname,log_file,ipt_dir,opt_dir,tmax_dir,tmin_dir,prcp_dir,inf_file,long_names(N_index),tmax_name,tmin_name,prec_name,lon_username,lat_username,time_username
      character*10  :: cmonth(13),ann(1),season(4),units(N_index)
      character(1)  :: sub2(2),sub
      real(DP),allocatable,dimension(:,:,:), save ::data_tmax,data_tmin,data_prcp
      integer :: Nlon,Nlat,Ntime
      real(DP),allocatable,save :: lon(:),lat(:) !,time(:)
      real(DP),allocatable,save :: lon2d(:,:),lat2d(:,:)
      integer(4),dimension(:),allocatable,save :: time ! time is YYYYMMDD, which is actually YYYY0000
!      real(4),dimension(:),allocatable,save :: time ! time is YYYYMMDD, which is actually YYYY0000
      logical   :: Tmax_miss,Tmin_miss,Prcp_miss, &  ! check if ALL data is missing...
                   save_thresholds, &                ! whether the threshold data should be saved.
                   sub_folder,is_rcm,is_thresfile                         ! put indices in different sub_folders or in same folder
			
      data MON    /31,28,31,30,31,30,31,31,30,31,30,31, &
                   31,29,31,30,31,30,31,31,30,31,30,31/
      data WINSIZE/5/
      data sub2/'/','\'/
      data folder/'FD','SU','ID','TR','GSL','TXx','TXn','TNx','TNn','DTR', &
            'R10mm','R20mm','Rnnmm','SDII','Rx1day','Rx5day','CDD','CWD','R95p','R99p', &
!            'PRCPTOT','TX10p','TX50p','TX90p','TN10p','TN50p','TN90p','prcpQC','tempQC','NASTAT', &
!            'PRCPTOT','TX10p','TX50p','TX90p','TN10p','TN50p','TN90p','PRCPmiss','TMAXmiss','TMINmiss', &
            'PRCPTOT','TN10p','TN50p','TN90p','TX10p','TX50p','TX90p','PRCPmiss','TMAXmiss','TMINmiss', &
            'WSDI','CSDI','thresholds'/
      data cMonth/'January','February','March','April','May','June','July','August','September','October','November','December','Annual'/
      data ann/'Annual'/
      data season/'DJF','MAM','JJA','SON'/
      data long_names/'Frost days','Summer days','Ice days','Tropical nights','Growing season Length','Max Tmax','Min Tmax','Max Tmin','Min Tmin','Diurnal temperature range', &
              'Number of days with PRCP >=10 mm','Number of days with PRCP>=20 mm','Number of days with PRCP>=nn mm','Simple daily intensity index', &
                       'Max 1-day precipitation amount','Max 5-day precipitation amount','Consecutive dry days','Consecutive wet days','Very wet days','Extremely wet days', &
              'Annual total wet-day precipitation','Cool nights','average nights','Warm nights','Cool days (Tmax<10th percentile)','average days','Warm days','PRCP mising days', &
                       'TMAX mising days','TMIN mising days', &
              'Warm spell duration indicator','Cold spell duration indicator','thresholds for Temperature'/
      data units/'days','days','days','days','days','degree','degree','degree','degree','degree', &
                'days','days','days','Mm/day','Mm','Mm','days','days','Mm','Mm', &
                'Mm','days','days','days','days','days','days', 'days','days','days',&
                'days','days','degree'/

      END MODULE COMM



module functions
      use COMM, only:DP, MISSING
      real(DP) :: rmiss
      public   :: ismiss,nomiss, leapyear

  contains 

      logical function ismiss(a)
      real(DP) :: a
      rmiss=Missing+1.
      if(a.gt.rmiss) then
        ismiss=.FALSE.
      else
        ismiss=.TRUE.
      endif
      end function ismiss


      logical function nomiss(a)
      real(DP) :: a
      rmiss=Missing+1.
      if(a.lt.rmiss) then
        nomiss=.FALSE.
      else
        nomiss=.TRUE.
      endif
      end function nomiss


! ************ check if it's a leap year ******************** !
      integer function leapyear(iyear)
      integer:: iyear

      if(mod(iyear,400)==0) then
        leapyear=1
      else
        if(mod(iyear,100)==0) then
          leapyear=0
        else
          if(mod(iyear,4)==0) then
            leapyear=1
          else
            leapyear=0
          endif
        endif
      endif

      end function leapyear

 end module functions
