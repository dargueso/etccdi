   subroutine print_head
   print*,'------ Note for netcdf data -------'
   print*,'ATENTION: This is the NetCDFversion of fclimdex'
   print*,'There is another version for station data'
   print*,'In this version, non-regular grids are allowed (RCMs)'
   print*,'and there is only one namelist (input.nml).'
   print*,'There are no restrictions in the name of the files or the variables,'
   print*,'you must specify them in the name list.'
   print*,'The only restriction is that the NetCDF files must follow the CF conventions'
   print*,'Most commonly used dataset follow these conventions.'
   print*,' -------- Now continue ... ---------'
   print*,''
   end subroutine print_head




! check if sub_dir exists, if not, create it.
      subroutine check_folders
      use COMM
      character(150) :: command,tempf
      integer        :: status,i,system

      if(OS=='windows') then
      	sub=sub2(2)
        do i=1,N_index
           tempf=trim(opt_dir)//trim(folder(i))
           if(opt_dir(2:2)==':') then
      	      command='if not exist '//trim(tempf)//sub//'nul '//opt_dir(1:2)//' && cd '//trim(opt_dir(3:))//' && mkdir '//trim(folder(i))
           else
      	      command='if not exist '//trim(tempf)//sub//'nul'//' cd '//trim(opt_dir)//' && mkdir '//trim(folder(i))
           endif
           status= system(command)
         enddo
      endif

      if(OS=='linux') then
      	sub=sub2(1)
        do i=1,N_index
            tempf=trim(opt_dir)//trim(folder(i))
            command='if [ ! -d '//trim(tempf)//' ]; then cd '//trim(opt_dir)//'; mkdir '//trim(folder(i))//'; fi'
            status= system(command)
        enddo
      endif
    print*,' folders checked ...'
      
      return
      end subroutine check_folders

! to convert different time format into standard YYYYMMDD format.
! only accept limited input format now.
! by H.Yang @ CCRC, UNSW.
    subroutine time_convert(n,hours,string)
    use functions
    implicit none
    integer :: n, i,i0,i1,S_year,S_mon,S_day,t0,DoYs !,leapyear
    integer :: hours(n),days(13,2),indx(1)
    character*120 :: string
    data days/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365,  &
              0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366/

    if(size(hours) /= n) stop 'Error in time_convert: different dimension lengths..'
    DoYs=365
    string=adjustl(string)
    if(len_trim(string) == 0) return  ! the hours are in unit of YYYYMMDD already.
    if((string(1:6)=='day as') .or. (string(1:7)=='days as')) return
    if(string(1:4)=='hour')    hours=ceiling(hours/24.)        ! days

    if((string(1:3)/='day') .and. (string(1:4)/='hour')) then
      print*,'======= Congratulations ! ======'
      print*,'You have found a new unit for time dimensions, which is'
      print*,trim(string)
      print*,'Current code only deals with units like'
      print*,'day(s) as..., day(s) after/since..., hour(s) after/since...'
      print*,'please change the subroutine time_convert and try again...'
      stop ' now the code stops '  
    endif

    i=0   ! now to get the start date from attribute string
    do
       i=i+1
       if((string(i:i)>='0') .and. (string(i:i)<='9')) exit
    enddo
    string=string(i:)
    i1=0
    do
       i1=i1+1
       if(string(i1:i1) =='-') exit
    enddo
    i=i1-1
    read(string(1:i1-1), '(i<i>)' )  S_year

    string=string(i1+1:)
    i1=0
    do
       i1=i1+1
       if(string(i1:i1) =='-') exit
    enddo
    i=i1-1
    read(string(1:i1-1), '(i<i>)' )  S_mon

    string=string(i1+1:)
    i1=0
    do
       i1=i1+1
       if(string(i1:i1) ==' ') exit
    enddo
    i1=min(i1,3)
    i=i1-1
    read(string(1:i1-1), '(i<i>)' )  S_day

    t0=0
    if(S_year >= 1) then
      do i=0,S_year-1
        t0=t0+DoYs
        if(leapyear(i)==1) t0=t0+1
      enddo
    endif
	print*,t0
	print*,S_day
	print*,S_mon
	print*,S_year
    t0=t0+days(S_mon,leapyear(S_year)+1)+S_day-1  ! total days before the start date
	print*,'t0 is: ',t0
    hours=hours+t0   ! days after 0000/00/00

   do i=1,n
     i1=int(hours(i)/365.25)  ! fist find the start Year value
     i0=int(hours(i)/365.)
     t0=0
     S_year=i1
     if(S_year>1) then
       do i1=0,S_year-1
         t0=t0+DoYs
         if(leapyear(i1)==1) t0=t0+1
       enddo
     endif

     if(t0 == hours(i)) then
       S_year=S_year-1
       S_mon=12
       S_day=31
     else
       do
         S_mon=t0+days(13,leapyear(i1)+1)
         if(hours(i)>=t0 .and. hours(i)<=S_mon) exit
         t0=S_mon
         i1=i1+1
         S_year=i1
         if(i1 >i0 ) exit
       enddo

       i1=leapyear(S_year)+1    ! then find month value
       hours(i)=hours(i)-t0
       indx=minloc(abs(days(:,i1)-hours(i)))
       if(hours(i) <= days(indx(1),i1)) indx(1)=indx(1)-1
       S_mon=indx(1)
       S_day=hours(i)-days(S_mon,i1)    ! find day value
     endif

     hours(i)=S_year*10000+S_mon*100+S_day  ! YYYYMMDD

   enddo

!stop
    return
    end subroutine time_convert



      subroutine percentile(x, length, nl, per, oout)
      use COMM
      use functions
      integer :: length,nl
      real(DP),dimension(length):: x,xtos
      real(DP),dimension(nl)    :: per,oout
      real(DP):: bb,cc
      integer :: nn,i

! This sentence is no use, because per=(0.1,0.5,0.9) is fixed, and nl=3 also fixed
! maybe it's better to only put them here.... 
      if(maxval(per)>1..or.minval(per)<0.) stop 'ERROR in percentile: per out of range [0.,1.] !'

      nn=0
      do i=1, length
        if(nomiss(x(i)))then
          nn=nn+1
          xtos(nn)=x(i)
        endif
      enddo

      if(nn.eq.0) then
        oout=MISSING
      else
       call sort(nn,xtos)
!        call sort2(nn,xtos)  ! HY
        do i=1,nl
          bb=nn*per(i)+per(i)/3.+1./3.
          cc=real(int(bb))
          if(int(cc).ge.nn) then
            oout(i)=xtos(nn)
          else
            oout(i)=xtos(int(cc))+(bb-cc)*	 &
               (xtos(int(cc)+1)-xtos(int(cc)))
          endif
        enddo
      endif

      return
      end subroutine percentile




!---Sorts an array arr(1:n) into ascending numerical order using the Quicksort
!    algorithm. n is inpu; arr is replace on output by its sorted rearrangement.
!    Parameters: M is the size of subarrays sorted by straight insertion
!    and NSTACK is the required auxiliary.
      SUBROUTINE sort(n,arr)
	  use COMM, only:DP
      INTEGER, PARAMETER :: M=7,NSTACK=50
      INTEGER  :: n, i,ir,j,jstack,k,l,istack(NSTACK)
      REAL(DP) :: arr(n)
      REAL(DP) :: a,temp

      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END SUBROUTINE sort





      subroutine threshold(idata, lev, nl, odata, flg)
      use COMM
      use functions
      integer  :: flg,nl
      real(DP) :: idata(BYRS,DoY+2*SS),odata(DoY,nl), lev(nl)
      real(DP) :: tosort(BYRS*WINSIZE),rtmp(nl)
      integer  :: nn,i,j,k,Icritical,SS2

      Icritical=int(BYRS*WINSIZE*NoMissingThreshold)
      SS2=2*SS

      do i=1,DoY
        nn=0
        do j=1,BYRS
          do k=i,i+SS2
!            if(j.eq.1.and.k.eq.1) print*,'##2##',idata(j,k),MISSING
            if(nomiss(idata(j,k))) then
              nn=nn+1
              tosort(nn)=idata(j,k)
            endif
          enddo
        enddo
        if(nn.lt.Icritical) then
!          print*,"##1##",nn
          flg=1
          return
        endif
        call percentile(tosort,nn,nl,lev,rtmp)
          odata(i,:)=rtmp(:)
      enddo

      return
      end subroutine threshold




! ******* to get a free unit ID for I/O file.
	  subroutine getID(ID)
	  integer :: ID,ios
          logical :: lopen

	  ID=10
          do 
              inquire(unit=ID,opened=lopen,iostat=ios)
              if(ios==0.and.(.not.lopen)) return
              ID=ID+1
          enddo

	  return
	  end subroutine getID





	subroutine get_time(Iflag,ID,current)
	integer      :: tarr(8),ii,Iflag, ID
        character*23 :: file
        character*45 :: now(3),printf
        double precision :: current
        data now/'Task starts from:','Task ends at:','Current time:'/

        ii=3
        if(Iflag==0) ii=1
        if(Iflag==1) ii=2
	call Date_and_Time(values=tarr)
	write(file,12) tarr(1),tarr(2),tarr(3),tarr(5),tarr(6),tarr(7),tarr(8)
12      format(i4.4,'.',i2.2,'.',i2.2,' ',i2.2,':',i2.2,':',i2.2,'.',i3.3)

       printf=trim(now(ii))//' '//file

        current= tarr(5)*3600.d0+tarr(6)*60.d0+tarr(7)+0.001d0*tarr(8)

	  write(*,'(/,a,/)') printf
	  if(ID > 0)  write(ID,'(/,a,/)') printf

	return
	end subroutine get_time





subroutine err_handle(status,string)
 use netcdf
 integer       :: status
 character*(*) :: string

 if(status/=NF90_noerr) then
   print*,'Error in Netcdf function: '//trim(string)
   print*,'Error code: ',status
!   stop
 endif

 return
 end subroutine err_handle