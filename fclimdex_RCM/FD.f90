      subroutine FD(oout) !FD, SU, ID, TR
      use COMM
      use functions

      integer      :: year, trno, kth, month,i,j,day, ID, Ky
      real         :: oout(YRS,4)
!  oout(,1)--FD, oout(,2)--SU, oout(,3)--ID, oout(,4)--TR      
!      data chrtmp/"FD","SU","ID","TR"/

      if(Tmax_miss .and. Tmin_miss) return

      trno=0
      oout=0
      do i=1,YRS
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        do month=1,12
          Kth=Mon(month,Ky)
          do day=1,kth
            trno=trno+1
            if(YMD(trno,3).ne.day) then
              print *, 'ERROR1 at FD!!!'
              stop
            endif
            if(nomiss(TMIN(trno)).and.TMIN(trno).lt.0) 	&
            !if(TMIN(trno).lt.0) 	&
               oout(i,1)=oout(i,1)+1
            if(nomiss(TMAX(trno)).and.TMAX(trno).gt.25) &
               oout(i,2)=oout(i,2)+1
            if(nomiss(TMAX(trno)).and.TMAX(trno).lt.0) 	&
               oout(i,3)=oout(i,3)+1
            if(nomiss(TMIN(trno)).and.TMIN(trno).gt.20) &
               oout(i,4)=oout(i,4)+1
          enddo
        enddo
      enddo

      do i=1,YRS
        if(YNASTAT(i,2)==1) then
                oout(i,2)=MISSING  ! SU
                oout(i,3)=MISSING  ! ID
        endif
        if(YNASTAT(i,3)==1) then
                oout(i,1)=MISSING  ! FD
                oout(i,4)=MISSING  ! TR
        endif
      enddo

      return
      end subroutine FD
