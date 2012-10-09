      subroutine GSL(oout)
      use COMM
      use functions
      character*80 :: ofile
      integer      :: year,cnt,kth,month,day,marks,marke,i,ID, Ky
      real         :: TG,oout(YRS),strt(YRS),ee(YRS)  ! strt and ee may be scalar only ! 

      if(Tmax_miss .and. Tmin_miss) return

      strt=MISSING
      ee=MISSING
      cnt=0

      do i=1,YRS
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        marks=0
        marke=0
        do month=1,6
           kth=Mon(month,Ky)
          do day=1,kth
            cnt=cnt+1
            if(YMD(cnt,1)*10000+YMD(cnt,2)*100+YMD(cnt,3).ne.  &  ! Do we still need this???
              year*10000+month*100+day) then
              print*, 'date count ERROR in GSL!'
              print*, YMD(cnt,1)*10000+YMD(cnt,2)*100+YMD(cnt,3), &
                     year*10000+month*100+day
              stop
            endif
            if(nomiss(TMAX(cnt)).and.nomiss(TMIN(cnt))) then
              TG=(TMAX(cnt)+TMIN(cnt))/2.
            else
              TG=MISSING
            endif
            if(LATITUDE.lt.0) then
              if(nomiss(TG).and.TG.lt.5.)then
                marke=marke+1
              else
                marke=0
              endif
              if(marke.ge.6.and.i.gt.1.and.ismiss(ee(i-1)))then
!                 ee(i-1)=cnt-5  !correction Jana to get up to 365/366days
                 ee(i-1)=cnt-6
              endif
              if(ismiss(ee(i-1)).and.month.eq.6.and.day.eq.kth) then !Growing season never end
                ee(i-1)=cnt
              endif
            else
              if(nomiss(TG).and.TG.gt.5.)then
                marks=marks+1
              else
                marks=0
              endif
              if(marks.ge.6.and.ismiss(strt(i)))then
!                strt(i)=cnt-5  !correction Jana to get up to 365/366days
                strt(i)=cnt-6
              endif
            endif
          enddo
        enddo 
!MD: this part not clear??
!!        if(LATITUDE.lt.0.and.i.gt.1) then
!!          if(ismiss(ee(i-1)).and.nomiss(strt(i-1))) then
!!            ee(i-1)=cnt
!!          endif
!!        endif

        marks=0
        marke=0
        do month=7,12
          do day=1,MON(month,1) ! Here it doesn't matter whether it's leap year. 
            cnt=cnt+1
            if(nomiss(TMAX(cnt)).and.nomiss(TMIN(cnt))) then
              TG=(TMAX(cnt)+TMIN(cnt))/2.
            else
              TG=MISSING
            endif
            if(LATITUDE.lt.0) then
              if(nomiss(TG).and.TG.gt.5.)then
                marks=marks+1
              else
                marks=0
              endif
              if(marks.ge.6.and.ismiss(strt(i)))then
!                strt(i)=cnt-5 !correction Jana to get up to 365/366days
                strt(i)=cnt-6
              endif
            else
              if(nomiss(TG).and.TG.lt.5.)then
                marke=marke+1
              else
                marke=0
              endif
              if(marke.ge.6.and.ismiss(ee(i)))then
!                ee(i)=cnt-5 !correction Jana to get up to 365/366days
                ee(i)=cnt-6
              endif
              if (ismiss(ee(i)).and.month.eq.12.and.day.eq.MON(12,1)) then
                 ee(i)=cnt
              endif
            endif
          enddo
        enddo
!MD: this part not clear??
!!        if(ismiss(ee(i)).and.nomiss(strt(i))) then
!!          ee(i)=cnt
!!        endif
      enddo

      do i=1,YRS

        if(nomiss(strt(i)).and.nomiss(ee(i)).and.	&  ! get the result - GSL index.
          YNASTAT(i,2).ne.1.and.YNASTAT(i,3).ne.1)then
          oout(i)=ee(i)-strt(i)
        elseif(ismiss(strt(i)).or.ismiss(ee(i))) then
          oout(i)=0.
        endif
        if(YNASTAT(i,2).eq.1.or.YNASTAT(i,3).eq.1) oout(i)=MISSING

      enddo

      return
      end subroutine GSL 
