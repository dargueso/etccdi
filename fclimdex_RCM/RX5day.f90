      subroutine RX5day(oout)  ! Rx1day, Rx5day (13 months)
      use COMM
      use functions
      character*80 :: ofile
      integer :: year, month, day,cnt,j,k,kth,i,ID, Ky
      real    :: oout(YRS,13,2),r1(YRS,13), r5(YRS,13), r5prcp

      if(Prcp_miss) return

      cnt=0
      r1=MISSING
      r5=MISSING
      do i=1,YRS
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        do month=1,12
          Kth=Mon(month,Ky)
          do day=1,kth
            cnt=cnt+1
!            if(year.eq.1904.and.month.eq.1) then
!              print *, year, month, day, PRCP(cnt)
!            endif
            if(cnt.gt.5)then
              r5prcp=0.
              do k=cnt-4,cnt
                if(nomiss(PRCP(k)))then
                  r5prcp=r5prcp+PRCP(k)
                endif
              enddo
            else
              r5prcp=MISSING
            endif
            if(nomiss(PRCP(cnt)).and.(ismiss(r1(i,month))	&
              .or.PRCP(cnt).gt.r1(i,month))) then
              r1(i,month)=PRCP(cnt)
            endif
            if(nomiss(PRCP(cnt)).and.r5prcp.gt.r5(i,month)) then
              r5(i,month)=r5prcp
            endif
          enddo
          if(MNASTAT(i,month,1).eq.1) then
            r1(i,month)=MISSING
            r5(i,month)=MISSING
          endif
          if(nomiss(r1(i,month)).and.(ismiss(r1(i,13))	  &
            .or.r1(i,month).gt.r1(i,13))) then
            r1(i,13)=r1(i,month)
          endif
          if(nomiss(r5(i,month)).and.(ismiss(r5(i,13))	  &
            .or.r5(i,month).gt.r5(i,13))) then
            r5(i,13)=r5(i,month)
          endif
        enddo
        if(YNASTAT(i,1).eq.1) then
          r1(i,13)=MISSING
          r5(i,13)=MISSING
        endif
      enddo
	  oout(:,:,1)=r1(:,:)
	  oout(:,:,2)=r5(:,:)

	  return
      end subroutine RX5day