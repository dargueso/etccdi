      subroutine R95p(oout)  ! r95out,r99out, prcpout
      use COMM
      use functions
      character*80:: ofile
      integer     :: year, month, day, kth,cnt,leng,i,ID, Ky
      real(DP)    :: p95,p99
      real(DP),dimension(YRS,3):: oout
!      real(DP),dimension(YRS):: r95out,r99out, prcpout
      real(DP),dimension(TOT):: prcptmp
      real(DP),dimension(2)  :: rlev,rtmp

      if(Prcp_miss) return

      cnt=0
      leng=0
      prcptmp=MISSING
      do i=1,YRS
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        do month=1,12
          kth=Mon(month,Ky)
          do day=1,kth
            cnt=cnt+1
            if(year.ge.BASESYEAR.and.year.le.BASEEYEAR.and.	&
              nomiss(PRCP(cnt)).and.PRCP(cnt).ge.1.)then
              leng=leng+1
              prcptmp(leng)=PRCP(cnt)
            endif
          enddo
        enddo
      enddo

      rlev(1)=0.95
      rlev(2)=0.99
      call percentile(prcptmp,leng,2,rlev,rtmp)
      p95=rtmp(1)
      p99=rtmp(2)
!      p95=percentile(prcptmp,leng,0.95)
!      p99=percentile(prcptmp,leng,0.99)

      cnt=0
!      r95out=0.
!      r99out=0.
!      prcpout=0.
      oout=0.
      do i=1,YRS
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        do month=1,12
           kth=Mon(month,Ky)
          do day=1,kth
            cnt=cnt+1
            if(PRCP(cnt).ge.1..and.nomiss(PRCP(cnt)))then
!              prcpout(i)=prcpout(i)+PRCP(cnt)
              oout(i,3)=oout(i,3)+PRCP(cnt)
!!              if(PRCP(cnt).gt.p95) oout(i,1)=oout(i,1)+PRCP(cnt)
!!              if(PRCP(cnt).gt.p99) oout(i,2)=oout(i,2)+PRCP(cnt)
              if(nomiss(p95).and.PRCP(cnt).gt.p95) oout(i,1)=oout(i,1)+PRCP(cnt)
              if(nomiss(p99).and.PRCP(cnt).gt.p99) oout(i,2)=oout(i,2)+PRCP(cnt)
            endif
          enddo
        enddo
        if(YNASTAT(i,1).eq.1) then
!          prcpout(i)=MISSING
!          r95out(i)=MISSING
!          r99out(i)=MISSING
          oout(i,:)=MISSING
        endif
      enddo

	 return
      end subroutine R95p 
