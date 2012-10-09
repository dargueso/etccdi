     subroutine Rnnmm(oout)	  ! R10mm,R20mm,Rnnmm,sdii
      use COMM
      use functions, only : leapyear
      character*80 :: ofile
      integer      :: year,month,day,kth,cnt,nn,i,k,ID, Ky
      real         :: oout(YRS,4) !,sdii(YRS)

      if(Prcp_miss) return

      cnt=0
      oout=0.
!      sdii=0.
      do i=1,YRS
        nn=0
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        do month=1,12
          Kth=Mon(month,Ky)
          do day=1,kth
            cnt=cnt+1
            if(PRCP(cnt).ge.1.) then
              oout(i,4)=oout(i,4)+PRCP(cnt)
              nn=nn+1
            endif
            if(PRCP(cnt).ge.10.) oout(i,1)=oout(i,1)+1.
            if(PRCP(cnt).ge.20.) oout(i,2)=oout(i,2)+1.
            if(PRCP(cnt).ge.PRCPNN) oout(i,3)=oout(i,3)+1.
          enddo
        enddo
        if(nn.gt.0) then
          oout(i,4)=oout(i,4)/nn
        endif
      enddo

      do i=1,YRS
        if(YNASTAT(i,1).eq.1) then
           oout(i,:)=MISSING
        endif
      enddo

	  return
      end subroutine Rnnmm 