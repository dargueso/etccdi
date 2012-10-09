     subroutine CDD(ocddwd)  ! get yearly index: ocdd and ocwd
      use COMM
      use functions

!       character*80 :: ofile
      integer      :: year, month, day, kth, cnt,i, Ky !,ID
!      real         :: ocdd(YRS), ocwd(YRS), nncdd, nncwd
      real         :: ocddwd(YRS,2), nncdd, nncwd

      if(Prcp_miss) return   ! if all PRCP missing, ignore this step

      cnt=0
      ocddwd=0.
      do i=1,YRS
        if(i==1)nncdd=0.
        if(i==1)nncwd=0.
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        do month=1,12
             kth=mon(month,Ky)
          do day=1,kth
            cnt=cnt+1
            if(ismiss(PRCP(cnt))) then
              if(nncwd.gt.ocddwd(i,2)) ocddwd(i,2)=nncwd ! from Markus
              if(nncdd.gt.ocddwd(i,1)) ocddwd(i,1)=nncdd ! from Markus
              nncdd=0.
              nncwd=0.
            elseif(PRCP(cnt).lt.1) then
              nncdd=nncdd+1.
              if(nncwd.gt.ocddwd(i,2)) ocddwd(i,2)=nncwd
              nncwd=0.
            else
              nncwd=nncwd+1.
              if(nncdd.gt.ocddwd(i,1)) ocddwd(i,1)=nncdd
              nncdd=0.
            endif
!            if(year.eq.1959.and.month.eq.12) then 
!                    print *, month,day,nncdd, ocdd(i)
!            endif
          enddo
        enddo

        if(ocddwd(i,2).lt.nncwd) then
          if(year.eq.EYEAR) then
                  ocddwd(i,2)=nncwd
          elseif(PRCP(cnt+1).lt.1..or.ismiss(PRCP(cnt+1)))then
                  ocddwd(i,2)=nncwd
          endif
        endif

        if(ocddwd(i,1).lt.nncdd) then
          if(year.eq.EYEAR) then
                  ocddwd(i,1)=nncdd
          elseif(PRCP(cnt+1).ge.1..or.ismiss(PRCP(cnt+1)))then
                  ocddwd(i,1)=nncdd
          endif
          if(ocddwd(i,1).eq.0) ocddwd(i,1)=MISSING
        endif

        if(YNASTAT(i,1).eq.1) then
          ocddwd(i,1)=MISSING
          ocddwd(i,2)=MISSING
        endif
      enddo

      return
      end subroutine CDD 
