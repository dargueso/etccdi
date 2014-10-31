      subroutine R95p(oout,thresOutpr,lon_i,lat_j)  ! r95out,r99out, prcpout
      use COMM
      use functions
      character*80:: ofile,file_thres
      integer     :: year, month, day, kth,cnt,leng,i,ID, Ky
      real(DP)    :: p95,p99
      real(DP),dimension(YRS,3):: oout
      real(DP),dimension(2):: thresOutpr
!      real(DP),dimension(YRS):: r95out,r99out, prcpout
      real(DP),dimension(TOT):: prcptmp
      real(DP),dimension(2)  :: rlev,rtmp

      if(Prcp_miss) return
      
      
      if (.not. is_thresfile) then
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
      else
        
        file_thres="thresholds.nc"
        call err_handle(nf90_open(file_thres,NF90_nowrite,ncidt), 'open file')
        call err_handle(NF90_INQ_VARID (ncidt, "thresanp95", ID_varnp95),'inquire p95 var ID')
        call err_handle(NF90_get_var(ncidt,ID_varnp95,p95,start=(/lon_i,lat_j/)),'get Var p95')
        
        call err_handle(NF90_INQ_VARID (ncidt, "thresanp99", ID_varnp99),'inquire p99 var ID')
        call err_handle(NF90_get_var(ncidt,ID_varnp99,p99,start=(/lon_i,lat_j/)),'get Var p99')
        
        call err_handle(nf90_close(ncidt), 'close file')
        
      end if  
      
      thresOutpr(1)=p95
      thresOutpr(2)=p99

      cnt=0

      oout=0.
      do i=1,YRS
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        do month=1,12
           kth=Mon(month,Ky)
          do day=1,kth
            cnt=cnt+1
            if(PRCP(cnt).ge.1..and.nomiss(PRCP(cnt))) then
              
              oout(i,3)=oout(i,3)+PRCP(cnt)

              if(nomiss(p95).and.PRCP(cnt).gt.p95) oout(i,1)=oout(i,1)+PRCP(cnt)
              if(nomiss(p99).and.PRCP(cnt).gt.p99) oout(i,2)=oout(i,2)+PRCP(cnt)
            endif
          enddo
        enddo
        if(YNASTAT(i,1).eq.1) then
          oout(i,:)=MISSING
        endif
      enddo

	 return
      end subroutine R95p 
