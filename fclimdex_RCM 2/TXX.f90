      subroutine TXX(oout) ! TXx,TXn,TNx,TNn,DTR (13 months)
      use COMM
      use functions
      character*80 :: ofile
 !     character*3  :: chrtmp(4)

      integer :: year,month,day,kth,cnt,nn,i,k,j,ID, Ky
      real    :: oout(YRS,13,5),yout(YRS,4), dtr(YRS,13)

!      data chrtmp/"TXx","TXn","TNx","TNn"/

      if(Tmax_miss .and. Tmin_miss) return

      oout=MISSING
      dtr=0.
      cnt=0
      do i=1,YRS
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        do month=1,12
          Kth=Mon(month,Ky)
          nn=0
          do day=1,kth
            cnt=cnt+1
            if(nomiss(TMAX(cnt)).and.nomiss(TMIN(cnt))) then
              dtr(i,month)=dtr(i,month)+(TMAX(cnt)-TMIN(cnt))
              nn=nn+1
            endif
            if(nomiss(TMAX(cnt)).and.(ismiss(oout(i,month,1)).or.  &
              TMAX(cnt).gt.oout(i,month,1))) then
              oout(i,month,1)=TMAX(cnt) ! TXX
            endif
            if(nomiss(TMAX(cnt)).and.(ismiss(oout(i,month,2)).or.  &
              TMAX(cnt).lt.oout(i,month,2))) then
              oout(i,month,2)=TMAX(cnt) ! TXN
            endif
            if(nomiss(TMIN(cnt)).and.(ismiss(oout(i,month,3)).or.  &
              TMIN(cnt).gt.oout(i,month,3))) then
              oout(i,month,3)=TMIN(cnt) ! TNX
            endif
            if(nomiss(TMIN(cnt)).and.(ismiss(oout(i,month,4)).or.  &
              TMIN(cnt).lt.oout(i,month,4))) then
              oout(i,month,4)=TMIN(cnt) ! TNN
            endif
          enddo 
          if(nn.gt.0.and.MNASTAT(i,month,2).eq.0.and.	  &
               MNASTAT(i,month,3).eq.0) then
            dtr(i,month)=dtr(i,month)/nn
          else
            dtr(i,month)=MISSING
          endif
          if(MNASTAT(i,month,2).eq.1)then
            oout(i,month,1)=MISSING
            oout(i,month,2)=MISSING
          endif
          if(MNASTAT(i,month,3).eq.1)then
            oout(i,month,3)=MISSING
            oout(i,month,4)=MISSING
          endif
        enddo
      enddo

      yout=MISSING
      do i=1,YRS
        nn=0
        do month=1,12
          if(nomiss(oout(i,month,1)).and.(ismiss(yout(i,1)).or.	 &
            oout(i,month,1).gt.yout(i,1))) then
            yout(i,1)=oout(i,month,1)
          endif
          if(nomiss(oout(i,month,2)).and.(ismiss(yout(i,2)).or.	 &
            oout(i,month,2).lt.yout(i,2))) then
            yout(i,2)=oout(i,month,2)
          endif
          if(nomiss(oout(i,month,3)).and.(ismiss(yout(i,3)).or.	 &
            oout(i,month,3).gt.yout(i,3))) then
            yout(i,3)=oout(i,month,3)
          endif
          if(nomiss(oout(i,month,4)).and.(ismiss(yout(i,4)).or.	 &
            oout(i,month,4).lt.yout(i,4))) then
            yout(i,4)=oout(i,month,4)
          endif
          if(nomiss(dtr(i,month))) then
            dtr(i,13)=dtr(i,13)+dtr(i,month)
            nn=nn+1
          endif
        enddo
        if(nn.gt.0.and.YNASTAT(i,2).eq.0.and.YNASTAT(i,3).eq.0) then
          dtr(i,13)=dtr(i,13)/nn
        else
          dtr(i,13)=MISSING
        endif
        if(YNASTAT(i,2).eq.1) then
          yout(i,1)=MISSING
          yout(i,2)=MISSING
        endif
        if(YNASTAT(i,3).eq.1) then
          yout(i,3)=MISSING
          yout(i,4)=MISSING
        endif
      enddo

     oout(:,13,1:4)=yout(:,:)
     oout(:,:,5)=dtr(:,:)

      return
      end subroutine TXX