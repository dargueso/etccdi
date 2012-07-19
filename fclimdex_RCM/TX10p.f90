      subroutine TX10p(Tout,wcsdi,thresOut)  ! tn10out,tn50out,tn90out,tx10out,tx50out,tx90out (13 months); wsdi,csdi
      use COMM
      use functions
       character*80:: ofile
       integer:: year, month, day, kth, cnt, nn,  missxcnt, missncnt,	&
             iter, cntx, cntn,i,byear,flgtn,flgtx,flg,j,ID, Ky !,idum

       real::Tout(YRS,13,6),wcsdi(YRS,2)
       real,dimension(DoY,BYRS,BYRS-1) :: thresbx50,thresbn50,thresbn10,thresbn90,thresbx90,thresbx10
       real,dimension(DoY,3)           :: threstmp
       real,dimension(DoY,6)           :: thresOut
       real,dimension(BYRS,DoY+2*SS)   :: txdata,tndata,tnboot,txboot
       real,dimension(DoY)             :: thresan10,thresan90,thresax10,thresax90,thresax50,thresan50
       real,dimension(BYRS,DoY)        :: txdtmp,tndtmp
       real,dimension(YRS,13)          :: tn50out,tx10out,tx50out,tx90out,tn10out,tn90out
       real,dimension(YRS)             :: wsdi,csdi
       real,dimension(3)               :: rlevs
       real :: BYRSm,rDoY,temp

       integer :: SmByear ! SYear-BaseYear

      data rlevs/0.1,0.5,0.9/

      if(Tmax_miss .and. Tmin_miss) return

      BYRSm=BYRS-1.0
      rDoY=100.0/DoY
      cnt=0
      nn=0
      txdtmp=MISSING
      tndtmp=MISSING

      SmByear=SYEAR-BASESYEAR
      do i=1,YRS
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        nn=0
        do month=1,12
          Kth=Mon(month,Ky)
          do day=1,kth
            cnt=cnt+1
            if(year.ge.BASESYEAR.and.year.le.BASEEYEAR.and.(month.ne.2 &
              .or.day.ne.29))then
              nn=nn+1
              txdtmp(i+SmByear,nn)=TMAX(cnt)
              tndtmp(i+SmByear,nn)=TMIN(cnt)
            endif
          enddo
        enddo
        if(year.ge.BASESYEAR.and.year.le.BASEEYEAR.and.nn.ne.DoY)then
          print *,"date count error in TX10p!", nn
          stop
        endif
      enddo

            tndata(1,1:SS)=tndtmp(1,1)  ! i=1
            txdata(1,1:SS)=txdtmp(1,1)
            tndata(2:BYRS,1:SS)=tndtmp(1:BYRS-1,DoY+1-SS:DoY+1) ! i /=1
            txdata(2:BYRS,1:SS)=txdtmp(1:BYRS-1,DoY+1-SS:DoY+1)
            tndata(:,1+SS:DoY+SS)=tndtmp(:,1:DoY)
            txdata(:,1+SS:DoY+SS)=txdtmp(:,1:DoY)
            tndata(BYRS,1+DoY+SS:DoY+SS*2)=tndtmp(BYRS,DoY)
            txdata(BYRS,1+DoY+SS:DoY+SS*2)=txdtmp(BYRS,DoY)
            tndata(1:BYRS-1,1+DoY+SS:DoY+SS*2)=tndtmp(2:BYRS,1:SS)
            txdata(1:BYRS-1,1+DoY+SS:DoY+SS*2)=txdtmp(2:BYRS,1:SS)

      flgtn=0
      flgtx=0
!      call threshold(tndata,.1,thresan10, flgtn)
      call threshold(tndata,rlevs,3,threstmp,flgtn)
        thresan10(:)=threstmp(:,1)-1e-5
        thresan50(:)=threstmp(:,2)+1e-5
        thresan90(:)=threstmp(:,3)+1e-5
!      thresan10=thresan10-1e-5
      if(flgtn.eq.1) then
!        write(6,*) "TMIN Missing value overflow in exceedance rate"
        tn10out=MISSING
        tn50out=MISSING
        tn90out=MISSING
!      else
!        call threshold(tndata,.5,thresan50, flgtn)
!        thresan50=thresan50+1e-5
!        call threshold(tndata,.9,thresan90, flgtn)
!        thresan90=thresan90+1e-5
      endif

!MD      call threshold(txdata,.1,thresax10, flgtx)
      call threshold(txdata,rlevs,3,threstmp,flgtx)
        thresax10(:)=threstmp(:,1)-1e-5
        thresax50(:)=threstmp(:,2)+1e-5
        thresax90(:)=threstmp(:,3)+1e-5
!      thresax10=thresax10-1e-5
      if(flgtx.eq.1) then
 !       write(6,*) "TMAX Missing value overflow in exceedance rate"
        tx10out=MISSING
        tx50out=MISSING
        tx90out=MISSING
!      else
!        call threshold(txdata,.5,thresax50, flgtx)
!        thresax50=thresax50+1e-5
!        call threshold(txdata,.9,thresax90, flgtx)
!        thresax90=thresax90+1e-5
      endif

     thresOut(:,1)=thresan10
     thresOut(:,2)=thresan50
     thresOut(:,3)=thresan90
     thresOut(:,4)=thresax10
     thresOut(:,5)=thresax50
     thresOut(:,6)=thresax90

!$OMP parallel do default(shared) private(i,nn,txboot,tnboot,iter,threstmp,flg)
      do i=1,BYRS   ! This part consumes most of the time, due to "threshold"...
        txboot=txdata
        tnboot=tndata
        nn=0
        do iter=1,BYRS
          if(iter.ne.i) then
            nn=nn+1
              if(flgtx.eq.0) txboot(i,:)=txboot(iter,:)
              if(flgtn.eq.0) tnboot(i,:)=tnboot(iter,:)
            if(flgtx.eq.0)then
              call threshold(txboot,rlevs,3,threstmp,flg)
                thresbx90(:,i,nn)=threstmp(:,3)+1e-5
                thresbx50(:,i,nn)=threstmp(:,2)+1e-5
                thresbx10(:,i,nn)=threstmp(:,1)-1e-5
            endif

            if(flgtn.eq.0) then
              call threshold(tnboot,rlevs,3,threstmp,flg)
                thresbn90(:,i,nn)=threstmp(:,3)+1e-5
                thresbn50(:,i,nn)=threstmp(:,2)+1e-5
                thresbn10(:,i,nn)=threstmp(:,1)-1e-5
            endif
          endif
        enddo
      enddo
!$OMP end parallel do

      if(flgtx.eq.0)then
        tx10out=0.
        tx50out=0.
        tx90out=0.
      endif
      if(flgtn.eq.0)then
        tn10out=0.
        tn50out=0.
        tn90out=0.
      endif

      cnt=0
      do i=1,YRS
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        byear=year-BASESYEAR+1
        nn=0
        do month=1,12
          missncnt=0
          missxcnt=0
          Kth=Mon(month,Ky)
          do day=1,kth
            if(month.ne.2.or.day.ne.29) nn=nn+1
            cnt=cnt+1
            if(nomiss(TMAX(cnt)))then
              if(year.lt.BASESYEAR.or.year.gt.BASEEYEAR) then
                if(TMAX(cnt).gt.thresax90(nn)) tx90out(i,month)= tx90out(i,month)+1
                if(TMAX(cnt).gt.thresax50(nn)) tx50out(i,month)= tx50out(i,month)+1
                if(TMAX(cnt).lt.thresax10(nn)) tx10out(i,month)= tx10out(i,month)+1
              else
                do iter=1,BYRS-1
!                    if(byear.gt.30.or.iter.gt.29) then
!                      print *, i, year,month,day
!                      stop
!                    endif
                  if(TMAX(cnt).gt.thresbx90(nn,byear,iter))   tx90out(i,month)=tx90out(i,month)+1
                  if(TMAX(cnt).gt.thresbx50(nn,byear,iter))   tx50out(i,month)=tx50out(i,month)+1
                  if(TMAX(cnt).lt.thresbx10(nn,byear,iter))   tx10out(i,month)=tx10out(i,month)+1
                enddo
              endif
            else
              missxcnt=missxcnt+1
            endif
            if(nomiss(TMIN(cnt)))then
              if(year.lt.BASESYEAR.or.year.gt.BASEEYEAR) then
                if(TMIN(cnt).gt.thresan90(nn)) tn90out(i,month)= tn90out(i,month)+1
                if(TMIN(cnt).gt.thresan50(nn)) tn50out(i,month)= tn50out(i,month)+1
                if(TMIN(cnt).lt.thresan10(nn)) tn10out(i,month)= tn10out(i,month)+1
              else
                do iter=1,BYRS-1
                  if(TMIN(cnt).gt.thresbn90(nn,byear,iter))  tn90out(i,month)=tn90out(i,month)+1
                  if(TMIN(cnt).gt.thresbn50(nn,byear,iter))  tn50out(i,month)=tn50out(i,month)+1
                  if(TMIN(cnt).lt.thresbn10(nn,byear,iter))  tn10out(i,month)=tn10out(i,month)+1
                enddo
              endif
            else
              missncnt=missncnt+1
            endif
          enddo ! do day=1,kth

!          if(year.ge.BASESYEAR.and.year.le.BASEEYEAR)then
!            print *, year,month,tx10out(i,month),tx90out(i,month),	&
!             tn10out(i,month),tn90out(i,month),missxcnt,missncnt
!          endif

          if(year.ge.BASESYEAR.and.year.le.BASEEYEAR)then
            tn90out(i,month)=tn90out(i,month)/BYRSm
            tn50out(i,month)=tn50out(i,month)/BYRSm
            tn10out(i,month)=tn10out(i,month)/BYRSm
            tx90out(i,month)=tx90out(i,month)/BYRSm
            tx50out(i,month)=tx50out(i,month)/BYRSm
            tx10out(i,month)=tx10out(i,month)/BYRSm
          endif

!          if(year.eq.1952) then
!            print *,year,month,tn10out(i,month),tn90out(i,month)
!          endif

          if(missxcnt.le.10.and.flgtx.eq.0)then
            temp=100./(kth-missxcnt)
            tx90out(i,13)=tx90out(i,13)+tx90out(i,month)
            tx90out(i,month)=tx90out(i,month)*temp
            tx50out(i,13)=tx50out(i,13)+tx50out(i,month)
            tx50out(i,month)=tx50out(i,month)*temp
            tx10out(i,13)=tx10out(i,13)+tx10out(i,month)
            tx10out(i,month)=tx10out(i,month)*temp
          else
            tx90out(i,month)=MISSING
            tx50out(i,month)=MISSING
            tx10out(i,month)=MISSING
          endif
          if(missncnt.le.10.and.flgtn.eq.0)then
            temp=100./(kth-missncnt)
            tn90out(i,13)=tn90out(i,13)+tn90out(i,month)
            tn90out(i,month)=tn90out(i,month)*temp
            tn50out(i,13)=tn50out(i,13)+tn50out(i,month)
            tn50out(i,month)=tn50out(i,month)*temp
            tn10out(i,13)=tn10out(i,13)+tn10out(i,month)
            tn10out(i,month)=tn10out(i,month)*temp
          else
            tn90out(i,month)=MISSING
            tn50out(i,month)=MISSING
            tn10out(i,month)=MISSING
          endif
        enddo ! do month=1,12
        if(YNASTAT(i,3).eq.1.or.flgtn.eq.1) then
          tn10out(i,13)=MISSING
          tn50out(i,13)=MISSING
          tn90out(i,13)=MISSING
        else
          tn10out(i,13)=tn10out(i,13)*rDoY
          tn50out(i,13)=tn50out(i,13)*rDoY
          tn90out(i,13)=tn90out(i,13)*rDoY
        endif
        if(YNASTAT(i,2).eq.1.or.flgtx.eq.1) then
          tx10out(i,13)=MISSING
          tx50out(i,13)=MISSING
          tx90out(i,13)=MISSING
        else
          tx10out(i,13)=tx10out(i,13)*rDoY
          tx50out(i,13)=tx50out(i,13)*rDoY
          tx90out(i,13)=tx90out(i,13)*rDoY
        endif
      enddo

130   continue

      cnt=0
      wsdi=0.
      csdi=0.
!      if(flg.eq.1) then
!        wsdi=MISSING
!        csdi=MISSING
!        goto 140
!      endif

      do i=1,YRS
        cntx=0
        cntn=0
        nn=0
        year=i+SYEAR-1
        Ky=leapyear(year)+1
        do month=1,12
          Kth=Mon(month,Ky)
          do day=1,kth
            if(month.ne.2.or.day.ne.29) nn=nn+1
            cnt=cnt+1
            if(TMAX(cnt).gt.thresax90(nn).and.nomiss(TMAX(cnt))) then
              cntx=cntx+1
              if(month.eq.12.and.day.eq.31.and.cntx.ge.6)   wsdi(i)=wsdi(i)+cntx
            elseif(cntx.ge.6)then
              wsdi(i)=wsdi(i)+cntx
              cntx=0
            else
              cntx=0
            endif
            if(TMIN(cnt).lt.thresan10(nn).and.nomiss(TMIN(cnt))) then
              cntn=cntn+1
              if(month.eq.12.and.day.eq.31.and.cntn.ge.6)   csdi(i)=csdi(i)+cntn
            elseif(cntn.ge.6)then
              csdi(i)=csdi(i)+cntn
              cntn=0
            else
              cntn=0
            endif
          enddo  ! day
        enddo    ! month
        if(YNASTAT(i,3).eq.1) csdi(i)=MISSING
        if(YNASTAT(i,2).eq.1) wsdi(i)=MISSING
      enddo      ! year

140   continue

     Tout(:,:,1)=tn10out(:,:)
     Tout(:,:,2)=tn50out(:,:)
     Tout(:,:,3)=tn90out(:,:)
     Tout(:,:,4)=tx10out(:,:)
     Tout(:,:,5)=tx50out(:,:)
     Tout(:,:,6)=tx90out(:,:)

      wcsdi(:,1)=wsdi(:)
      wcsdi(:,2)=csdi(:)

 !     if(flgtx.eq.0)then
  !      what to do if flgtx NE 0 ?  ! H.Yang @ 2011.7.5
 !     endif

      return
      end subroutine TX10p 