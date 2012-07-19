!   (C) Copr. 1986-92 Numerical Recipes Software &#5,.

      subroutine qc(oout)  ! PRCPmiss,TMAXmiss,TMINmiss  (13 months)
      use COMM
      use functions
      integer  :: rno, i,j !,ios,ID,IDotmp,IDofil !, tmpymd(MaxYear*DoY,3), ith
      integer  :: kth,month,k,trno,ymiss(3),mmiss(3),	&  !,stdcnt(3),tmpcnt,now,now1
                   Ky
      real :: oout(YRS,13,3)

      rno=Ntime

      oout=0.
      trno=0
      MNASTAT=0  ! This was used globally
      YNASTAT=0  ! This was used globally
      do i=SYEAR,EYEAR
        Ky=leapyear(i)+1
        ymiss=0
   !     stdcnt=0
        do month=1,12
          mmiss=0
          Kth=Mon(month,Ky)
          do k=1,kth
            trno=trno+1

            if(TMAX(trno).lt.TMIN(trno).and.nomiss(TMAX(trno)) &
                        .and.nomiss(TMIN(trno))) then
              TMAX(trno)=MISSING
              TMIN(trno)=MISSING
   !           write(IDofil, *) i*10000+month*100+k, "TMAX<TMIN!!"
            endif

     !       if(month.ne.2.or.k.ne.29) then
     !         stdcnt=stdcnt+1
     !         stddata(stdcnt, i-SYEAR+1, 1)=PRCP(trno)
     !         stddata(stdcnt, i-SYEAR+1, 2)=TMAX(trno)
     !         stddata(stdcnt, i-SYEAR+1, 3)=TMIN(trno)
     !       endif

            if((TMAX(trno).lt.-100..or.TMAX(trno).gt.100.).and.  &
              nomiss(TMAX(trno))) then
              TMAX(trno)=MISSING
   !           write(IDofil, *) i*10000+month*100+k, "TMAX over bound!!"
            endif
            !!if((TMIN(trno).lt.-120..or.TMIN(trno).gt.100.).and.  &
            if((TMIN(trno).lt.-150..or.TMIN(trno).gt.100.).and.  &
              nomiss(TMIN(trno))) then
              TMIN(trno)=MISSING
   !           write(IDofil, *) i*10000+month*100+k, "TMIN over bound!!"
            endif
            if(PRCP(trno).lt.0.and.nomiss(PRCP(trno))) then
              PRCP(trno)=MISSING
   !           write(IDotmp,*) i*10000+month*100+k, "PRCP less then 0!!"
            endif
            if(ismiss(PRCP(trno))) then
              mmiss(1)=mmiss(1)+1
              ymiss(1)=ymiss(1)+1
            endif
            if(ismiss(TMAX(trno))) then
              mmiss(2)=mmiss(2)+1
              ymiss(2)=ymiss(2)+1
            endif
            if(ismiss(TMIN(trno))) then
              mmiss(3)=mmiss(3)+1
              ymiss(3)=ymiss(3)+1
            endif
          enddo

          do k=1,3
            oout(i-SYEAR+1,month,k)=mmiss(k)
            if (mmiss(k).gt.3) then
              MNASTAT(i-SYEAR+1,month,k)=1
            endif
          enddo
        enddo  ! month
        do k=1,3
          oout(i-SYEAR+1,13,k)=ymiss(k)
          if(ymiss(k).gt.15) then
            YNASTAT(i-SYEAR+1,k)=1
          endif
        enddo
      enddo


!      do i=1,YRS
!      print *,(YNASTAT(i,k),k=1,3)
!      enddo

!  Calculate STD for PRCP, TMAX and TMIN; then figure out outliers
!       stdval=0.   ! We don't use this part currently... H.Yang
!       m1=0.
!       do i=1,DoY
!         stdcnt=0
!         do j=1,YRS
!           do k=2,3
!             if(nomiss(stddata(i,j,k))) then
!               stdcnt(k)=stdcnt(k)+1
!               m1(i,k)=m1(i,k)+stddata(i,j,k)
!             endif
!           enddo
!         enddo
!         do k=2,3
!           if(stdcnt(k).gt.0) then
!             m1(i,k)=m1(i,k)/real(stdcnt(k))
!           endif
!         enddo
!         stdtmp=0.
!         do j=1,YRS
!           do k=2,3
!             if(stdcnt(k).gt.2.and.nomiss(stddata(i,j,k))) then
!                stdval(i,k)=stdval(i,k)+			&
!               (stddata(i,j,k)-m1(i,k))**2./(real(stdcnt(k))-1.)
!             endif
!           enddo
!         enddo
! 
!         do k=2,3
!           if(stdcnt(k).gt.2) then
!             stdval(i,k)=stdval(i,k)**0.5
!           else 
!             stdval(i,k)=MISSING
!           endif
!         enddo
!       enddo  ! i loop

!       trno=0  ! This part will be added later ...  H.Yang
!       do i=SYEAR,EYEAR
!         Ky=leapyear(i)+1
!         tmpcnt=0
!         do month=1,12
!           Kth=Mon(month,Ky)
!           do k=1, kth
!             trno=trno+1
!             if(month.ne.2.or.k.ne.29) tmpcnt=tmpcnt+1
!             if(nomiss(stdval(tmpcnt,2)))then
!               if(abs(TMAX(trno)-m1(tmpcnt,2)).gt.		  &
!                stdval(tmpcnt,2)*STDSPAN.and.nomiss(TMAX(trno)))   &
!            write(IDofil, *) "Outlier: ", i, month, k, "TMAX: ", TMAX(trno), &
!                  "Lower limit:",m1(tmpcnt,2)-stdval(tmpcnt,2)*STDSPAN,	 &
!                  "Upper limit:",m1(tmpcnt,2)+stdval(tmpcnt,2)*STDSPAN	 
!             endif
!             if(nomiss(stdval(tmpcnt,3)))then
!               if(abs(TMIN(trno)-m1(tmpcnt,3)).gt.				   &
!                stdval(tmpcnt,3)*STDSPAN.and.nomiss(TMIN(trno)))	   &
!            write(IDofil, *) "Outlier: ", i, month, k, "TMIN: ", TMIN(trno),	&
!                  "Lower limit:",m1(tmpcnt,3)-stdval(tmpcnt,3)*STDSPAN,	&
!                  "Upper limit:",m1(tmpcnt,3)+stdval(tmpcnt,3)*STDSPAN	
!             endif
!           enddo ! end do day
!         enddo !end do month
!       enddo ! end do year

!      open(20,file=trim(omissf))
!       call getID(ID)
!       open(ID,file=omissf)			! H.Yang
!       do i=SYEAR,EYEAR
!         do k=1,3
!           write(ID,'(i4,2x,a8,13i4)')	 &
!          i,title(k),(missout(i+1-SYEAR,j,k),j=1,13)
!         enddo
!       enddo
!       close(ID)
!       close(IDotmp)
!       close(IDofil)

!   QC  part finished, prepared data set: YMD(3), PRCP, TMAX & TMIN
!   and NASTAT dataset for missing values monthly and annual
      return
      end subroutine qc 
