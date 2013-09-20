        subroutine meta_all(mode,nobj,xin,yout)
        use dimKrig
        implicit none
        integer, intent(in) :: nobj,mode
        double precision, dimension(ndim), intent(in)  :: xin
        double precision, dimension(nobj+1+3*(nfunc-1)), intent(out) :: yout
        double precision, dimension(nfunc) :: ycon
        integer :: i,k,imode,ict
        double precision :: yhat,RMSE,EI,pena,yhatprime(ndim)

        yout = 0.d0

        ict  = 0
        do k=1,nfunc-1
          imode = 0
          do i=1,nOPT
            if(nfOPT(i).eq.k)then
               imode = 1
               ict   = ict + 1
            end if
          end do
          call meta_call(k,imode,xin,yhat,yhatprime,RMSE,EI)
          if(imode.eq.1)then
            if(nobj.eq.nOPT*2)then
              yout((ict-1)*2+1) = EI
              yout((ict-1)*2+2) = yhat
            else if(nobj.eq.nOPT)then
              if(mode.eq.15)then
                yout((ict-1)*1+1) = RMSE
              else
                yout((ict-1)*1+1) = EI
              end if
            end if 
          end if
          ycon(k) = yhat
          yout(nobj+1+(k-1)*3+1) = yhat
          yout(nobj+1+(k-1)*3+2) = RMSE
          yout(nobj+1+(k-1)*3+3) = EI
        end do

        call constraint(ycon,pena)
        yout(nobj+1) = pena

        end subroutine meta_all
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine meta_call(k,mode,xin,yhat,yhatprime,RMSE,EI)
        use dimKrig
        implicit none
        integer, intent(in) :: k,mode
        double precision, dimension(ndim), intent(in) :: xin
        double precision, intent(out) :: yhat,RMSE,EI
        double precision, dimension(ndim), intent(out) :: yhatprime

        call meta_estimation(mode,xin,yhat,yhatprime,RMSE,EI,       &
             mEI,mLC,id_proc,ndim,kreg_orig,iscf,mode_dck,          &
             ccrf,                                                  &
             ipv(1,k),ipv(2,k),ipv(3,k),ipv(4,k),ipv(5,k),          &
             dak(iak( 1,k),k),dak(iak( 2,k),k),dak(iak( 3,k),k),    &
             dak(iak( 4,k),k),dak(iak( 5,k),k),dak(iak( 6,k),k),    &
             dak(iak( 7,k),k),dak(iak(10,k),k),dak(iak(11,k),k),    &
             dak(iak(12,k),k),dak(iak(13,k),k),dak(iak(14,k),k),    &
             dak(iak(15,k),k),dak(iak(16,k),k),                     &
             jbl(ibl( 1,k),k),jbl(ibl( 2,k),k),                     &
             ccm(icm( 1,k),k),                                      &
             ncr1,ncr2,dadd(1,k),dmul(1,k)  )

        end subroutine meta_call
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine meta_estimation(mode,xin,yhat,yhatprime,RMSE,EI,     &
                                   mEI,mLC,id_proc,ndim,kreg_orig,      &
                                   iscf,mode_dck,                       &
                                   ccrf,                                &
                                   rsample,nsize,kreg,lmax,tdim,        &
                                   train,mean,devi,deviratio,           &
                                   theta,power,Rinv,                    &
                                   RYFB,FR,FRFinv,ymax,ymin,dxadd,hvec, &
                                   ict_sample,nparent,                  &
                                   inf,                                 &
                                   ncr1,ncr2,dad,dmu)
        implicit none
! appear in correct, meta, search_krig, tool
! mode=0 : only for yhat
! mode=1 : yhat,RMSE,EI
! mode=2 : only for yhat and yhatprime
! else   : yhat,yhatprime,RMSE,EI
        integer, intent(in) :: mode,id_proc,mEI,mLC
        double precision, dimension(ndim), intent(in) :: xin
        double precision, intent(out) :: yhat,RMSE,EI

        double precision, dimension(ndim), intent(out) :: yhatprime

        integer,          intent(in) :: ndim,kreg,kreg_orig
        integer,          intent(in) :: rsample,nsize,iscf,mode_dck
        integer,          intent(in) :: lmax,tdim
        double precision, intent(in) :: devi,ymin,ymax,ccrf
        double precision, dimension(lmax),           intent(in) :: deviratio
        double precision, dimension(rsample,ndim),   intent(in) :: train
        double precision, dimension(ndim),           intent(in) :: theta,power
        double precision, dimension(kreg,1),         intent(in) :: mean
        double precision, dimension(nsize,1),        intent(in) :: RYFB
        double precision, dimension(kreg,nsize),     intent(in) :: FR
        double precision, dimension(kreg,kreg),      intent(in) :: FRFinv
        double precision, dimension(nsize,nsize),    intent(in) :: Rinv
        double precision, dimension(rsample),        intent(in) :: dxadd
        integer,          dimension(0:10,0:10),      intent(in) :: ict_sample
        integer,          dimension(rsample),        intent(in) :: nparent
        character(len=6), dimension(rsample),        intent(in) :: inf
        integer,                                     intent(in) :: ncr1,ncr2
        double precision, dimension(ncr1,ncr2),      intent(in) :: dad,dmu
        double precision, dimension(rsample,ndim,2), intent(in) :: hvec

        double precision, dimension(kreg_orig)   :: yout
        double precision, dimension(1,kreg)      :: freg
        double precision, dimension(1,1)         :: y1,y2
        double precision, dimension(kreg,1)      :: U
        double precision, dimension(1,kreg)      :: UFRF
        double precision, dimension(1,1)         :: rRr,UFRFU
        double precision, dimension(1,nsize)     :: rRinv
        double precision, dimension(nsize,1)     :: r,x
        double precision, dimension(nsize,kreg)  :: RFFRF
        double precision, dimension(1,1)         :: xr,myu,FRr,FRr1

        integer :: i,j
        double precision :: eps,z,cdf,pdf,ycor
        eps  = -1.d-3
        yhat =  0.d0
        RMSE =  0.d0
        EI   =  0.d0

        call makeMatrixPetitR(0,0,0, &
             ndim,nsize,rsample,iscf,mode_dck, &
             lmax,tdim,ict_sample, &
             train,xin,nparent,dxadd,inf,hvec, &
             theta,power,ccrf,r)
        call get_regression(0,0,0,xin,yout)
        freg = 0.d0
        do i=1,kreg_orig
          freg(1,i) = yout(i)
        end do

! yhat
        y1   = matmul(freg,mean)               ! (1*k)x(k*1)
        y2   = matmul(transpose(r),RYFB)       ! (1*ns)x(ns*1)
        yhat = y1(1,1) + y2(1,1)
!       write(*,'(3e15.5)')y1(1,1),y2(1,1),yhat
!       write(*,'(99f8.3)')(r(i,1),i=1,nsize)

! yhatprime

        if (mode .ge. 2) then

           do i=1,ndim
              call makeMatrixPetitR(1,i,0, &
                   ndim,nsize,rsample,iscf,mode_dck, &
                   lmax,tdim,ict_sample, &
                   train,xin,nparent,dxadd,inf,hvec, &
                   theta,power,ccrf,r)
!!$              call get_regression(1,i,0,xin,yout)
!!$              freg = 0.d0
!!$              do j=1,kreg_orig
!!$                 freg(1,j) = yout(j)
!!$              end do
!!$              y1 = matmul(freg,mean)
!!$              print *,i,y1(1,1)
              y2 = matmul(transpose(r),RYFB) 
              yhatprime(i)= y2(1,1) !+ y1(1,1)
           end do

!!$           write(*,*) xin,yhat,yhatprime
!!$           call calcdf(xin,ndim,8,yhatprime)
!!$           write(*,*) yhatprime

        end if

        if(mLC.ge.1) &
        call local_correction(xin,yhat,                                &
                              mEI,mLC,id_proc,ndim,kreg,iscf,mode_dck, &
                              rsample,nsize,                           &
                              train,mean,devi,theta,power,Rinv,        &
                              RYFB,FR,FRFinv,ymax,ymin,dxadd,          &
                              nparent,inf,ncr1,ncr2,dad,dmu )
        if(mode.eq.0 .or. mode .eq.2) return

        if(mEI.eq.0)then ! Real RMSE
          U      = matmul(FR,r)                  ! (k*ns)x(ns*1) FtR-1r
          U(:,1) = U(:,1) - freg(1,:)            ! (k*1) U=(FtR-1r-f)
          rRinv  = matmul(transpose(r),Rinv)     ! (1*ns)x(ns*ns)
          rRr    = matmul(rRinv,r)               ! (1*ns)x(ns*1)
          UFRF   = matmul(transpose(U),FRFinv)   ! (1*k)x(k*k)
          UFRFU  = matmul(UFRF,U)                ! (1*k)x(k*1)
          RMSE   = devi*( 1.d0 - rRr(1,1) + UFRFU(1,1) )
          if(RMSE.gt.0)then
             RMSE = dsqrt(RMSE)
          else if(RMSE.ge.eps.and.RMSE.le.0.d0)then
             RMSE = 0.d0
          else
!            if(id_proc.eq.0) &
!            write(*,'(11x,a,e12.3)')'>> Negative RMSE = ',RMSE
          end if
        else if(mEI.ge.1.and.mEI.le.4)then ! Quasi RMSE by Distance
          call Quasi_RMSE(mEI,ndim,rsample,xin,train,nparent,inf,ymin,ymax,RMSE)
        else
          stop'unknown mEI'
        end if
! EI
        if(RMSE.le.0.d0)then
          EI = 0.d0
        else
          z  = (ymin-yhat)/RMSE
          call get_pdfcdf(z,pdf,cdf)
          EI = (ymin-yhat)*cdf + RMSE*pdf
          if(EI.ge.eps.and.EI.le.0.d0)then
             EI = 0.d0
          else if(EI.lt.eps)then
           if(id_proc.eq.0)then
           write(*,'(11x,a,4e12.3)') &
           '>> Negative  EI  = ',EI,RMSE,yhat,ymin
           end if
          end if
        end if

        end subroutine meta_estimation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine Quasi_RMSE(mEI,ndim,rsample, &
                              xin,train,nparent,inf,ymin,ymax,uncty)
        implicit none
        integer,                                     intent(in) :: mEI
        integer,                                     intent(in) :: ndim,rsample
        double precision, dimension(ndim),           intent(in) :: xin
        double precision, dimension(rsample,ndim),   intent(in) :: train
        integer,          dimension(rsample),        intent(in) :: nparent
        double precision,                            intent(in) :: ymin,ymax
        character(len=6), dimension(rsample),        intent(in) :: inf
        double precision, intent(out) :: uncty

        integer :: i,k,ict
        double precision :: d,dd,d1min,d2min
        double precision :: dscale,fscale

        d1min = 1.d10
        d2min = 1.d11
        ict   = 0
        do 100 i=1,rsample
          if(nparent(i).ne.0)go to 100 ! only for real points
          if(inf(i)(1:2).ne.'1_')go to 100
          ict = ict + 1
          dd = 0.d0
          do 110 k=1,ndim
             d  = xin(k)-train(i,k)
             if(abs(d).ge.d2min)go to 100
             dd = dd + d**2
110       continue
          dd= dsqrt(dd)
          if(dd.lt.d1min)then
            d2min = d1min
            d1min = dd
          else if(dd.ge.d1min.and.dd.lt.d2min)then
            d2min = dd
          end if
100     continue

        if(mEI.eq.1.or.mEI.eq.2)then
          dscale = 1.d0
        else if(mEI.eq.3.or.mEI.eq.4)then
          dscale = (1.d0/dble(ict))**(1.d0/dble(ndim))
        end if
        if(mEI.eq.1.or.mEI.eq.3)then
          fscale = 1.d0
        else if(mEI.eq.2.or.mEI.eq.4)then
          fscale = dabs(ymax-ymin)*0.5d0
        end if

        d1min =  d1min / dscale
        d2min =  d2min / dscale
        uncty =  d1min * fscale * d2min 

!       write(*,*)d1min,fscale,d2min
        if(uncty.lt.0.d0)stop'uncertainty<0 in Quasi_RMSE'
        end subroutine Quasi_RMSE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine get_pdfcdf(zin,pdf,cdf)
        implicit none
        double precision, intent(in)  :: zin
        double precision, intent(out) :: pdf,cdf
        integer :: num,i
        double precision :: pi,zst,zed,xst,xed,x1,x2,dx,pdf1,pdf2
! zin E [-10:10] is enough
! zin = -inf, pdf = 0.d0, cdf = 0.d0
! zin =  0.0, pdf ~ 0.4,  cdf = 0.50
! zin =  inf, pdf = 0.d0, cdf = 1.d0

        pi  = 4.d0*datan(1.d0)
        pdf = 0.d0
        cdf = 0.d0

        zst = -10.d0
        zed = +10.d0
        num = 101
        if(zed.le.zst)stop'zst,zed in EI'

        if(zin.lt.zst)then
          pdf = 0.d0
          cdf = 0.d0
        else if(zin.gt.zed)then
          pdf = 0.d0
          cdf = 1.d0
        else
          pdf = exp(-0.5d0*zin*zin)/dsqrt(2.d0*pi)
          if(zin.le.0.d0)then
            cdf = 0.d0
            xst = zst
            xed = zin
          else
            cdf = 0.5d0
            xst = 0.d0
            xed = zin
          end if
          if(xed.le.xst)stop'xst,xed in EI'
          dx = (xed-xst)/dble(num-1)
          do i=1,num-1
            x1   = xst + dble(i-1)*dx
            x2   = xst + dble(i)*dx
            pdf1 = exp(-0.5d0*x1*x1)/dsqrt(2.d0*pi)
            pdf2 = exp(-0.5d0*x2*x2)/dsqrt(2.d0*pi)
            cdf  = cdf + (pdf1+pdf2)*dx*0.5d0
          end do
          if(dabs(pdf).le.1.d-10)pdf = 0.d0
          if(dabs(cdf).le.1.d-10)cdf = 0.d0
        end if

        end subroutine get_pdfcdf
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine meta_quadratic(xin,yhat,RMSE,EI)
        use dimKrig
        implicit none
        double precision, dimension(ndim), intent(in) :: xin
        double precision, intent(out) :: yhat,RMSE,EI
        integer :: i,j,id
        double precision, dimension(1,1)       :: d1y,d2y
        double precision, dimension(1,ndim)    :: dx
        double precision, dimension(ndim,1)    :: g
        double precision, dimension(ndim,ndim) :: h

        if(IDopt(1).lt.0.or.IDopt(1).gt.nsample)then
          stop'no optimal for meta_quadratic'
        end if
        if(info(IDopt(1)).ne.'1_FGH ')then
          stop'no hessian for optimal in meta_quadratic'
        end if
!       write(*,'(1x,a,i5,f15.8)')'>> Optimal',IDopt(1),func(IDopt(1),1)
!       call stop_all
        id = IDopt(1)

        do i=1,ndim
          dx(1,i) = xin(i) - sample(id,i)
           g(i,1) = gfunc(id,1,i)
           do j=1,ndim
             h(i,j) = hfunc(id,1,i,j)
           end do
        end do

        d1y  = matmul( dx,g )                       ! (1*n)x(n*1)
        d2y  = matmul( matmul(dx,h),transpose(dx) ) ! (1*n)x(n*n)x(n*1)
        yhat = func(id,1) + d1y(1,1) + 0.5d0*d2y(1,1)
        RMSE = 0.d0
        EI   = 0.d0

        end subroutine meta_quadratic
