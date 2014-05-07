        subroutine calc_correction_term
        use dimKrig
        implicit none
        integer :: i,j,k,l
        double precision :: yhat,RMSE,EI
        double precision, dimension(ndim) :: xin,dy,yhatprime
        double precision, dimension(ndim,ndim) :: d2y
        integer :: ict,mode,inform
        integer :: iv,ig,id,ih
        double precision :: vv,vg,vd,vh

        ncr1 = 1+ndim+1+ndim+nhes+1+ndim+nhes
        ncr2 = nsample
        allocate( dadd(ncr2*ncr1,nfunc-1) )
        allocate( dmul(ncr2*ncr1,nfunc-1) )
        dadd  = 0.d0
        dmul  = 0.d0
        if(mLC.eq.0)return

        if(id_proc.eq.0) &
        write(*,*)'>> Make local correction Databese for metamodel'

        do 100 k=1,nfunc-1
        mode = 0
        do i=1,nCOK
          if(k.eq.nfCOK(i))mode = 1
        end do
        iv = 0
        ig = 0
        id = 0
        ih = 0
        vv = 0.d0
        vg = 0.d0
        vd = 0.d0
        vh = 0.d0
        do 110 i=1,nsample
          xin(:) = sample(i,:)
          call meta_estimation(0,xin,yhat,yhatprime,RMSE,EI,          &
!              mEI,mLC,id_proc,ndim,kreg_orig,iscf,mode_dck,          &
               mEI,  0,id_proc,ndim,kreg_orig,iscf,mode_dck,          &
               ccrf,                                                  &
               ipv(1,k),ipv(2,k),ipv(3,k),ipv(4,k),ipv(5,k),          &
               dak(iak( 1,k),k),dak(iak( 2,k),k),dak(iak( 3,k),k),    &
               dak(iak( 4,k),k),dak(iak( 5,k),k),dak(iak( 6,k),k),    &
               dak(iak( 7,k),k),dak(iak(10,k),k),dak(iak(11,k),k),    &
               dak(iak(12,k),k),dak(iak(13,k),k),dak(iak(14,k),k),    &
               dak(iak(15,k),k),dak(iak(16,k),k),                     &
               jbl(ibl( 1,k),k),jbl(ibl( 2,k),k),                     &
               ccm(icm( 1,k),k),                                      &
               ncr1,ncr2,dadd(1,k),dmul(1,k)  ) ! still dummy here
          call diff_estimation(xin,dy,d2y, &
               id_proc,ndim,kreg_orig,iscf,mode_dck,ccrf,          &
               ipv(1,k),ipv(2,k),ipv(3,k),ipv(4,k),ipv(5,k),       &
               dak(iak( 1,k),k),dak(iak( 2,k),k),dak(iak( 3,k),k), &
               dak(iak( 4,k),k),dak(iak( 5,k),k),dak(iak( 6,k),k), &
               dak(iak(10,k),k),dak(iak(15,k),k),dak(iak(16,k),k), &
               jbl(ibl( 1,k),k),jbl(ibl( 2,k),k),                  &
               ccm(icm( 1,k),k) )

          ! info
          if(func(i,nfunc).ne.0.d0)then
            dadd((i-1)*ncr1+1,k) = -1.d0
            dmul((i-1)*ncr1+1,k) = -1.d0
            go to 110
          else
            if(mode.eq.0)then
              dadd((i-1)*ncr1+1,k) = 0.d0
              dmul((i-1)*ncr1+1,k) = 0.d0
            else
              if(     info(i).eq.'1_F   ')then
                inform = 0
              else if(info(i).eq.'1_FG  ')then
                inform = 1
              else if(info(i).eq.'1_FGH ')then
                inform = 3
              else
                inform = 0
              end if
              dadd((i-1)*ncr1+1,k) = dble(inform)
              dmul((i-1)*ncr1+1,k) = dble(inform)
            end if
          end if
          ! dv
          do j=1,ndim
            dadd((i-1)*ncr1+1+j,                 k) = xin(j)
            dmul((i-1)*ncr1+1+j,                 k) = xin(j)
          end do
          ! func
          if(yhat.eq.0.d0)stop'multiplicative correlation with flo = 0'
          dadd(  (i-1)*ncr1+1+ndim+1,            k) = func(i,k) - yhat
          dmul(  (i-1)*ncr1+1+ndim+1,            k) = func(i,k) / yhat
          dadd(  (i-1)*ncr1+1+ndim+1+ndim+nhes+1,k) = func(i,k)
          dmul(  (i-1)*ncr1+1+ndim+1+ndim+nhes+1,k) = func(i,k)
          iv = iv + 1
          vv = vv + (func(i,k)-yhat)**2
          if(mode.eq.0)go to 110
          ! grad/hess
          ict = 0
          do j=1,ndim
            dadd((i-1)*ncr1+1+ndim+1+j,k) = gfunc(i,k,j) - dy(j)
            dmul((i-1)*ncr1+1+ndim+1+j,k) = gfunc(i,k,j) / yhat &
                                          - dy(j) * func(i,k) / (yhat**2)
            dadd((i-1)*ncr1+1+ndim+1+ndim+nhes+1+j,k) = gfunc(i,k,j)
            dmul((i-1)*ncr1+1+ndim+1+ndim+nhes+1+j,k) = gfunc(i,k,j)
            do l=j,ndim
               ict = ict + 1
               dadd((i-1)*ncr1+1+ndim+1+ndim+ict,k)           &
                    = hfunc(i,k,j,l) - d2y(j,l)
               dmul((i-1)*ncr1+1+ndim+1+ndim+ict,k)           &
                    = hfunc(i,k,j,l) / yhat                   &
                    - d2y(j,l) * func(i,k) / (yhat**2)        &
                    + 2.d0*func(i,k)*dy(j)*dy(l)/(yhat**3)    &
                    - dy(j)*gfunc(i,k,l) / (yhat**2)          &
                    - dy(l)*gfunc(i,k,j) / (yhat**2)
               dadd((i-1)*ncr1+1+ndim+1+ndim+nhes+1+ndim+ict,k) = hfunc(i,k,j,l)
               dmul((i-1)*ncr1+1+ndim+1+ndim+nhes+1+ndim+ict,k) = hfunc(i,k,j,l)
            end do
          end do
          ! rmse
          if(info(i).eq.'1_FG  '.or.info(i).eq.'1_FGH ')then
            do j=1,ndim
               ig = ig + 1
               vg = vg + (gfunc(i,k,j)-dy(j))**2
            end do
          end if
          if(info(i).eq.'1_FH  '.or.info(i).eq.'1_FGH ')then
            do j=1,ndim
            do l=j,ndim
              ih = ih + 1
              vh = vh + (hfunc(i,k,j,l)-d2y(j,l))**2
              if(j.eq.l)then
                id = id + 1
                vd = vd + (hfunc(i,k,j,l)-d2y(j,l))**2
              end if
            end do
            end do
          end if
110     continue ! sample loop
        vv = vv / dble(iv)
        vg = vg / dble(ig)
        vd = vd / dble(id)
        vh = vh / dble(ih)
        if(id_proc.eq.0)then
        if(ndebug.eq.1)then
          write(*,'(1x,a,i2)')'>> Check in calc_corr_term for Func-',k
          if(iv.ne.0) &
          write(*,'(6x,a,e15.5,a,i5,a)') &
          '>> RMSE for Func = ',vv,' among',iv,' data'
          if(ig.ne.0) &
          write(*,'(6x,a,e15.5,a,i5,a)') &
          '>> RMSE for Grad = ',vg,' among',ig,' data'
          if(id.ne.0) &
          write(*,'(6x,a,e15.5,a,i5,a)') &
          '>> RMSE for DHes = ',vd,' among',id,' data'
          if(ih.ne.0) &
          write(*,'(6x,a,e15.5,a,i5,a)') &
          '>> RMSE for Hess = ',vh,' among',ih,' data'
        end if
        end if
100     continue ! function loop

        if(ndebug.eq.1)then
          do i=1,ncr1*ncr2
             write(*,*)i,dadd(i,1)
          end do
          do i=1,ncr1*ncr2
             write(*,*)i,dmul(i,1)
          end do
          call stop_all
        end if
        end subroutine calc_correction_term
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine local_correction(xin,yhat, &
                   mEI,mLC,id_proc,ndim,kreg,iscf,mode_dck, &
                   rsample,nsize,                           &
                   train,mean,devi,theta,power,Rinv,        &
                   RYFB,FR,FRFinv,ymax,ymin,dxadd,          &
                   nparent,inf,ncr1,ncr2,dad,dmu )

        implicit none
        integer, intent(in) :: mEI,mLC,id_proc,ndim,kreg,iscf,mode_dck
        integer, intent(in) :: rsample,nsize,ncr1,ncr2
        double precision, dimension(ndim), intent(in) :: xin
        double precision, intent(inout) :: yhat

        double precision, intent(in) :: devi,ymin,ymax
        double precision, dimension(rsample,ndim),   intent(in) :: train
        double precision, dimension(ndim),           intent(in) :: theta,power
        double precision, dimension(kreg,1),         intent(in) :: mean
        double precision, dimension(nsize,1),        intent(in) :: RYFB
        double precision, dimension(kreg,nsize),     intent(in) :: FR
        double precision, dimension(kreg,kreg),      intent(in) :: FRFinv
        double precision, dimension(nsize,nsize),    intent(in) :: Rinv
        double precision, dimension(rsample),        intent(in) :: dxadd
        integer,          dimension(rsample),        intent(in) :: nparent
        character(len=4), dimension(rsample),        intent(in) :: inf
        double precision, dimension(ncr1,ncr2),      intent(in) :: dad,dmu

        integer :: i,k
        integer :: loc1,loc2,nhes
        double precision :: yhap,RMSE,EI,r_add,r_mul
        double precision :: dd,d1min,d2min,d1,d2,dp,co
        double precision :: yadd1,yadd2,ymul1,ymul2
        double precision :: yadp1,yadp2,ymup1,ymup2,yhip,gam
        double precision :: yad11,yad12,ymu11,ymu12
        double precision :: yad21,yad22,ymu21,ymu22
        double precision :: ydtc
        double precision, dimension(ndim) :: xp,yhatprime

! Find Target Sample
        loc1  = 0
        loc2  = 0
        d1min = 1.d11
        d2min = 1.d10
        do 100 i=1,ncr2 !nsample
          if(dad(1,i).lt.0.d0)go to 100
          dd = 0.d0
          do 110 k=1,ndim
             dd = dd + (xin(k)-dad(1+k,i))**2
             if(dd.ge.d2min)go to 100
110       continue
          if(dd.lt.d1min)then
            d2min = d1min
            loc2  = loc1
            d1min = dd
            loc1  = i
          else if(dd.lt.d2min.and.dd.ge.d1min)then
            d2min = dd
            loc2  = i
          end if
100     continue
        if(loc1.eq.0.or.loc2.eq.0)then
          return
        end if
!       loc1 = 2
!       loc2 = 4

! Additive/Multiplicative Correction
        call correction(loc1,ndim,xin,yhat,ncr1,ncr2,dad,dmu, &
                        yadd1,yadd2,ymul1,ymul2)
        if(     mLC.eq.1)then ! additive correction
          yhat = yadd2
        else if(mLC.eq.2)then ! multiplicative correction
          yhat = ymul2
        else if(mLC.eq.3)then ! to agree at 2nd nearest sample
          do k=1,ndim
            xp(k) = dad( 1+k,loc2 )
          end do  
          call correction(loc1,ndim,xp,yhat,ncr1,ncr2,dad,dmu, &
                          yadp1,yadp2,ymup1,ymup2)
          if(yadp2.eq.ymup2)then
            gam = 0.5d0
          else
            nhes = (ndim*ndim-ndim)/2 + ndim
            yhip = dad(1+ndim+1+ndim+nhes+1,loc2)
            gam  = (yhip-ymup2)/(yadp2-ymup2)
          end if
          yhat = gam*yadd2 + (1.d0-gam)*ymul2
        else if(mLC.ge.4)then
          d1 = 0.d0
          d2 = 0.d0
          dp = 0.d0
          do k=1,ndim
            xp(k) = (dad(1+k,loc2)+dad(1+k,loc1))*0.5d0
            d1 = d1 + (xin(k)-dad(1+k,loc1))**2
            d2 = d2 + (xp( k)-dad(1+k,loc1))**2
            dp = dp + (xin(k)-dad(1+k,loc1))*(xp(k)-dad(1+k,loc1))
          end do  
          d1 = sqrt(d1)
          d2 = sqrt(d2)
          if(d1.eq.0.d0)then
            gam = 0.5d0
            yhat = gam*yadd2 + (1.d0-gam)*ymul2
          else
            co = dp / d1 / d2
            if(co.eq.0.d0)then
              write(*,*)'*co in loc',(xin(k),k=1,ndim)
              co = dcos(89.9d0*4.d0*datan(1.d0)/180.d0)
            end if
            do k=1,ndim
              xp(k) = dad(1+k,loc1) + (xin(k)-dad(1+k,loc1))/d1 * (d2/co)
            end do
            dd = 0.d0
            do k=1,ndim
              dd = dd + (xp( k)-dad(1+k,loc1))**2
            end do
            dd = sqrt(dd)
            if(dd.le.0.d0)stop'|xp-L1| =< 0 in local'
            call meta_estimation(0,xp,yhap,yhatprime,RMSE,EI, &
!                               mEI,mLC,id_proc,ndim,kreg,iscf,mode_dck, &
                                mEI,  0,id_proc,ndim,kreg,iscf,mode_dck, &
                                rsample,nsize, &
                                train,mean,devi,theta,power,Rinv, &
                                RYFB,FR,FRFinv,ymin,ymin,dxadd,       &
                                nparent,inf, &        ! dummy
                                ncr1,ncr2,dad,dmu )   ! dummy
!           call dutch_interp(ydtc,loc1,loc2,ndim,ncr1,ncr2,dad,dmu ) 
            call correction(loc1,ndim,xp,yhap,ncr1,ncr2,dad,dmu, &
                            yad11,yad12,ymu11,ymu12)
            call correction(loc2,ndim,xp,yhap,ncr1,ncr2,dad,dmu, &
                            yad21,yad22,ymu21,ymu22)
            if(mLC.eq.4)then ! agreement at xp with same gam
              if(yad12.eq.ymu12)then
                gam = 0.5d0
              else
                gam = (yhap-ymu12)/(yad12-ymu12)
!               gam = (ydtc-ymu12)/(yad12-ymu12)
              end if
              yhat = gam*yadd2 + (1.d0-gam)*ymul2
            else if(mLC.eq.5)then ! agreement at xp only by additive
              if(yad12.eq.0.d0)stop'yad12 for mLC=5'
              r_add = yhap/yad12
!             r_add = ydtc/yad12
              yhat  = yadd2 * ( (r_add-1.d0)/dd/dd*d1*d1 + 1.d0 )
            else if(mLC.eq.6)then ! agreement at xp only by multiplicative
              if(ymu12.eq.0.d0)stop'ymu12 for mLC=6'
              r_mul = yhap/ymu12
!             r_mul = ydtc/ymu12
              yhat  = ymul2 * ( (r_mul-1.d0)/dd/dd*d1*d1 + 1.d0 )
            end if
          end if
!         write(100,'(99e15.5)')(xin(k),k=1,ndim),gam,(xp(k),k=1,ndim)
        else
          stop'unknown mLC in Local Correction'
        end if

        end subroutine local_correction
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine dutch_interp(yout,loc1,loc2,ndim,ncr1,ncr2,dad,dmu)
        implicit none
        integer, intent(in) :: loc1,loc2,ncr1,ncr2,ndim
        double precision, dimension(ncr1,ncr2), intent(in) :: dad,dmu
        double precision, intent(out) :: yout
        integer :: i,j,k,l,ict
        integer :: io1,io2,io,ns,nhes,loc,kd
        double precision :: f,yd,yd1,yd2,prd,ank
        double precision :: d21,dc1,d2c
        double precision, dimension(1,1) :: xg,xhx
        double precision, dimension(1,ndim) :: xh
        double precision, dimension(ndim,1) :: x,g
        double precision, dimension(ndim) :: xc,x1,x2
        double precision, dimension(ndim,ndim) :: h

        io1 = int(dad(1,loc1))
        if( int(dad(1,loc1)).eq.3 ) io1 = 2
        io2 = int(dad(1,loc2))
        if( int(dad(1,loc1)).eq.3 ) io2 = 2
        do i=1,ndim
          x1(i) = dad(1+i,loc1)
          x2(i) = dad(1+i,loc2)
          xc(i) = 0.5d0*( dad(1+i,loc1)+dad(1+i,loc2) )
        end do
        nhes = (ndim*ndim-ndim)/2 + ndim
        ns   = 1+ndim+1+ndim+nhes

        do l = 1,2
          if(l.eq.1)then
            loc = loc1
            io  = io1
          else if(l.eq.2)then
            loc = loc2
            io  = io2
          end if
          do i=1,ndim
            x(i,1) = xc(i) - dad(1+i,loc)
          end do
          f = dad(ns+1,loc)
          do i=1,ndim
            if(io.ne.0) g(i,1) = dad(ns+1+i,loc)
            ict = 0
            do j=i,ndim
              ict = ict + 1
              if(io.eq.2) h(i,j) = dad(ns+1+ndim+ict,loc)
              h(j,i) = h(i,j)
            end do
          end do

          yd = 0.d0
          do k=0,io
            kd = 1
            do i=1,k
              kd = kd * i
            end do
            ank = 1.d0 - dble(k)/dble(io+1)
            if(k.eq.0)then
              prd = f
            else if(k.eq.1)then
              xg  = matmul( transpose(x),g )
              prd = xg(1,1)
            else if(k.eq.2)then
              xh  = matmul( transpose(x),h )
              xhx = matmul( xh,x )
              prd = xhx(1,1)
            end if
            yd  = yd + ank/dble(kd) * prd
          end do

          if(l.eq.1)then
            yd1 = yd
          else if(l.eq.2)then
            yd2 = yd
          end if
        end do ! l loop

        d21 = 0.d0
        d2c = 0.d0
        dc1 = 0.d0
        do i=1,ndim
           d21 = d21 + (x2(i)-x1(i))**2
           d2c = d2c + (x2(i)-xc(i))**2
           dc1 = dc1 + (xc(i)-x1(i))**2
        end do
        d21 = sqrt(d21)
        d2c = sqrt(d2c)
        dc1 = sqrt(dc1)
        yout = d2c/d21*yd1 + dc1/d21*yd2

        end subroutine dutch_interp
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine correction(nd,ndim,xin,yhat,ncr1,ncr2,dad,dmu, &
                              yad1,yad2,ymu1,ymu2)
        implicit none
        integer, intent(in) :: nd,ndim,ncr1,ncr2
        double precision, intent(in) :: yhat
        double precision, dimension(ndim), intent(in) :: xin
        double precision, dimension(ncr1,ncr2), intent(in) :: dad,dmu
        double precision, intent(out) :: yad1,yad2,ymu1,ymu2
        integer :: j,k
        integer :: mode,ict
        double precision :: c1,c2
        double precision, dimension(1,1) :: cor1,cor2
        double precision, dimension(ndim,1) :: x,dc1,dc2,ct
        double precision, dimension(ndim,ndim) :: d2c1,d2c2

        mode = int(dad(1,nd))
        c1   = dad(1+ndim+1,nd)
        c2   = dmu(1+ndim+1,nd)
        do j=1,ndim
          x(j,1) = xin(j) - dad(1+j,nd)
        end do
        if(mode.eq.0)then
           dc1  = 0.d0
           d2c1 = 0.d0
           dc2  = 0.d0
           d2c2 = 0.d0
        else if(mode.eq.1)then
           do j=1,ndim
             dc1(j,1) = dad(1+ndim+1+j,nd)
             dc2(j,1) = dmu(1+ndim+1+j,nd)
           end do
           d2c1 = 0.d0
           d2c2 = 0.d0
        else if(mode.eq.3)then
           ict = 0
           do j=1,ndim
             dc1(j,1) = dad(1+ndim+1+j,nd)
             dc2(j,1) = dmu(1+ndim+1+j,nd)
             do k=j,ndim
               ict = ict + 1
               d2c1(j,k) = dad(1+ndim+1+ndim+ict,nd)
               d2c2(j,k) = dmu(1+ndim+1+ndim+ict,nd)
               d2c1(k,j) = d2c1(j,k)
               d2c2(k,j) = d2c2(j,k)
             end do
           end do
        else
          stop'unknown mode of nd in correction'
        end if

        cor1 = matmul(transpose(dc1),x)  ! (1*n)x(n*1)
        ct   = matmul(d2c1,x)            ! (n*n)x(n*1)
        cor2 = matmul(transpose(x),ct)   ! (1*n)x(n*1)
!       write(*,'(a,i5,4f12.7)')'*add  ',nd,yhat,c1,cor1(1,1),cor2(1,1)
        yad1 = yhat + (c1 + cor1(1,1))
        yad2 = yhat + (c1 + cor1(1,1) + 0.5d0*cor2(1,1))

        cor1 = matmul(transpose(dc2),x)  ! (1*n)x(n*1)
        ct   = matmul(d2c2,x)            ! (n*n)x(n*1)
        cor2 = matmul(transpose(x),ct)   ! (1*n)x(n*1)
!       write(*,'(a,i5,4f12.7)')'*mul  ',nd,yhat,c2,cor1(1,1),cor2(1,1)
        ymu1 = yhat * (c2 + cor1(1,1))
        ymu2 = yhat * (c2 + cor1(1,1) + 0.5d0*cor2(1,1))

        end subroutine correction
