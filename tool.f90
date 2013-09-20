        subroutine allocate_krig(igrad)
        use dimKrig
        implicit none
        integer, intent(in) :: igrad
        integer :: iflag,istat

        if(igrad.eq.0.or.igrad.eq.2)then ! Kriging
          nsize = rsample
        else if(igrad.eq.1)then          ! Direct
          if(mode_dck.eq.1)then
            nsize = ict_sample(0,1)*(          1) &
                  + ict_sample(0,2)*(     ndim+1) &
                  + ict_sample(0,3)*(ndim     +1) &
                  + ict_sample(0,4)*(ndim+ndim+1) &
                  + ict_sample(0,5)*(ndim     +1) &
                  + ict_sample(0,6)*(ndim+ndim+1)
          else
            nsize = ict_sample(0,1)*(          1) &
                  + ict_sample(0,2)*(     ndim+1) &
                  + ict_sample(0,3)*(nhes     +1) &
                  + ict_sample(0,4)*(nhes+ndim+1) &
                  + ict_sample(0,5)*(ndim     +1) &
                  + ict_sample(0,6)*(ndim+ndim+1)
          end if
        else if(igrad.eq.-1)then         ! already given before

        else
          stop'unknown igrad in allocate_krig'
        end if

        tdim  = lmax*(lmax+1)/2
        iflag = 0
        allocate( Rij(nsize,nsize),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( Rinv(nsize,nsize),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( yy(nsize),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( r(nsize,1),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( d_theta(ndim,tdim),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( d_power(ndim,tdim),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( RYFB(nsize,1),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( iRij(rsample,4),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1

        if(igrad.eq.-1)then    ! already given
        else if(lmax.eq.1)then ! Single-Fidelity
          kreg = kreg_orig
        else                   ! Variable-Fidelity
          if(kreg_orig.ne.1)stop'not yet for unversal Krig with VF'
          kreg = kreg_orig * lmax
        end if
          
        allocate( mean(kreg,1),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( FB(nsize,kreg),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( FR(kreg,nsize),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( FRFinv(kreg,kreg),stat=istat )
                            if(istat.ne.0)iflag = iflag + 1
        allocate( deviratio(lmax),stat=istat)
                            if(istat.ne.0)iflag = iflag + 1
        allocate( deviratio_best(lmax),stat=istat)
                            if(istat.ne.0)iflag = iflag + 1
      
        if(iflag.ne.0)stop'*Allocation Miss in Krig'
        
        end subroutine allocate_krig
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine deallocate_krig(igrad)
        use dimKrig
        implicit none
        integer, intent(in) :: igrad

        deallocate(Rij,Rinv,r)
        deallocate(d_theta,d_power)
        deallocate(FB,RYFB,FR,FRFinv)
        deallocate(yy,mean)
        deallocate(deviratio,deviratio_best)
        deallocate(iRij)

        end subroutine deallocate_krig
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine make_training_data(igrad)
        use dimKrig
        implicit none
        integer, intent(in) :: igrad
        integer :: i,j,k,l,ict,hct
        double precision, dimension(ndim)      :: x
        double precision, dimension(kreg_orig) :: y
        character(len=2) :: cl

! yy [func, gfunc, hfunc, hvec]_fid1, []_fid2, ....
        ict = 0
        hct = 0
        do 100 l=1,lmax
         write(cl,101)l
101      format(i1,'_')
         do 110 i=1,rsample
            if(inf(i)(1:2).ne.cl)go to 110
            ict = ict + 1
            yy(ict) = fun(i)
110      continue
         if(igrad.ne.1) go to 199 ! for Direct

         do 120 i=1,rsample
            if(inf(i)(1:2).ne.cl)go to 120
            if(inf(i)(3:6).eq.'FG  '.or.inf(i)(3:6).eq.'FGH '.or. &
               inf(i)(3:6).eq.'FGHv')then
               do j=1,ndim
                 ict = ict + 1
                 yy(ict) = gfun(i,j)
               end do
            end if
120      continue
         do 130 i=1,rsample
            if(inf(i)(1:2).ne.cl)go to 130
            if(inf(i)(3:6).eq.'FH  '.or.inf(i)(3:6).eq.'FGH ')then
              if(mode_dck.eq.1)then
                do j=1,ndim
                 ict = ict + 1
                 yy(ict) = hfun(i,j,j)
                end do
              else
                do j=1,ndim
                do k=j,ndim
                 ict = ict + 1
                 yy(ict) = hfun(i,j,k)
                end do
                end do
              end if
            end if
130      continue
         do 140 i=1,rsample
            if(inf(i)(1:2).ne.cl)go to 140
            if(inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv')then
               do j=1,ndim
                 ict = ict + 1
                 yy(ict) = hvec(i,j,2)
               end do
            end if
140      continue

199      continue
         if(mode_dck.eq.1)then
            hct = hct + ict_sample(l,1)*(1) &
                      + ict_sample(l,2)*(1+ndim) &
                      + ict_sample(l,4)*(1+ndim+ndim) &
                      + ict_sample(l,6)*(1+ndim+ndim)
         else
            hct = hct + ict_sample(l,1)*(1) &
                      + ict_sample(l,2)*(1+ndim) &
                      + ict_sample(l,4)*(1+ndim+nhes) &
                      + ict_sample(l,6)*(1+ndim+ndim)
         end if
         if(ict.ne.hct)stop'ict.ne.hct in make_training_data'
100     continue ! lmax loop

        if(ict.ne.nsize)stop'ict!=nsize in make_training_data'
        if(id_proc.eq.0) &
        write(filenum,'(1x,a,i6)')'>> Size of Correlation Matrix = ',nsize

! FB regression
        FB  = 0.d0 ! (nsize,kreg)
        ict = 0
        hct = 0
        do 200 l=1,lmax
          write(cl,101)l
          do 210 i=1,rsample
            if(inf(i)(1:2).ne.cl)go to 210
            x(:) = sampl(i,:)
            call get_regression(0,0,0,x,y) ! yout(kreg_orig)
            ict = ict + 1
            do j=1,kreg_orig
              FB(ict,kreg_orig*(l-1)+j)  = y(j)
            end do
210       continue
          if(mode_dck.eq.1)then
            hct = hct + ict_sample(l,1)*(1) &
                      + ict_sample(l,2)*(1+ndim) &
                      + ict_sample(l,4)*(1+ndim+ndim) &
                      + ict_sample(l,6)*(1+ndim+ndim)
          else
            hct = hct + ict_sample(l,1)*(1) &
                      + ict_sample(l,2)*(1+ndim) &
                      + ict_sample(l,4)*(1+ndim+nhes) &
                      + ict_sample(l,6)*(1+ndim+ndim)
          end if
          ict = hct
200     continue ! lmax loop
        return
 
        do l=1,lmax
          write(*,'(998f8.2)')(FB(i,l),i=1,nsize)
        end do
        call stop_all

        end subroutine make_training_data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine likelihood(theta,power,llfd,verr)
        use dimKrig
        implicit none
        double precision, dimension(ndim,tdim), intent(in) :: &
        theta,power
        double precision, intent(out) :: llfd,verr
        integer :: k,l,iflag
        double precision :: det,detlog

        call makeMatrixRij(theta,power)
        call LUDecomposition(1,0,nsize,Rij,Rinv,det,detlog,verr)
        call makeMeanDevi

        iflag = 0
        if(devi.le.0.d0) iflag = 1
        do l=2,lmax
        if(deviratio(l).le.0.d0)iflag = 1
        end do

        if(iflag.eq.1)then
          if(id_proc.eq.0)then
            write(filenum,'(a,99e12.3)') &
            '*negative devi = ',devi,(deviratio(l),l=2,lmax)
            if(ndebug.eq.1)then
             write(filenum,'(a,6i5    )')'*',ndim,nsample,rsample,nsize,lmax
             do l=1,lmax
             write(filenum,'(a,100f8.3)')'*theta',(theta(k,l),k=1,ndim)
             end do
             write(filenum,'(a,100f8.3)')'*mean ',(mean(k,1),k=1,kreg)
             write(filenum,'(a,100f8.3)')'*Rinv ',det,detlog,verr
            end if
          end if
          llfd = -1.d10
          verr =  10.d0
        else
          llfd = -0.5d0*(nsize)*dlog(devi) - 0.5d0*detlog
          if(llfd.ge.llfd_best.and.verr.lt.1.d-3)then
            llfd_best      = llfd      ! update in each processor
            deviratio_best = deviratio ! and then gather at evaluation
          end if
        end if

        end subroutine likelihood
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine update_deviratio    ! update for all processors
        use dimKrig
        implicit none
        include 'mpif.h'
        integer :: l,ierr,id,idg
        double precision :: lbglb

        call MPI_ALLREDUCE(llfd_best,lbglb,1,MPI_DOUBLE_PRECISION, &
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        id = -1
        if(lbglb.eq.llfd_best) id = id_proc
        call MPI_ALLREDUCE(id,idg,1,MPI_INTEGER, &
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        if(idg.lt.0.or.idg.ge.num_proc)stop'idg for deviratio'

        call MPI_BCAST(deviratio_best,lmax,MPI_DOUBLE_PRECISION,idg, &
                       MPI_COMM_WORLD,ierr)
        llfd_best = lbglb

!       if(id_proc.eq.0) &
!       write(*,'(11x,a,99e15.5)') &
!       '>> Update Deviratio, ',llfd_best,(deviratio_best(l),l=1,lmax)

        end subroutine update_deviratio
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine makeMeanDevi
        use dimKrig 
        implicit none
! dimension mean(kreg,1)
!           devi
!           deviratio(lmax)
        double precision, dimension(nsize,1)    :: YD,YFB,Yp,Ym
        double precision, dimension(kreg,kreg)  :: mpar,mparinv
        double precision, dimension(kreg,1)     :: mchr
        double precision, dimension(kreg,nsize) :: TMP1
        double precision, dimension(1,nsize)    :: TMP2
        double precision, dimension(1,1)        :: pRp,pRfbm,chld
! tmp
        integer :: i,j,l
        double precision :: tdet,tdetlog,terr

! mean
        call makeYD(nsize,lmax,ndim,mode_dck,deviratio_best,ict_sample,yy,YD)
                                             ! make with previous
        TMP1 = matmul( transpose(FB),Rinv )  ! (k*ns)x(ns*ns)
        mpar = matmul( TMP1,FB )             ! (k*ns)x(ns*k)
        mchr = matmul( TMP1,YD )             ! (k*ns)x(ns*1)
        if(kreg.eq.1)then
          if(mpar(1,1).eq.0.d0)stop'mpar=0 in makeMean'
          mean(1,1) = mchr(1,1)/mpar(1,1)
        else
          call LUDecomposition(0,0,kreg,mpar,mparinv,tdet,tdetlog,terr)
          mean = matmul(mparinv,mchr)        ! (k*k)x(k*1)
        end if

! deviratio
        deviratio(1) = 1.d0 ! dummy
        if(lmax.ne.1)then
          if(kreg.ne.lmax)stop'kreg.ne.lmax in makeMeanDevi for VFM'
          YFB = matmul(FB,mean)               ! (ns*k)x(k*1)
          do l=2,lmax
            call makeYpm(l,nsize,ndim,mode_dck,ict_sample,yy,YD,Yp,Ym)
            YFB(:,1) = YFB(:,1) - Ym(:,1)
            TMP2  = matmul(transpose(Yp),Rinv)! (1*ns)x(ns*ns)
            pRp   = matmul(TMP2,Yp)           ! (1*ns)x(ns*1)
            pRfbm = matmul(TMP2,YFB)          ! (1*ns)x(ns*1)
            if(pRp(1,1).eq.0.d0)stop'pRp=0 in deviratio'
            deviratio(l) = pRfbm(1,1) / pRp(1,1)

            if(deviratio(l).lt.devmin)deviratio(l) = devmin
            if(dabs(deviratio(l)).le.1.d-10)deviratio(l) = 1.d-10
            if(deviratio(l).gt.devmax)deviratio(l) = devmax
!           deviratio(l) = 1.0d0
          end do
          call makeYD(nsize,lmax,ndim,mode_dck,deviratio,ict_sample,yy,YD)
                                              ! make with current
        end if
! devi
        Yp       = matmul(FB,mean)             ! (ns*k)x(k*1), tmp
        YFB(:,1) = YD(:,1) - Yp(:,1)           ! (ns*1), Y-FB
        TMP2     = matmul(transpose(YFB),Rinv) ! (1*ns)x(ns*ns)
        chld     = matmul(TMP2,YFB)            ! (1*ns)x(ns*1)
        devi     = chld(1,1)/dble(nsize)

        end subroutine makeMeanDevi
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine makeYD(nsize,lmax,ndim,mode_dck,deviratio,ict,yy,YD)
        implicit none
        integer, intent(in) :: nsize,lmax,ndim,mode_dck
        integer, dimension(0:10,0:10), intent(in) :: ict
        double precision, dimension(lmax), intent(in) :: deviratio
        double precision, dimension(nsize), intent(in) :: yy
        double precision, dimension(nsize,1), intent(out) :: YD
        integer :: i,l,is,ie,nh

        nh = ndim*(ndim+1)/2
        if(mode_dck.eq.1)nh = ndim

        is = 1
        ie = 0
        do 100 l=1,lmax
           ie = ie + ict(l,0) + ict(l,2)*(ndim) &
                              + ict(l,3)*(     nh) &
                              + ict(l,4)*(ndim+nh) &
                              + ict(l,5)*(        ndim) &
                              + ict(l,6)*(ndim   +ndim)
           do 110 i=is,ie
             if(l.eq.1)then
               YD(i,1) = yy(i)
             else
               YD(i,1) = yy(i) * deviratio(l)
             end if
110        continue
!          write(*,'(a,4i5)')'*YD',l,is,ie,nsize
           is = ie + 1
100     continue
        if(ie.ne.nsize)stop'ie.ne.nsize in makeYD'

        end subroutine makeYD
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine makeYpm(lin,nsize,ndim,mode_dck,ict,yy,YD,Yp,Ym)
        implicit none
        integer, intent(in) :: lin,nsize,ndim,mode_dck
        integer, dimension(0:10,0:10), intent(in) :: ict
        double precision, dimension(nsize),   intent(in)  :: yy
        double precision, dimension(nsize,1), intent(in)  :: YD
        double precision, dimension(nsize,1), intent(out) :: Yp,Ym
        integer :: i,l,is,ie,nh

        nh = ndim*(ndim+1)/2
        if(mode_dck.eq.1)nh = ndim

        Yp(:,1) = 0.d0
        Ym(:,1) = YD(:,1)
        is = 1
        ie = 0
        do 100 l=1,lin
           ie = ie + ict(l,0) + ict(l,2)*(ndim) &
                              + ict(l,3)*(     nh) &
                              + ict(l,4)*(ndim+nh) &
                              + ict(l,5)*(        ndim) &
                              + ict(l,6)*(ndim   +ndim)
           if(l.ne.lin)is = ie + 1
100     continue

        do 200 i=is,ie
           Yp(i,1) = yy(i)
           Ym(i,1) = 0.d0
200     continue

        return
        do i=1,nsize
          write(*,'(2i5,4f12.7)')lin,i,yy(i),YD(i,1),Yp(i,1),Ym(i,1)
        end do
        call stop_all
        end subroutine makeYpm
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine get_regression(mode,idif,jdif,xin,yout)
        use dimKrig
        implicit none
        integer, intent(in) :: mode,idif,jdif
        double precision, dimension(ndim), intent(in)  :: xin
        double precision, dimension(kreg_orig), intent(out) :: yout
        integer :: i,j,ip

        yout = 1.d0
        if(mode.eq.0)then
          do i=1,kreg_orig
             do j=1,ndim
                ip = mreg(i,j)
                if(ip.gt.abs(kord))stop'ip>korder'
                yout(i) = yout(i) * ((xin(j))**(ip))
             end do
          end do
        else if(mode.eq.1)then
          do i=1,kreg_orig
             do j=1,ndim
                ip = mreg(i,j)
                if(ip.gt.abs(kord))stop'ip>korder'
                if(idif.eq.j)then
                  if(ip.eq.0)then
                    yout(i) = yout(i) * 0.d0
                  else
                    yout(i) = yout(i) * (dble(ip)*(xin(j))**(ip-1))
                  end if
                else
                  yout(i) = yout(i) * ((xin(j))**(ip))
                end if
             end do
          end do
        else if(mode.eq.2)then
          do i=1,kreg_orig
             do j=1,ndim
                ip = mreg(i,j)
                if(ip.gt.abs(kord))stop'ip>korder'
                if(      idif.eq.j.and.jdif.eq.j )then
                  if(ip.le.1)then
                    yout(i) = yout(i) * 0.d0
                  else
                    yout(i) = yout(i) * (dble(ip*(ip-1))*(xin(j))**(ip-2))
                  end if
                else if((idif.NE.j.and.jdif.eq.j).or. &
                        (idif.eq.j.and.jdif.NE.j))then
                  if(ip.eq.0)then
                    yout(i) = yout(i) * 0.d0
                  else
                    yout(i) = yout(i) * (dble(ip)*(xin(j))**(ip-1))
                  end if
                else if( idif.NE.j.and.jdif.NE.j )then
                  yout(i) = yout(i) * ((xin(j))**(ip))
                end if
             end do
          end do
        end if

        end subroutine get_regression
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine Make_Mreg
        use dimKrig
        implicit none
        integer :: i,j,k,isum,ict

        mreg = 0    ! mreg(kreg_orig,ndim)
        if(kord.eq.0)then

        else if(kord.lt.0)then ! diagonal basis
         ict = 1
         do 300 k=1,abs(kord)
         do 300 i=1,ndim
            ict = ict + 1
            mreg(ict,i) = k
300      continue
         if(ict.ne.kreg_orig)stop'kreg_orig for diagonal_basis'

        else                   ! general basis
         do 100 i=1,kreg_orig+5
           isum = 0
           do j=1,ndim
             isum = isum + mreg(i,j)
           end do
           if(isum.ne.kord)then
             mreg(i+1,:) = mreg(i,:)
             mreg(i+1,1) = mreg(i+1,1) + 1
             go to 100
           else
             do j=1,ndim
               if(mreg(i,j).ne.0)then
                 if(j.eq.ndim)go to 200
                 mreg(i+1,:  ) = mreg(i,:)
                 mreg(i+1,j  ) = 0
                 mreg(i+1,j+1) = mreg(i+1,j+1) + 1
                 go to 100
               end if
             end do
             stop'No target j in Make_Mreg'
           end if
100      continue
         write(*,'(3i6)')kord,kreg_orig,ndim
         stop'Error in Make_Mreg'
200      continue
         if(i.ne.kreg_orig)then
           write(*,*)i,kreg_orig,kord
           stop'i.ne.kreg_orig in Make_Mreg'
         end if
        end if

        return
        write(*,'(a,3i6)')'*ndim,kord,kreg_orig = ',ndim,kord,kreg_orig
        do i=1,kreg_orig
         write(*,'(i5,a,99i4)')i,':',(mreg(i,j),j=1,ndim)
        end do
        call stop_all
        end subroutine Make_Mreg
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine read_krig(mode,ifunc)
        use dimKrig
        implicit none
! mode=0 : usual
! mode=1 : only to get theta and power
! mode=2 : check the size of matrix
        integer, intent(in) :: mode,ifunc
        character(len=10) :: Cfile
        integer :: mdim,msample,qsample,msize,mmeg,mscf,mmeg_orig
        integer :: mmode_dck,mlmax,mtdim,idum
        integer :: i,j,k,l
        double precision :: dum

        if(ifunc.le.9)then
          write(Cfile,11)ifunc
        else if(ifunc.ge.10.and.ifunc.le.99)then
          write(Cfile,12)ifunc
        else
          stop'ifunc.gt.100'
        end if
11      format('meta0',i1,'.krg')
12      format('meta',i2,'.krg')


        if(mode.eq.1)then
          open(10,file=Cfile,form='unformatted',status='old')
          read(10) mdim,msample,qsample,msize
          read(10) mmeg_orig,mmeg,mscf,mmode_dck
          read(10) mlmax,mtdim
          read(10) ((idum,i=0,10),j=0,10)
          if(ndim.ne.mdim)then
            write(*,'(2a,2i5)')'*ndim is different in ',Cfile,ndim,mdim
            call stop_all
          else if(mscf.ne.iscf)then
            write(*,'(2a,2i5)')'*iscf is different in ',Cfile,iscf,mscf
            call stop_all
          else if(mlmax.ne.lmax.or.mtdim.ne.tdim)then
            write(*,'(2a,2i5)')'*lmax is different in ',Cfile,lmax,mlmax
            call stop_all
          end if
          read(10)((dum,k=1,ndim),i=1,qsample) !sampl
          read(10) (dum,k=1,mmeg)              !(mean(k,1),k=1,kreg)
          if(mlmax.eq.1)then
          read(10)  dum                        !devi
          else
          read(10)  dum,(dum,l=2,lmax)        !devi,ratio
          end if
          read(10) ((d_theta(k,l),k=1,ndim),l=1,tdim)
          if(iscf.eq.0.or.iscf.eq.1) &
          read(10) ((d_power(k,l),k=1,ndim),l=1,tdim)
          close(10)
          return
        else if(mode.eq.2)then
          open(10,file=Cfile,form='unformatted',status='old')
          read(10) mdim,msample,qsample,msize
          read(10) mmeg_orig,mmeg,mscf,mmode_dck
          read(10) mlmax,mtdim
          close(10)
          if(ndim.ne.mdim)then
            write(*,'(2a,2i5)')'*ndim is different in ',Cfile,ndim,mdim
            call stop_all
          else if(msample.ne.nsample)then
            write(*,'(2a,2i5)')'*nsample is different in ',Cfile,nsample,msample
            call stop_all
          else if(mscf.ne.iscf)then
            write(*,'(2a,2i5)')'*iscf is different in ',Cfile,iscf,mscf
            call stop_all
          else if(mmeg_orig.ne.kreg_orig)then
            write(*,'(2a,2i5)')'*kreg_orig is different in ', &
            Cfile,kreg_orig,mmeg_orig
            call stop_all
          end if
          rsample = qsample
          nsize   = msize
          kreg    = mmeg
          lmax    = mlmax
          tdim    = mtdim
          return
        end if

        ! conventional reading
        open(10,file=Cfile,form='unformatted',status='unknown')
        read(10) mdim,msample,qsample,msize
        read(10) mmeg_orig,mmeg,mscf,mmode_dck
        read(10) mlmax,mtdim
        read(10) ((ict_sample(i,j),j=0,10),i=0,10)

        if(ndim.ne.mdim)then
          write(*,'(2a,2i5)')'*ndim is different in ',Cfile,ndim,mdim
          call stop_all
        else if(msample.ne.nsample)then
          write(*,'(2a,2i5)')'*nsample is different in ',Cfile,nsample,msample
          call stop_all
        else if(qsample.ne.rsample)then
          write(*,'(2a,2i5)')'*rsample is different in ',Cfile,rsample,qsample
          call stop_all
        else if(msize.ne.nsize)then
          write(*,'(2a,2i5)')'*nsize is different in ',Cfile,nsize,msize
          call stop_all
        else if(mmeg_orig.ne.kreg_orig)then
          write(*,'(2a,2i5)')'*kreg_orig is different in ', &
                             Cfile,kreg_orig,mmeg_orig
          call stop_all
        else if(mscf.ne.iscf)then
          write(*,'(2a,2i5)')'*iscf is different in ',Cfile,iscf,mscf
          call stop_all
        else if(mode_dck.ne.mmode_dck)then
          write(*,'(2a,2i5)')'*mode_dck is different in ', &
                             Cfile,mode_dck,mmode_dck
          call stop_all
        else if(lmax.ne.mlmax.or.tdim.ne.mtdim)then
          write(*,'(2a,4i5)')'*lmax is different in ', &
                             Cfile,lmax,mlmax,tdim,mtdim
          call stop_all
        end if

        kreg = mmeg
        if(kreg.ne.kreg_orig.and.kreg.ne.kreg_orig*mlmax) &
        stop'kreg? in read_krig'

        read(10)((sampl(i,k),k=1,ndim),i=1,rsample)
        read(10) (mean(k,1),k=1,kreg)
        if(lmax.eq.1)then
          read(10)  devi
        else
          read(10)  devi,(deviratio(l),l=2,lmax)
        end if
        read(10) ((d_theta(k,l),k=1,ndim),l=1,tdim)
        if(iscf.eq.0.or.iscf.eq.1) &
        read(10) ((d_power(k,l),k=1,ndim),l=1,tdim)

        read(10)((Rinv(i,j),i=1,nsize),j=1,nsize)
        read(10) (yy(i),i=1,nsize)
        read(10)((FB(i,k),k=1,kreg),i=1,nsize)
        read(10) (RYFB(i,1),i=1,nsize)
        read(10)((FR(i,j),j=1,nsize),i=1,kreg)
        read(10)((FRFinv(i,j),j=1,kreg),i=1,kreg)
        read(10) (nparent(i),i=1,rsample)
        read(10) (dxadd(i),i=1,rsample)
        read(10) (inf(i),i=1,rsample)
        read(10)(((hvec(i,j,k),k=1,2),j=1,ndim),i=1,rsample)
        close(10)

        end subroutine read_krig
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine write_krig(ifunc,theta,power)
        use dimKrig
        implicit none
        integer, intent(in) :: ifunc
        double precision, dimension(ndim,tdim), intent(in) :: theta,power

        integer :: i,j,k,l
        double precision :: tdet,tdetlog,terr
        character(len=10) :: Cfile
        double precision, dimension(nsize,1)   :: YD,YFB,FBm
        double precision, dimension(kreg,kreg) :: FRF

        call makeYD(nsize,lmax,ndim,mode_dck,deviratio,ict_sample,yy,YD)
        FBm      = matmul(FB,mean)             ! (ns*k)x(k*1)
        YFB(:,1) = YD(:,1) - FBm(:,1)          ! (ns*1), Y-FB
        RYFB     = matmul(Rinv,YFB)            ! (ns*ns)x(ns*1)
        FR       = matmul(transpose(FB),Rinv)  ! (k*ns)x(ns*ns)
        FRF      = matmul(FR,FB)               ! (k*ns)x(ns*k)

        if(kreg.eq.1)then
          if(FRF(1,1).eq.0.d0)stop'FRF=0 in write_krig'
          FRFinv(1,1) = 1.d0/FRF(1,1)
        else
          call LUDecomposition(0,0,kreg,FRF,FRFinv,tdet,tdetlog,terr)
        end if

        if(ifunc.le.9)then
          write(Cfile,11)ifunc
        else if(ifunc.ge.10.and.ifunc.le.99)then
          write(Cfile,12)ifunc
        else
          stop'ifunc.gt.100'
        end if
11      format('meta0',i1,'.krg')
12      format('meta',i2,'.krg')

        if(id_proc.eq.0)then
          write(filenum,'(1x,2a)')'>> Output to ',Cfile
          open( 10,file=Cfile,form='unformatted',status='unknown')
          write(10) ndim,nsample,rsample,nsize
          write(10) kreg_orig,kreg,iscf,mode_dck
          write(10) lmax,tdim
          write(10)((ict_sample(i,j),j=0,10),i=0,10)

          write(10)((sampl(i,k),k=1,ndim),i=1,rsample)
          write(10) (mean(k,1),k=1,kreg)
          if(lmax.eq.1)then
            write(10) devi
          else
            write(10) devi,(deviratio(l),l=2,lmax)
          end if
          write(10)((theta(k,l),k=1,ndim),l=1,tdim)
          if(iscf.eq.0.or.iscf.eq.1) &
          write(10)((power(k,l),k=1,ndim),l=1,tdim)

          write(10)((Rinv(i,j),i=1,nsize),j=1,nsize)
          write(10) (yy(i),i=1,nsize)
          write(10)((FB(i,k),k=1,kreg),i=1,nsize)
          write(10) (RYFB(i,1),i=1,nsize)
          write(10)((FR(i,j),j=1,nsize),i=1,kreg)
          write(10)((FRFinv(i,j),j=1,kreg),i=1,kreg)
          write(10) (nparent(i),i=1,rsample)
          write(10) (dxadd(i),i=1,rsample)
          write(10) (inf(i),i=1,rsample)
          write(10) (((hvec(i,j,k),k=1,2),j=1,ndim),i=1,rsample)
          close(10)
        end if

        end subroutine write_krig
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine check_krig(ifunc,theta,power)
        use dimKrig
        implicit none
        integer, intent(in) :: ifunc
        double precision, dimension(ndim,tdim), intent(in) :: theta,power
        double precision, dimension(ndim) :: xin,dy,yhatprime
        double precision, dimension(ndim,ndim) :: d2y
        double precision, dimension(rsample,ndim) :: train
        integer :: i,j,k,itarg
        integer :: ir,ia,il,ig,ih,id,iv
        double precision :: yexa,ymin,yhat,RMSE,EI
        double precision :: verr,verr_r,verr_a,verr_l
        double precision :: vg,vgmax,vgmin,vh,vhmin,vhmax,diff
        double precision :: vd,vdmax,vdmin,vv,vvmin,vvmax
        double precision, dimension(ndim,1) :: V,HV

        d_theta = theta
        d_power = power
        do i=1,rsample
           train(i,:) = sampl(i,:)
        end do
        ncr1 = 1+ndim+1+ndim+nhes
        ncr2 = nsample
        allocate(dadd(ncr2*ncr1,nfunc-1)) ! dummy
        allocate(dmul(ncr2*ncr1,nfunc-1)) ! dummy

        itarg  = 0
        do i=1,nOPT
           if(ifunc.eq.nfOPT(i))itarg = i
        end do
        if(itarg.ne.0)then
          ymin = OFopt(itarg)
        else
          ymin = fun(rsample) ! dummy
        end if

        verr   = 0.d0
        verr_r = 0.d0
        verr_a = 0.d0
        verr_l = 0.d0
        ir = 0
        ia = 0
        il = 0
        do 100 i=1,rsample
           xin(:) = sampl(i,:)
           call meta_estimation(0,xin,yhat,yhatprime,RMSE,EI, &
                                mEI,0,id_proc,ndim,kreg_orig, &
                                iscf,mode_dck, &
                                ccrf, &
                                rsample,nsize,kreg,lmax,tdim, &
                                train,mean,devi,deviratio, &
                                d_theta,d_power,Rinv, &
                                RYFB,FR,FRFinv,ymin,ymin,dxadd,hvec,  &
                                ict_sample,nparent, &
                                inf, & 
                                ncr1,ncr2,dadd,dmul) ! dummy
           yexa   = fun(i)
           verr   = verr + (yexa-yhat)**2
           if(nparent(i).eq.0)then ! real pt
             if(inf(i)(1:2).ne.'1_')then
               il = il + 1
               verr_l = verr_l + (yexa-yhat)**2
             else
               ir = ir + 1
               verr_r = verr_r + (yexa-yhat)**2
             end if
           else
             ia = ia + 1
             verr_a = verr_a + (yexa-yhat)**2
           end if
!          if(id_proc.eq.0) &
!          write(*,'(i5,99e15.5)')i,(xin(k),k=1,ndim),yexa,yhat
!          write(*,'(i5,2e15.5)')i,yexa,yhat
100     continue
        verr   = dsqrt(verr)   / dble(rsample)
        verr_r = dsqrt(verr_r) / dble(ir)
        verr_a = dsqrt(verr_a) / dble(ia)
        verr_l = dsqrt(verr_l) / dble(il)
        if(id_proc.eq.0)then
        write(filenum,'(6x,a,e15.5,a,i6,a)') &
        '>> RMSE on All Sample Points    =',verr,' among',rsample,' pts'
        if(ir.ne.0.and.ia.ne.0)then
        write(filenum,'(6x,a,2e15.5)') &
        '>> RMSE on Real / Add Points    =',verr_r,verr_a
        end if
        if(il.ne.0.and.il.ne.rsample)then
        write(filenum,'(6x,a,e15.5,a,i6,a)') &
        '>> RMSE on High Fid.  Points    =',verr_r,' among',ir,' pts'
        write(filenum,'(6x,a,e15.5,a,i6,a)') &
        '>> RMSE on Low  Fid.  Points    =',verr_l,' among',il,' pts'
        end if
        end if

        ig = 0
        id = 0
        ih = 0
        iv = 0
        vg    = 0.d0
        vgmax = 0.d0
        vgmin = 1.d9
        vd    = 0.d0
        vdmax = 0.d0
        vdmin = 1.d9
        vh    = 0.d0
        vhmax = 0.d0
        vhmin = 1.d9
        vv    = 0.d0
        vvmax = 0.d0
        vvmin = 1.d9
        do 200 i=1,rsample
           if(rsample.eq.nsize)then
             if(nparent(i).ne.0)go to 200
           else
             if(inf(i)(1:2).ne.'1_')go to 200
             if(inf(i).eq. '1_F   ')go to 200
           end if

           xin(:) = sampl(i,:)
           call diff_estimation(xin,dy,d2y, &
                                id_proc,ndim,kreg_orig, &
                                iscf,mode_dck, &
                                ccrf, &
                                rsample,nsize,kreg,lmax,tdim, &
                                train,mean,devi,deviratio, &
                                d_theta,d_power, &
                                RYFB,dxadd,hvec, &
                                ict_sample,nparent, &
                                inf)
           if(inf(i)(3:6).eq.'FG  '.or.inf(i)(3:6).eq.'FGH '.or. &
              inf(i)(3:6).eq.'FGHv')then
             do 210 j=1,ndim
               ig = ig + 1
               diff =  abs(gfun(i,j)-dy(j))
               vg = vg + diff**2
               if(diff.gt.vgmax)vgmax = diff
               if(diff.lt.vgmin)vgmin = diff
!              write(*,'(a,2i5,2e15.5)')'*G',i,j,gfun(i,j),dy(j)
210          continue
           end if
           if(inf(i)(3:6).eq.'FH  '.or.inf(i)(3:6).eq.'FGH ')then
             do 220 j=1,ndim
             do 220 k=j,ndim
               ih = ih + 1
               diff =  abs(hfun(i,j,k)-d2y(j,k))
               vh = vh + diff**2
               if(diff.gt.vhmax)vhmax = diff
               if(diff.lt.vhmin)vhmin = diff

               if(j.eq.k)then
                 id = id + 1
                 diff =  abs(hfun(i,j,k)-d2y(j,k))
                 vd = vd + diff**2
                 if(diff.gt.vdmax)vdmax = diff
                 if(diff.lt.vdmin)vdmin = diff
               end if
!              write(*,'(a,3i5,2e15.5)')'*H',i,j,k,hfun(i,j,k),d2y(j,k)
220          continue
           end if
           if(inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv')then
             V(:,1) = hvec(i,:,1)
             HV = matmul(d2y,V)
             do 230 j=1,ndim
               iv = iv + 1
               diff = abs(hvec(i,j,2)-HV(j,1))
               vv = vv + diff**2
               if(diff.gt.vvmax)vvmax = diff
               if(diff.lt.vvmin)vvmin = diff
230          continue
           end if
200     continue
        if(id_proc.eq.0)then
          if(ig.ne.0) then
            vg = dsqrt(vg) / dble(ig)
            write(filenum,'(6x,a,3e15.5,a,i5,a)') &
            '>> RMSE for Grad =',vg,vgmin,vgmax,' among',ig,' data'
          end if
          if(ih.ne.0) then
            vd = dsqrt(vd) / dble(id)
            write(filenum,'(6x,a,3e15.5,a,i5,a)') &
            '>> RMSE for DHes =',vd,vdmin,vdmax,' among',id,' data'
            vh = dsqrt(vh) / dble(ih)
            write(filenum,'(6x,a,3e15.5,a,i5,a)') &
            '>> RMSE for Hess =',vh,vhmin,vhmax,' among',ih,' data'
          end if
          if(iv.ne.0) then
            vv = dsqrt(vv) / dble(iv)
            write(filenum,'(6x,a,3e15.5,a,i5,a)') &
            '>> RMSE for HVec =',vv,vvmin,vvmax,' among',iv,' data'
          end if
        end if

        deallocate(dadd)
        deallocate(dmul)
        end subroutine check_krig
