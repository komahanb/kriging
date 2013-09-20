        subroutine search_krig
        use dimKrig
        implicit none
        include 'mpif.h'
        integer :: ierr
        integer :: ip,nmv
        integer :: iswitch,ister
        character(len=10), dimension(100,2) :: Cswtc
        double precision,  dimension(100)   :: Pswtc
        double precision, dimension(ndim+1) :: bound
        double precision, dimension(ndim)   :: xin,xout
        character(len=30) :: Cfile
        integer :: i,mode,ndes,iadd
        double precision :: vrange
        double precision, allocatable, dimension(:,:) :: design

! Switching Strategy
        if(Cmode(:16).eq.'Search_by_Switch')then
          Cmode = 'Search_by_GA2_Local '
          call read_switch(id_proc,ip,nmv,iswitch,ister,Cswtc,Pswtc)

          if(Cswtc(nmv,1).eq.'MaxVar')then
            if( (ip.gt.iswitch).or. &
                (nmv.ne.1.and.Cswtc(nmv,2).eq.'Go') )then
              Cmode = 'Search_by_MaxVar '
            end if
          end if
          if(id_proc.eq.0)then
            if( (ip.gt.iswitch).or. &
                (nmv.ne.1.and.Cswtc(nmv,2).eq.'Go') )then
              open(10,file='info_switch.dat', &
                      form='formatted',status='unknown')
              write(10,'(i5)')nmv+1
              close(10)
            else
              open(10,file='info_switch.dat', &
                      form='formatted',status='unknown')
              write(10,'(i5)')1 ! clean-up
              close(10)
            end if
            write(filenum,'(1x,3a,2i5)') &
            '>> Strategy of Switching : ',Cmode,' by ',ip,nmv
          end if
          call MPI_BCAST(Cmode,20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        end if

!  Pre-Process for Each Strategy
        if(     Cmode(:14).eq.'Search_by_GA1 ')then
          call read_iran(id_proc,iran)
          Cfile = 'paramGAEI.inp'
          mode  = 11
          ndes  = nOPT
        else if(Cmode(:14).eq.'Search_by_GA2 ')then
          call read_iran(id_proc,iran)
          Cfile = 'paramGAEI.inp'
          mode  = 12
          ndes  = nOPT*2
        else if(Cmode(:19).eq.'Search_by_GA1_Local')then
          call read_iran(id_proc,iran)
          call read_vran(id_proc,vrange)
          Cfile = 'paramGAEI.inp'
          mode  = 13
          ndes  = nOPT
        else if(Cmode(:19).eq.'Search_by_GA2_Local')then
          call read_iran(id_proc,iran)
          call read_vran(id_proc,vrange)
          Cfile = 'paramGAEI.inp'
          mode  = 14
          ndes  = nOPT*2
        else if(Cmode(:20).eq.'Search_by_GA2_and_GM')then
          call read_iran(id_proc,iran)
          call read_vran(id_proc,vrange)
          Cfile = 'paramGAEI.inp'
          mode  = 14
          ndes  = nOPT*2
        else if(Cmode(:16).eq.'Search_by_MaxVar')then
          call read_iran(id_proc,iran)
          call read_vran(id_proc,vrange)
          Cfile = 'paramGAEI.inp'
          mode  = 15
          ndes  = nOPT
        else if(Cmode(:13).eq.'Search_by_GM1 ')then
          ndes  = nOPT
        else if(Cmode(:13).eq.'Search_by_GM2 ')then
          ndes  = nOPT
        else
          if(id_proc.eq.0)then
            write(*,*)'Unknown Cmode = ',Cmode
          end if
          call stop_all
        end if

        allocate( design(ndes,ndim) )
        design = -1.d0
        call find_Optimal
        call read_all_krig
        bound         = 1.d0 ! dum
        bound(ndim+1) = 0.d0 ! dum

        if(id_proc.eq.0)then
           write(filenum,'(1x,a,a)')'>> ',Cmode
           if(mode.eq.13.or.mode.eq.14.or.mode.eq.15)then
           write(filenum,'(6x,a,2(e15.5,a))')'>> Search around Optimal, [', &
                -.5*vrange,' :',.5*vrange,' ]'
           end if
        end if
        if(     Cmode(:12).eq.'Search_by_GA')then
          call Genetic_Algorithm(mode,id_proc,num_proc,       &
                                 ndim,tdim,iscf,iran,         &
                                 ndebug,Cfile,                &
                                 bound,                       &
                                 nOPT,nfunc,vrange,           &
                                 ndes,ndim,design )

          if(Cmode(:20).eq.'Search_by_GA2_and_GM')then
            do i=1,ndes
              if(mod(i,2).eq.1)then
                iadd = 20 + (i+1)/2
              else
                iadd = 10 + (i  )/2 
              end if
              xin(:) = design(i,:)
              call BFGS(iadd,id_proc,num_proc, &
                        ndim,tdim,iscf,ndebug,50, &
                        bound, &
                        xin,xout )
              design(i,:) = xout(:)
            end do
          end if

        else if(Cmode(:12).eq.'Search_by_GM')then
          if(Cmode(13:13).eq.'1')then
            iadd = 10
          else if(Cmode(13:13).eq.'2')then
            iadd = 20
          else
            stop'unknown Cmode for GM'
          end if
          do i=1,nOPT
          xin(:) = sample(IDopt(i),:)
          call BFGS(iadd,id_proc,num_proc, &
                    ndim,tdim,iscf,ndebug,50, &
                    bound, &
                    xin,xout )
          design(i,:) = xout(:)
          end do

        else if(Cmode(:16).eq.'Search_by_MaxVar')then
          call Genetic_Algorithm(mode,id_proc,num_proc,       &
                                 ndim,tdim,iscf,iran,         &
                                 ndebug,Cfile,                &
                                 bound,                       &
                                 nOPT,nfunc,vrange,           &
                                 ndes,ndim,design )

        end if

        call Output_Design(mode,ndes,ndim,design)

        deallocate(design)
        call deallocate_all_krig
        end subroutine search_krig
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine read_all_krig
        use dimKrig
        implicit none
        integer :: ifunc,igrad,rmax,nmax,kmax,hmax,tmax
        integer :: i,j,k,l,dct,jct,cct,indx,itarg,ict
        integer :: iflag,istat
        double precision :: verr,ymin,ymax,yhat,rmse,EI,yexa
        double precision, dimension(ndim) :: xin,yhatprime

        i_dim = 20
        allocate( ipv(i_dim,nfunc-1),stat=istat )
        allocate( iak(i_dim,nfunc-1),stat=istat )
        allocate( ibl(i_dim,nfunc-1),stat=istat )
        allocate( icm(i_dim,nfunc-1),stat=istat )
        ipv   = 0
        iak   = 0
        ibl   = 0
        icm   = 0

        rmax = 0
        nmax = 0
        kmax = 0
        hmax = 0
        tmax = 0
        do ifunc = 1,nfunc-1
          call read_krig(2,ifunc)
          rmax = max(rmax,rsample)
          nmax = max(nmax,nsize)
          kmax = max(kmax,kreg)
          hmax = max(hmax,lmax)
          tmax = max(tmax,tdim)
          ipv(1,ifunc) = rsample
          ipv(2,ifunc) = nsize
          ipv(3,ifunc) = kreg
          ipv(4,ifunc) = lmax
          ipv(5,ifunc) = tdim
        end do
        if(tmax.ne.(hmax*(hmax+1)/2))stop'tmax.and.hmax in all_krig'

        d_dim = ndim*rmax + kmax                  & ! sampl,mean
              + 1         + hmax                  & ! devi,deviratio
              + ndim*tmax + ndim*tmax             & ! theta,power
              + nmax*nmax + nmax                  & ! Rinv,yy
              + kmax*nmax + nmax      + nmax*kmax & ! FB,RYFB,FR
              + kmax*kmax + 2         + rmax      & ! FRFinv,ymax/min,dxadd
              + rmax*ndim*2                         ! hvec
        j_dim = 11*11     + rmax                    ! ict_sample,nparent
        c_dim = rmax                                ! info

        iflag = 0
        allocate( dak(d_dim,nfunc-1),stat=istat ) ! with iak
        if(istat.ne.0)iflag = iflag + 1
        allocate( jbl(j_dim,nfunc-1),stat=istat ) ! with ibl
        if(istat.ne.0)iflag = iflag + 1
        allocate( ccm(j_dim,nfunc-1),stat=istat ) ! with icm
        if(istat.ne.0)iflag = iflag + 1
        if(nmax.gt.nsample)then
          deallocate(sampl,nparent,dxadd,inf)
          allocate( sampl(nmax,ndim),stat=istat )
          if(istat.ne.0)iflag = iflag + 1
          allocate( nparent(nmax),stat=istat )
          if(istat.ne.0)iflag = iflag + 1
          allocate( dxadd(nmax),stat=istat )
          if(istat.ne.0)iflag = iflag + 1
          allocate( inf(nmax),stat=istat )
          if(istat.ne.0)iflag = iflag + 1
          sampl = 0.d0
          dxadd = 0.d0
          nparent = 0
          inf     = '      '
        end if
        dak = 0.d0
        jbl = 0
        ccm = '      '
        if(iflag.ne.0)stop'allocation in read_all_krig'

        rsample = rmax
        nsize   = nmax
        kreg    = kmax
        lmax    = hmax
        call allocate_krig(-1)
        if(id_proc.eq.0)then
        write(filenum,*)'>> Make "dak" Database for metamodel'
        end if

        do ifunc = 1,nfunc-1
          rsample =  ipv(1,ifunc)
          nsize   =  ipv(2,ifunc)
          kreg    =  ipv(3,ifunc)
          lmax    =  ipv(4,ifunc)
          tdim    =  ipv(5,ifunc)
          call read_krig(0,ifunc)

! make databse of dak
          dct = 0
          jct = 0
          cct = 0
! sampl
          iak(1,ifunc) = dct + 1
          do k=1,ndim
          do i=1,rsample
             dct = dct + 1
             dak(dct,ifunc) = sampl(i,k)
          end do
          end do
          if(dct.gt.d_dim)stop'dct > d_dim1'
! mean
          iak(2,ifunc) = dct + 1
          do k=1,kreg
             dct = dct + 1
             dak(dct,ifunc) = mean(k,1)
          end do
          if(dct.gt.d_dim)stop'dct > d_dim2'
! devi
          iak(3,ifunc) = dct + 1
          dct = dct + 1
          dak(dct,ifunc) = devi
          if(dct.gt.d_dim)stop'dct > d_dim3'
! deviratio
          iak(4,ifunc) = dct + 1
          if(lmax.ne.1)then
           do k=1,lmax
            dct = dct + 1
            dak(dct,ifunc) = deviratio(k)
           end do
          end if
          if(dct.gt.d_dim)stop'dct > d_dim4'
! theta
          iak(5,ifunc) = dct + 1
          do l=1,tdim
          do k=1,ndim
             dct = dct + 1
             dak(dct,ifunc) = d_theta(k,l)
          end do
          end do
          if(dct.gt.d_dim)stop'dct > d_dim5'
! power
          iak(6,ifunc) = dct + 1
          if(iscf.eq.0.or.iscf.eq.1)then
            do l=1,tdim
            do k=1,ndim
             dct = dct + 1
             dak(dct,ifunc) = d_power(k,l)
            end do
            end do
            if(dct.gt.d_dim)stop'dct > d_dim6'
          end if
! Rinv
          iak(7,ifunc) = dct + 1
          do k=1,nsize
          do i=1,nsize
             dct = dct + 1
             dak(dct,ifunc) = Rinv(i,k)
          end do
          end do
          if(dct.gt.d_dim)stop'dct > d_dim7'
! yy
          iak(8,ifunc) = dct + 1
          do i=1,nsize
             dct = dct + 1
             dak(dct,ifunc) = yy(i)
          end do
          if(dct.gt.d_dim)stop'dct > d_dim8'
! FB
          iak(9,ifunc) = dct + 1
          do k=1,kreg
          do i=1,nsize
             dct = dct + 1
             dak(dct,ifunc) = FB(i,k)
          end do
          end do
          if(dct.gt.d_dim)stop'dct > d_dim9'
! RYFB
          iak(10,ifunc) = dct + 1
          do i=1,nsize
             dct = dct + 1
             dak(dct,ifunc) = RYFB(i,1)
          end do
          if(dct.gt.d_dim)stop'dct > d_dim10'
! FR
          iak(11,ifunc) = dct + 1
          do k=1,nsize
          do i=1,kreg
             dct = dct + 1
             dak(dct,ifunc) = FR(i,k)
          end do
          end do
          if(dct.gt.d_dim)stop'dct > d_dim11'
! FRFinv
          iak(12,ifunc) = dct + 1
          do k=1,kreg
          do i=1,kreg
             dct = dct + 1
             dak(dct,ifunc) = FRFinv(i,k)
          end do
          end do
          if(dct.gt.d_dim)stop'dct > d_dim12'
! ymin/ymax
          ymax = -1.d10
          ict  =  0
          do i=1,rsample
            if(inf(i)(1:2).eq.'1_')then
              ict = ict + 1
              if(yy(ict).gt.ymax)ymax = yy(ict)
            end if
          end do
          iak(13,ifunc) = dct + 1
          dct = dct + 1
          dak(dct,ifunc) = ymax

          iak(14,ifunc) = dct + 1
          dct = dct + 1
          itarg = 0
          do i=1,nOPT
            if(ifunc.eq.nfOPT(i))itarg = i
          end do
          if(itarg.ne.0)then
            dak(dct,ifunc) = OFopt(itarg)
          else
            dak(dct,ifunc) = yy(ict) ! dummy
          end if
          if(dct.gt.d_dim)stop'dct > d_dim14'
! dxadd
          iak(15,ifunc) = dct + 1
          do i=1,rsample
             dct = dct + 1
             dak(dct,ifunc) = dxadd(i)
          end do
          if(dct.gt.d_dim)stop'dct > d_dim15'
! hvec
          iak(16,ifunc) = dct + 1
          do k=1,2
          do j=1,ndim
          do i=1,rsample
             dct = dct + 1
             dak(dct,ifunc) = hvec(i,j,k)
          end do
          end do
          end do
          if(dct.gt.d_dim)stop'dct > d_dim16'

! integer info
! ict_sample
          ibl(1,ifunc) = jct + 1
          do j=0,10
          do i=0,10
             jct = jct + 1
             jbl(jct,ifunc) = ict_sample(i,j)
          end do
          end do
! nparent
          ibl(2,ifunc) = jct + 1
          do i=1,rsample
             jct = jct + 1
             jbl(jct,ifunc) = nparent(i)
          end do
! character info
! inf
          icm(1,ifunc) = cct + 1
          do i=1,rsample
             cct = cct + 1
             ccm(cct,ifunc) = inf(i)
          end do

          if(dct.gt.d_dim)stop'dct.gt.d_dim in read_all_krig'
          if(jct.gt.j_dim)stop'jct.gt.j_dim in read_all_krig'
          if(cct.gt.c_dim)stop'jct.gt.j_dim in read_all_krig'
        end do ! loop of ifunc
! ipv( 1) = rsample
!      2) = nsize
!      3) = kreg 
!      4) = lmax 
!      5) = tdim
! iak( 1) = index for sampl
!      2) = index for mean
!      3) = index for devi
!      4) = index for deviratio
!      5) = index for theta
!      6) = index for power
!      7) = index for Rinv
!      8) = index for yy
!      9) = index for FB
!     10) = index for RYFB
!     11) = index for FR
!     12) = index for FRFinv
!     13) = index for ymax
!     14) = index for ymin
!     15) = index for dxadd
!     16) = index for hvec
! ibl( 1) = index for ict_sample
!      2) = index for nparent
! icm( 1) = index for inf
!     to add new info, you should modify read/write_krig,
!     and i_dim/d_dim/j_dim should be increased.
!     Then make new connection with meta_estimation

        call deallocate_krig(-1)
        call calc_correction_term ! dadd/dmul
        return

! check RMSE
        do i = 1,nfunc-1
          rsample = iak(1,i)
          verr    = 0.d0
          ict     = 0
          do 900 j=1,rsample
            if(ccm(icm(1,i)+(j-1),i)(1:2).ne.'1_')go to 900
            ict = ict + 1
            do k=1,ndim
              indx   = iak(1,i) + (k-1)*rsample + (j-1)
              xin(k) = dak(indx,i)
            end do
            ! here ymin is dummy
            call meta_estimation(0,xin,yhat,yhatprime,RMSE,EI,          &
!                mEI,mLC,id_proc,ndim,kreg_orig,iscf,mode_dck,          &
                 mEI,  0,id_proc,ndim,kreg_orig,iscf,mode_dck,          &
                 ccrf,                                                  &
                 ipv(1,i),ipv(2,i),ipv(3,i),ipv(4,i),ipv(5,i),          &
                 dak(iak( 1,i),i),dak(iak( 2,i),i),dak(iak( 3,i),i),    &
                 dak(iak( 4,i),i),dak(iak( 5,i),i),dak(iak( 6,i),i),    &
                 dak(iak( 7,i),i),dak(iak(10,i),i),dak(iak(11,i),i),    &
                 dak(iak(12,i),i),dak(iak(13,i),i),dak(iak(14,i),i),    &
                 dak(iak(15,i),i),dak(iak(16,i),i),                     &
                 jbl(ibl( 1,i),i),jbl(ibl( 2,i),i),                     &
                 ccm(icm( 1,i),i),                                      &
                 ncr1,ncr2,dadd(1,i),dmul(1,i)  )
            yexa = dak(iak(8,i)+ict-1,i)
            verr = verr + (yexa-yhat)**2
!           write(*,'(i5,2e15.5)')j,yexa,yhat
900       continue
          verr = dsqrt(verr)/dble(ict)
          write(filenum,'(6x,a,e15.5,a,i6,a)') &
          '>> RMSE with dak = ',verr,' among',ict,' pts'
        end do
        call deallocate_all_krig
        call Deallocating
        call stop_all

        end subroutine read_all_krig
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine deallocate_all_krig
        use dimKrig
        implicit none
        deallocate(iak,dak)
        deallocate(ibl,jbl)
        deallocate(icm,ccm)
        deallocate(dadd,dmul)
        end subroutine deallocate_all_krig
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine DV2DESV2(xin,vrange,xout)
        use dimKrig
        implicit none
        double precision, dimension(ndim), intent(in)  :: xin
        double precision,                  intent(in)  :: vrange
        double precision, dimension(ndim), intent(out) :: xout
        integer :: i

        do i=1,ndim
          if(nOPT.eq.1.and.IDopt(1).ge.1.and.IDopt(1).le.nsample)then
             xout(i) = sample(IDopt(1),i) + vrange*(xin(i)-0.5d0)
          else
             xout(i) = xin(i)
          end if
          if(xout(i).lt.0.d0)xout(i) = 0.d0
          if(xout(i).gt.1.d0)xout(i) = 1.d0
        end do
        
        end subroutine DV2DESV2
