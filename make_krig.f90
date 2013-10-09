        subroutine Make_Krig
        use dimKrig
        implicit none
        include 'mpif.h'
        integer :: ierr,ip,nmv
        integer :: ifunc,igrad
        double precision :: llfd,verr
        double precision, dimension(ndim,tdim) :: theta,power
        double precision, dimension(2*ndim) :: dinit
        double precision, dimension(ndim+1) :: bound
        integer :: iswitch,ister
        character(len=10), dimension(100,2) :: Cswtc
        double precision, dimension(100)    :: Pswtc
        integer :: i,l,idum,ndes,ndv,nlim,mode
        double precision :: dum,dlim
        double precision, allocatable, dimension(:,:) :: design
        character(len=30) :: Cfile

! Switching Strategy
        dlim = 0.d0 ! Then nothing happen
        nlim = 0
        if(Cmode(:14).eq.'Make_by_Switch')then
!         Cmode = 'Make_by_GA '
          Cmode = 'Make_by_GA1 '
!         Cmode = 'Make_with_Adjust '
          call read_switch(id_proc,ip,nmv,iswitch,ister,Cswtc,Pswtc)

          if(Cswtc(nmv,1).eq.'Local')then
            if( (ip.gt.iswitch).or. &
                (nmv.ne.1.and.Cswtc(nmv,2).eq.'Go') )then
              dlim = Pswtc(nmv)
              if(id_proc.eq.0)then
                write(filenum,'(1x,1a,2i5,f8.4)') &
                '>> Strategy of Switching : Local Kriging by ',ip,nmv,dlim
              end if
            end if
          else if(Cswtc(nmv,1).eq.'Limit')then
            if( (ip.gt.iswitch).or. &
                (nmv.ne.1.and.Cswtc(nmv,2).eq.'Go') )then
              Cmode = 'Make_by_GA '
              nlim  = int(Pswtc(nmv))
              if(nlim.ge.1000)    stop'too big   nlim'
              if(id_proc.eq.0)then
                write(filenum,'(1x,1a,2i5,i5)') &
                '>> Strategy of Switching : Limit Kriging by ',ip,nmv,nlim
              end if
            end if
          end if
          call MPI_BCAST(Cmode,20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
          ! vrange will be updated at reduce_sample
        end if

        ndes = 1
        ndv  = tdim*ndim
        allocate( design(ndes,ndv) )

! Pre-Process for Each Strategy
        if(     Cmode(:11).eq.'Make_by_GA ')then
          call read_iran(id_proc,iran)
          Cfile = 'paramGALN.inp'
        else if(Cmode(:11).eq.'Make_by_GA1')then
          call read_iran(id_proc,iran)
          Cfile = 'paramGALN1.inp'
        else if(Cmode(:11).eq.'Make_by_GM ')then
        else if(Cmode(:17).eq.'Make_by_GA_and_GM')then
          call read_iran(id_proc,iran)
          Cfile = 'paramGALN.inp'
        else if(Cmode(:18).eq.'Make_from_Previous')then
        else if(Cmode(:18).eq.'Make_with_Previous')then
        else if(Cmode(:18).eq.'Make_with_Constant')then
        else if(Cmode(:16).eq.'Make_with_Adjust')then
        else
          if(id_proc.eq.0)then
            write(*,*)'*Unknown Cmode = ',Cmode
          end if
          call stop_all
        end if

        call find_Optimal

        do 100 ifunc=1,nfunc-1

           call reduce_data(ifunc,igrad,dlim,nlim)
           if(igrad.eq.2)call indirect(ifunc)
           call allocate_krig(igrad)
           call make_training_data(igrad)
           call make_index_Rij

           bound              = 10.d0
           if(nlim.ne.0)bound = 50.d0
           bound(ndim+1)      = 0.01d0
           deviratio_best     = 1.d0
           llfd_best          =-1.d10

           if(     Cmode(:11).eq.'Make_by_GA '.or. &
                   Cmode(:11).eq.'Make_by_GA1')then
!            call read_krig(1,ifunc)
!            call make_init(id_proc,ndim,tdim,iscf,nsize,bound, &
!                           d_theta,d_power)
             cpena = 0.d0 ! x2 at the end
             mode  = 1    ! for level*n dvs
             mode  = 2    ! for       n dvs
             if(Cmode(:11).eq.'Make_by_GA1')mode = 3 ! for 1dv
             call Genetic_Algorithm(mode,id_proc,num_proc,            &
                                    ndim,tdim,iscf,iran,ndebug,Cfile, &
                                    bound,                            &
                                    nOPT,nfunc,dum,                   &
                                    ndes,ndv,design)
             if(iscf.eq.1)then ! General Gaussian
               do l=1,tdim
               do i=1,ndim
                  theta(i,l) = design(1,          (l-1)*ndim+i)
                  power(i,l) = design(1,tdim*ndim+(l-1)*ndim+i)
               end do
               end do
             else
               do l=1,tdim
               do i=1,ndim
                  theta(i,l) = design(1,(l-1)*ndim+i)
                  power(i,l) = 2.d0
               end do
               end do
             end if

           else if(Cmode(:11).eq.'Make_by_GM ')then
             cpena = 10.d0 ! adjust the weight of constraint
             cpena = -1.d0 ! anyway proceed
151          continue
             call read_thpw(0,ndim,tdim,iscf,theta,power)
             mode = 1 ! for level*n dvs
             mode = 2 ! for       n dvs
             mode = 3 ! for 1 dv
             call BFGS(mode,id_proc,num_proc, &
                       ndim,tdim,iscf,ndebug,50, &
                       bound, &
                       theta,power )

           else if(Cmode(:17).eq.'Make_by_GA_and_GM')then
             cpena = 0.d0 ! x2 at the end
             mode = 1 ! for level*n dvs
             mode = 2 ! for       n dvs
             mode = 3 ! for 1 dv

             call Genetic_Algorithm(mode,id_proc,num_proc,            &
                                    ndim,tdim,iscf,iran,ndebug,Cfile, &
                                    bound,                            &
                                    nOPT,nfunc,dum,                   &
                                    ndes,ndv,design)
152          continue
             if(iscf.eq.1)then
               do l=1,tdim
               do i=1,ndim
                 theta(i,l) = design(1,          (l-1)*ndim+i)
                 power(i,l) = design(1,tdim*ndim+(l-1)*ndim+i)
               end do
               end do
             else
               do l=1,tdim
               do i=1,ndim
                 theta(i,l) = design(1,(l-1)*ndim+i)
                 power(i,l) = 2.d0
               end do
               end do
             end if
             cpena = 10.d0 ! adjust the weight of constraint
             cpena = -1.d0 ! anyway proceed
             mode = 1 ! for level*n dvs
             mode = 2 ! for       n dvs
             mode = 3 ! for 1 dv
             call BFGS(mode,id_proc,num_proc, &
                       ndim,tdim,iscf,ndebug,50, &
                       bound, &
                       theta,power )

           else if(Cmode(:18).eq.'Make_from_Previous')then
             cpena = 10.d0 ! adjust the weight of constraint
             cpena = -1.d0 ! anyway proceed
153          continue
             call read_krig(1,ifunc)
             theta = d_theta
             power = d_power
             call BFGS(mode,id_proc,num_proc, &
                       ndim,tdim,iscf,ndebug,50, &
                       bound, &
                       theta,power )

           else if(Cmode(:18).eq.'Make_with_Previous')then
             call read_krig(1,ifunc)
             theta   = d_theta
             power =   d_power

           else if(Cmode(:18).eq.'Make_with_Constant')then
             cpena = -1.d0 ! awyway proceed
             call read_thpw(1,ndim,tdim,iscf,theta,power)

           else if(Cmode(:16).eq.'Make_with_Adjust')then
             cpena = -2.d0 ! adjust by addint 0.1/0.5
             theta = 0.1d0
             power = 2.0d0

           else
             if(id_proc.eq.0)then
               write(*,*)'*Unknown Cmode = ',Cmode
             end if
             call stop_all
           end if

110        continue
           call likelihood(theta,power,llfd,verr)
           if(id_proc.eq.0)then
            do l=1,tdim
             write(filenum,'(1x,a,99f8.3)')'>> Theta      = ',(theta(i,l),i=1,ndim)
             if(iscf.eq.0.or.iscf.eq.1) &
             write(filenum,'(1x,a,99f8.3)')'>> Power      = ',(power(i,l),i=1,ndim)
            end do
            write(filenum,'(1x,a,1e15.5)')'>> Deviation  = ',devi
            do l=2,lmax
            write(filenum,'(1x,a,i1,a,2e15.5)')'>> Deviratio',l,' = ', &
                                    deviratio(l),deviratio_best(l)
            end do
            write(filenum,'(1x,2(a,e15.5))') &
            '>> Likelihood = ',llfd,' with Error = ',verr
           end if

           if(verr.ge.1.d-3)then ! Consideration for Constraint
             if(cpena.eq.0.d0)then ! for GA
               do l=1,tdim
               do i=1,ndim
                 if(theta(i,l).gt.bound(i))go to 120 ! fin
               end do
               end do
               theta   = theta   * 2.d0
               if(id_proc.eq.0) &
               write(*,'(6x,a)')'>> Wrong Model, x2 for Theta'
               go to 110

             else if(cpena.eq.-1.d0)then ! for grad/constant
               if(id_proc.eq.0) &
               write(filenum,'(6x,a)')'>> Wrong Model, but proceeds'

             else if(cpena.eq.-2.d0)then ! adjust
               if(dlim.eq.0.d0.and.theta(1,1).ge.3.d0)go to 120
               if(                 theta(1,1).ge.3.d0)go to 120
               if(theta(1,1).lt.1.5d0)then
                 theta   = theta   + 0.1d0
                 if(id_proc.eq.0) &
                 write(*,'(6x,a)')'>> Wrong Model, +0.1 for Theta'
               else
                 theta   = theta   + 0.5d0
                 if(id_proc.eq.0) &
                 write(*,'(6x,a)')'>> Wrong Model, +0.5 for Theta'
               end if
               go to 110

             else ! for gradient
               if(cpena.ge.1.d6)go to 120 ! fin
               cpena = cpena * 10.d0
               if(id_proc.eq.0) &
               write(*,'(6x,a,e12.3)')'>> Wrong Model, x10 for Cpena',cpena
               if(     Cmode(:10).eq.'Make_by_GM')then
                 go to 151
               else if(Cmode(:17).eq.'Make_by_GA_and_GM')then
                 go to 152
               else if(Cmode(:18).eq.'Make_from_Previous')then
                 go to 153
               else
                 if(id_proc.eq.0)then
                   write(*,*)'*Unknown Cmode for Cpena = ',Cmode
                 end if
                 call stop_all
               end if

             end if
           end if
120        continue

           call write_krig(ifunc,theta,power)
           call check_krig(ifunc,theta,power)
           call deallocate_krig(igrad)

100     continue ! ifunc loop

        deallocate(design)
        end subroutine Make_Krig
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine read_thpw(mode,ndim,tdim,iscf,theta,power)
        implicit none
        integer, intent(in) :: mode,ndim,tdim,iscf
        double precision, dimension(ndim,tdim), intent(out) :: theta,power
! mode=0 : for BFGS
! mode=1 : for Decision
        integer :: i,l

        open(100,file='thetainit.dat',form='formatted',status='old',err=200)
        do l=1,tdim
          read(100,*,end=200,err=200)(theta(i,l),i=1,ndim)
          if(iscf.eq.1)then
            read(100,*,end=200,err=200)(power(i,l),i=1,ndim)
          end if
        end do
        close(100)
        return

200     continue
!       write(*,*)'*No thetainit.dat'
        theta = 1.d0
        if(iscf.eq.1)then
          if(mode.eq.0)then
            power = 1.d0
          else if(mode.eq.1)then
            power = 2.d0
          else
            stop'unknown mode in read_thpw'
          end if
        else
          power = 2.d0
        end if

        end subroutine read_thpw
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine read_iran(id_proc,iran)
        implicit none
        include 'mpif.h'
        integer, intent(in)  :: id_proc
        integer, intent(out) :: iran
        integer :: ierr

        if(id_proc.eq.0)then
          open(10,file='iran.dat',form='formatted',status='old',err=100)
          read(10,*)iran
          close(10)
          go to 200
100       continue
          iran = 1
200       continue
          open(10,file='iran.dat',form='formatted',status='unknown')
          write(10,*)iran+1
          close(10)
        end if
        call MPI_BCAST(iran,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        end subroutine read_iran
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine read_vran(id_proc,vrange)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: id_proc
        double precision, intent(out) :: vrange
        integer :: ierr

        if(id_proc.eq.0)then
          open(10,file='vrange.inp',form='formatted',status='old')
          read(10,*)vrange
          close(10)
          if(vrange.le.0.d0.or.vrange.ge.2.d0)vrange = 2.d0
        end if
        call MPI_BCAST(vrange,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        end subroutine read_vran
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine DV2DESV1(neva,ndim,tdim,ndv,bound,xin,theta,power)
        implicit none
        integer, intent(in) :: neva,ndim,tdim,ndv
        double precision, dimension(ndim+1),    intent(in) :: bound
        double precision, dimension(ndv),       intent(in) :: xin
        double precision, dimension(ndim,tdim), intent(out) :: theta,power
        integer :: i,l,ilog
        double precision :: fac
 
        ilog = 1
        if(ilog.eq.1)then
         ! log distribution
         do i=1,ndim
           fac = log10( bound(i)/bound(ndim+1) )
           if(fac.lt.1.or.fac.gt.5.)stop'fac in DV2DESV1'
           if(neva.eq.1.or.neva.eq.2)then      ! for general cases
             do l=1,tdim
               theta(i,l) = bound(ndim+1) * ( 10.d0**(xin((l-1)*ndim+i)*fac) )
               if(neva.eq.1)then
                 power(i,l) = 2.d0
               else if(neva.eq.2)then ! general Gaussian
                 power(i,l) = 1.d0 + xin(ndim*tdim+(l-1)*ndim+i)*(2.d0-1.d0)
               end if
             end do
           else if(neva.eq.3.or.neva.eq.4)then ! for simplified VFM
             do l=1,tdim
               theta(i,l) = bound(ndim+1) * ( 10.d0**(xin(i)*fac) )
               if(neva.eq.3)then
                 power(i,l) = 2.d0
               else if(neva.eq.4)then ! general Gaussian
                 power(i,l) = 1.d0 + xin(ndim+i)*(2.d0-1.d0)
               end if
             end do
           else if(neva.eq.5.or.neva.eq.6)then ! for 1DV cases
             do l=1,tdim
               theta(i,l) = bound(ndim+1) * ( 10.d0**(xin(1)*fac) )
               if(neva.eq.5)then
                 power(i,l) = 2.d0
               else if(neva.eq.6)then ! general Gaussian
                 power(i,l) = 1.d0 + xin(2)*(2.d0-1.d0)
               end if
             end do
           else
             stop'unknown set of ndim/ndv in DV2DESV1'
           end if
         end do
        else
         ! general distribution
         do i=1,ndim
           fac = bound(i)-bound(ndim+1)
           if(neva.eq.1.or.neva.eq.2)then      ! for general cases
             do l=1,tdim
               theta(i,l) = bound(ndim+1) + xin((l-1)*ndim+i)*fac
               if(neva.eq.1)then
                 power(i,l) = 2.d0
               else if(neva.eq.2)then ! general Gaussian
                 power(i,l) = 1.d0 + xin(ndim*tdim+(l-1)*ndim+i)*(2.d0-1.d0)
               end if
             end do
           else if(neva.eq.3.or.neva.eq.4)then ! for simplified VFM
             do l=1,tdim
               theta(i,l) = bound(ndim+1) + xin(i)*fac
               if(neva.eq.1)then
                 power(i,l) = 2.d0
               else if(neva.eq.2)then ! general Gaussian
                 power(i,l) = 1.d0 + xin(ndim+i)*(2.d0-1.d0)
               end if
             end do
           else if(neva.eq.5.or.neva.eq.6)then ! for 1DV cases
             do l=1,tdim
               theta(i,l) = bound(ndim+1) + xin(1)*fac
               if(neva.eq.1)then
                 power(i,l) = 2.d0
               else if(neva.eq.2)then ! general Gaussian
                 power(i,l) = 1.d0 + xin(2)*(2.d0-1.d0)
               end if
             end do
           else
             stop'unknown set of ndim/ndv in DV2DESV1'
           end if
         end do
        end if

        end subroutine DV2DESV1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine make_init(id_proc,ndim,tdim,iscf,nsize,bound, &
                             theta,power)
        implicit none
        integer, intent(in) :: id_proc,ndim,iscf,nsize,tdim
        double precision, dimension(ndim+1), intent(in) :: bound
        double precision, dimension(ndim,tdim), intent(in) :: &
        theta,power
        integer :: k,ilog,l
        double precision :: fac
        double precision, dimension(ndim,tdim) :: th,po

        ilog = 1
        do l=1,tdim
        do k=1,ndim
          if(theta(  k,l).lt.bound(ndim+1).or. &
             theta(  k,l).gt.bound(k)     ) go to 100
          if(ilog.eq.1)then
            fac    = log10( bound(  k)/bound(ndim+1) )
            th(k,l) = log10( theta(k,l)/bound(ndim+1) ) / fac
          else
            th(k,l) = (theta(k,l)-bound(ndim+1))/(bound(k)-bound(ndim+1))
          end if
        end do
        end do

        if(iscf.eq.0.or.iscf.eq.1)then
         do k=1,ndim
         do l=1,tdim
          if(power(k,l).lt.1.d0.or. &
             power(k,l).gt.2.d0        ) go to 100
          po(k,l) = (power(k,l)-1.d0)/(2.d0-1.d0)
         end do
         end do
        end if

        if(id_proc.eq.0)then
          open(11,file='initLN.inp',form='formatted',status='unknown')
          if(iscf.eq.1)then
            write(11,'(999f15.8)') &
            ((th(k,l),k=1,ndim),l=1,tdim),((po(k,l),k=1,ndim),l=1,tdim)
          else
            write(11,'(999f15.8)') &
            ((th(k,l),k=1,ndim),l=1,tdim)
          end if
          close(11)
        end if
100     continue
        return
        end subroutine make_init
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine read_switch(id_proc,ip,nmv,iswitch,ister,Cswtc,Pswtc)
        implicit none
        integer, intent(in)  :: id_proc
        integer, intent(out) :: iswitch,ister,ip,nmv
        character(len=10), dimension(100,2), intent(out) :: Cswtc
        double precision, dimension(100), intent(out)    :: Pswtc
        integer :: i

        open(10,file='switch.inp',form='formatted',status='old')
        read(10,*)iswitch
        read(10,*)ister
        if(ister.ge.100)stop'ister.ge.100'
        do i=1,ister
          if(i.eq.1)then
            read(10,*)Cswtc(i,1)
            if(Cswtc(i,1).eq.'Local'.or.Cswtc(i,1).eq.'Limit')then
              backspace 10
              read(10,*)Cswtc(i,1),Pswtc(i)
            end if
          else
            read(10,*)Cswtc(i,1),Cswtc(i,2)
            if(Cswtc(i,1).eq.'Local'.or.Cswtc(i,1).eq.'Limit')then
              backspace 10
              read(10,*)Cswtc(i,1),Cswtc(i,2),Pswtc(i)
            end if
          end if
          if(Cswtc(i,1).ne.'Usual'.and. &
             Cswtc(i,1).ne.'Local'.and. &
             Cswtc(i,1).ne.'Limit'.and. &
             Cswtc(i,1).ne.'MaxVar') then
             if(id_proc.eq.0)write(*,'(2a)')'*Cswtc1 = ',Cswtc(i,1)
             call stop_all
          end if
          if(i.ne.1.and.Cswtc(i,2).ne.'Go'.and.Cswtc(i,2).ne.'Stop') then
             if(id_proc.eq.0)write(*,'(2a)')'*Cswtc2 = ',Cswtc(i,2)
             call stop_all
          end if
          if(Cswtc(i,1).eq.'Local')then
             if(Pswtc(i).le.0.d0.or.Pswtc(i).ge.1.d0) &
             stop'Pswtc for Local'
          else if(Cswtc(i,1).eq.'Limit')then
             if(Pswtc(i).le.1) &
             stop'Pswtc for Limit'
          end if
        end do
        close(10)

        open(10,file='appears.dat',form='formatted',status='old',err=11)
        read(10,*)ip
        go to 12
11      continue
        ip = 0
12      continue
        close(10)
        open(10,file='info_switch.dat',form='formatted',status='old',err=13)
        read(10,*,err=13,end=13)nmv
        go to 14
13      continue
        nmv = 1
14      continue
        close(10)
        if(nmv.lt.1.or.nmv.gt.ister+1)stop'undesirable nmv'
        if(nmv.eq.ister+1)nmv = 1

        end subroutine read_switch
