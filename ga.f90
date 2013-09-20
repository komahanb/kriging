      subroutine Genetic_Algorithm(mode,i_proc,n_proc,          & ! common
                                   ndim,tdim,iscf,iran,         & ! common
                                   idebug,Cfile,                & ! common
                                   bound,                       & ! likeli
                                   nOPT,nfunc,vrangein,         & ! newpts
                                   ndes,nsize,design)             ! output
      use dimGA
      use dimKrig, only:filenum
      implicit none
      include 'mpif.h'
      integer, intent(in) :: mode,i_proc,n_proc
      integer, intent(in) :: ndim,tdim,iscf,iran,idebug
      character(len=30), intent(in) :: Cfile
      double precision, dimension(ndim+1), intent(in) :: bound
      integer, intent(in) :: nOPT,nfunc
      double precision, intent(in) :: vrangein
      integer, intent(in) :: ndes,nsize
      double precision, dimension(ndes,nsize), intent(out)  :: design
! mode=1-3 : Likelihood
! other    : New Sample Search
      integer :: generation,istop,i,j,l,g1,n1
      double precision :: ran
      double precision, allocatable, dimension(:) :: xnew

!     call MPI_START
!     call getarg(1,Cfile)
!     if(Cfile.eq.'help'.or.Cfile.eq.'Help'.or.Cfile.eq.'HELP') &
!     call help

      id_proc  = i_proc
      num_proc = n_proc
      ntmp     = ndim
      ttmp     = tdim
      do i=1,abs(iran)
        call random_number(ran)
      end do

      if(mode.ge.1.and.mode.le.3)then
        if(     mode.eq.1)then ! LM by  (n*l) DVs
         if(iscf.ne.1)then
           ndv   = ndim*tdim
           neva  = 1
         else if(iscf.eq.1)then
           ndv   = ndim*2*tdim
           neva  = 2
         end if
        else if(mode.eq.2)then ! LM by   n DVs for VFM
         if(iscf.ne.1)then
           ndv   = ndim
           neva  = 3
         else if(iscf.eq.1)then
           ndv   = ndim*2
           neva  = 4
         end if
        else if(mode.eq.3)then ! LM by   1 DV
         if(iscf.ne.1)then
           ndv   = 1
           neva  = 5
         else if(iscf.eq.1)then
           ndv   = 2
           neva  = 6
         end if
        end if
        nobj    = 1
        Cobj(1) = 'Max'
        ndat    = ndv+nobj+1
        ndebug  = idebug
        allocate(dv_t(ndim,tdim))   ! theta
        allocate(dv_p(ndim,tdim))   ! power
        allocate(bd_t(ndim+1))
        bd_t = bound

      else if(mode.ge.11.and.mode.le.20)then
        ndv  = ndim
        neva = mode
        if(mode.eq.11.or.mode.eq.13)then
          nobj = nOPT
          do i=1,nobj
            Cobj(i) = 'Max' ! EI
          end do
        else if(mode.eq.12.or.mode.eq.14)then
          nobj = nOPT*2
          do i=1,nobj
            if(mod(i,2).eq.1)then
               Cobj(i) = 'Max' !EI
            else
               Cobj(i) = 'Min' !Yhat
            end if
          end do
        else if(mode.eq.15)then
          nobj = nOPT
          do i=1,nobj
            Cobj(i) = 'Max'
          end do
        end if
        ndat    = ndv + nobj + 1 + (nfunc-1)*3
        ndebug  = idebug
        vrange  = vrangein
        allocate(dv_t(ndim,1))
      else
        stop'unknown mode in GA'
      end if

      if(ndes.ne.nobj)stop'ndes/nobj'
      if(nsize.lt.ndv)stop'nsize<ndv'

      call read_param(Cfile)

      generation = 0
      call make_initial_pop
      if(id_proc.eq.0) &
      write(filenum,*)'>> Calling evaluation in GA'
      call evaluation(generation)

! ----------- Main Loop of GA --------------- !
      if(id_proc.eq.0) &
      write(filenum,'(1x,a,i5)')'>> Starting Main Loop of GA',iran
      do 100 generation = 1,ngen
         
         call ranking(generation)
         call sharing(generation)

         call convergence_check(generation,istop)
         if(istop.eq.1)go to 200

         call mating(generation)
         call pregnancy(generation)
         if(mod(generation,ngen_grad).ne.0) &
         call evaluation(generation)

100   continue
      call ranking(ngen+1)
      call sharing(ngen+1)
200   continue
      if(id_proc.eq.0) &
      write(filenum,*)'>> Finishing Main Loop of GA'
! ----------- Main Loop of GA --------------- !

      call output


! Post-Process
      allocate(xnew(ndv))
      if(mode.ge.1.and.mode.le.3)then !----------------- Likelihood
        g1 = elite(1,1)
        n1 = elite(1,2)
        do i=1,ndv
          xnew(i) = d(g1,n1,i)
        end do
        call DV2DESV1(neva,ntmp,ttmp,ndv,bd_t,xnew,dv_t,dv_p)
        if(id_proc.eq.0) &
        write(filenum,'(6x,a,99f8.3)')'>> Best Design = ', &
                                (d(g1,n1,i),i=ndv+1,ndat)
        if(neva.eq.1.or.neva.eq.3.or.neva.eq.5)then 
          do l=1,ttmp
            if(id_proc.eq.0) &
            write(filenum,'(6x,a,99f8.3)')'>> Theta = ', &
                                    (dv_t(i,l),i=1,ndim)
            do i=1,ndim
              design(1,(l-1)*ndim+i) = dv_t(i,l)
            end do
          end do
        else if(neva.eq.2.or.neva.eq.4.or.neva.eq.6)then
          do l=1,ttmp
            if(id_proc.eq.0) &
            write(filenum,'(6x,a,99f8.3)')'>> Theta/Power = ', &
                               (dv_t(i,l),i=1,ndim),(dv_p(i,l),i=1,ndim)
            do i=1,ndim
              design(1,          (l-1)*ndim+i) = dv_t(i,l)
              design(1,ndim*tdim+(l-1)*ndim+i) = dv_p(i,l)
            end do
          end do
        else
          stop'unknown neva at post of GA'
        end if
        deallocate(dv_t,dv_p,bd_t)

      else if(mode.ge.11.and.mode.le.20)then !--------- NewPts
        do i=1,nobj
          g1 = elite(i,1)
          n1 = elite(i,2)
          do j=1,ndv
            xnew(j) = d(g1,n1,j)
          end do
          if(neva.eq.11.or.neva.eq.12)then
            dv_t(:,1) = xnew(:)
          else if(neva.eq.13.or.neva.eq.14)then
            call DV2DESV2(xnew,vrange,dv_t)
          else if(neva.eq.15)then
            call DV2DESV2(xnew,vrange,dv_t)
          end if
          if(id_proc.eq.0)then
          write(filenum,'(6x,a,i2,a,99f8.3)')'>> Best Design (',i,' ) = ', &
          (dv_t(j,1),j=1,ndim),(d(g1,n1,j),j=ndv+1,ndv+nobj+1)
          end if
          do j=1,ndv
            design(i,j) = dv_t(j,1)
          end do
!         do j=1,nobj+1
!           design(i,ndv+j) = d(g1,n1,ndv+j)
!         end do
        end do
        deallocate(dv_t)

      else
        stop'unknown mode at the last of GA'
      end if
      deallocate(xnew)

!     call stop_all
      deallocate(d,rank,Frank,Rrank,parent,pool)
      end subroutine Genetic_Algorithm
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine evaluation(gen)
        use dimGA
        implicit none
        include 'mpif.h'
        integer, intent(in) :: gen
        integer :: ista,iend
        integer, dimension(0:num_proc-1,1000)      :: ispl
        integer, dimension(0:num_proc-1)           :: ispm
        double precision, dimension(npop,ndv)      :: dv
        double precision, dimension(npop,ndat-ndv) :: res
        double precision, dimension(ndv)           :: xnew
        double precision, dimension(ndat-ndv)      :: ynew
        integer :: i,j,id,ierr

        ista = 1
        iend = npop

! decomposition of npop evaluation on nodes
        if(id_proc.eq.0) &
        call spliting(npop,num_proc,ista,ispl,ispm)

        call MPI_BCAST(ispl(0,1),num_proc*1000,MPI_INTEGER,0, &
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ispm(0),num_proc,MPI_INTEGER,0, &
                       MPI_COMM_WORLD,ierr)
!       if(id_proc.eq.num_proc-1)then
!         do i=0,num_proc-1
!           write(*,'(a,i3,a,100i3)') &
!           'Node',i,' :',(ispl(i,j),j=1,ispm(i))
!         end do
!       end if

! Accord of DV in All Nodes
        if(id_proc.eq.0)then
         do i=1,npop
          dv(i,:) = d(gen,i,:)
         end do
        end if
        call MPI_BCAST( &
        dv(1,1),npop*ndv,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! Evaluation
        do i=1,ispm(id_proc)
           id = ispl(id_proc,i)
           if(id.lt.ista.or.id.gt.iend)stop'wrong id in ispl'
           xnew(:) = dv(id,:)
           call eva(1,xnew,ynew)
           res(id,:) = ynew(:)
          !write(*,'(2i3,99f6.2)')id_proc,id,(xnew(j),j=1,ndv),ynew(1)
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

! Sharing of Results
        do i=0,num_proc-1
          if(ispm(i).ne.0)then
            id = ispl(i,1)
            do j=1,ndat-ndv
              call MPI_BCAST( res(id,j), ispm(i), &
              MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
            end do 
!           if(id_proc.eq.0) &
!           write(*,'(i3,100f6.2)')i,(res(j,1),j=ista,iend)
          end if
        end do

! Update dimension d[]
        do i=1,npop
          do j=1,ndat-ndv
            d(gen,i,ndv+j) = res(i,j)
          end do
        end do

! Update deviratio_best
        if(neva.ge.1.and.neva.le.6)then
          call update_deviratio
          call MPI_Barrier(MPI_COMM_WORLD, ierr)
        end if

        end subroutine evaluation
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine eva_1dsearch(xini,wlam,xdir,alp,alpnew,ynew,fnew)
        use dimGA
        implicit none
        double precision, dimension(ndat),     intent(in)  :: xini
        double precision, dimension(nobj,1),   intent(in)  :: wlam
        double precision, dimension(ndv),      intent(in)  :: xdir
        double precision,                      intent(in)  :: alp
        double precision,                      intent(out) :: alpnew,fnew
        double precision, dimension(ndat),     intent(out) :: ynew
        double precision, dimension(ndv) :: dnew
        double precision, dimension(ndat-ndv) :: yeva
        integer :: i
        double precision :: alptmp

        if(alp.lt.0.d0)stop'alp in 1d_search'
        alpnew = alp
        do i=1,ndv
          alptmp  = alp
          dnew(i) = xini(i) + alp*xdir(i)
          if(dnew(i).lt.0.d0)then
             dnew(i) = 0.d0
             alptmp  = (0.d0-xini(i))/(xdir(i))
          else if(dnew(i).gt.1.d0)then
             dnew(i) = 1.d0
             alptmp  = (1.d0-xini(i))/(xdir(i))
          end if
          alpnew = min(alpnew,alptmp)
        end do
        if(alpnew.lt.0.d0.and.alpnew.ge.-1.d-8)alpnew = 0.d0
        if(alpnew.gt.alp.or.alpnew.lt.0.d0)then
           write(*,*)alpnew,alp
           stop'alpnew in eva'
        end if

        call eva(3,dnew,yeva)

! fw
        fnew = 0.d0
        do i=1,nobj
          if(Cobj(i).eq.'Min'.or.Cobj(i).eq.'min')then
             fnew = fnew + wlam(i,1)*yeva(i)
          else if(Cobj(i).eq.'Max'.or.Cobj(i).eq.'max')then
             fnew = fnew - wlam(i,1)*yeva(i)
          end if
        end do
! ynew
        do i=1,ndv
           ynew(    i) = dnew(i)
        end do
        do i=1,ndat-ndv
           ynew(ndv+i) = yeva(i)
        end do

        end subroutine eva_1dsearch
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine evagrad(xnew,gnew)
        use dimGA
        implicit none
        double precision, dimension(ndat),     intent(in)  :: xnew
        double precision, dimension(nobj,ndv), intent(out) :: gnew
        double precision, dimension(ndv)      :: dnew
        double precision, dimension(ndat-ndv) :: ynew
        integer :: i,j,k
        double precision :: dx,ff,fo

        dx   = 1.d-5
        gnew = 0.d0
        do i=1,ndv
           do k=1,ndv
              dnew(k) = xnew(k)
           end do
           dnew(i) = dnew(i) + dx
           call eva(2,dnew,ynew)
           do j=1,nobj
             fo = xnew(ndv+j)
             ff = ynew(j)
             gnew(j,i) = (ff-fo)/dx
           end do
        end do

        end subroutine evagrad
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine eva(mode,xnew,ynew)
        use dimGA
! mode=1 : usual evaluation of GA
! mode=2 : FD for gradient
! mode=3 : 1D search
        implicit none
        integer, intent(in) :: mode
        double precision, dimension(ndv),      intent(in)  :: xnew
        double precision, dimension(ndat-ndv), intent(out) :: ynew
        integer :: i
! likelihood
        double precision :: llfd,verr,pena
! Within Here, You have to Get All Objectives and Constraint [ynew]

        do i=1,ndv
          if(mode.eq.1)then
           if(xnew(i).lt.0.d0.or.xnew(i).gt.1.d0) &
           write(*,'(a,f12.7)')'*Outside of Design Space (1)??',xnew(i)
          else
           if(xnew(i).lt.-0.001d0.or.xnew(i).gt.1.001d0) &
           write(*,'(a,f12.7)')'*Outside of Design Space (2)??',xnew(i)
          end if
        end do
        if(mode.le.0.or.mode.ge.4)stop'unknown mode in eva'
        ict_eva(mode) = ict_eva(mode) + 1

        ynew = 0.d0
        if(neva.lt.0)then
          call testfunction(xnew,ynew)

        else if(neva.ge.1.and.neva.le.6)then
          call DV2DESV1(neva,ntmp,ttmp,ndv,bd_t,xnew,dv_t,dv_p)
          call likelihood(dv_t,dv_p,llfd,verr)
          ynew(1) = llfd
          if(verr.lt.1.d-3)then ! to be same with the value in likelihood
            ynew(2) = 0.d0
          else
            ynew(2) = verr
          end if

        else if(neva.ge.11.and.neva.le.20)then
          if(neva.eq.11.or.neva.eq.12)then
            dv_t(:,1) = xnew(:)
          else if(neva.eq.13.or.neva.eq.14)then
            call DV2DESV2(xnew,vrange,dv_t)
          else if(neva.eq.15)then
            call DV2DESV2(xnew,vrange,dv_t)
          end if
          call meta_all(neva,nobj,dv_t,ynew)
!         write(*,'(3f8.3,99e12.3)')vrange,(dv_t(i),i=1,2),(ynew(i),i=1,nobj+1)

        else
          stop'unknown neva in eva'
        end if

        end subroutine eva
