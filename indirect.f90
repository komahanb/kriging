        subroutine indirect(ifunc)
        use dimKrig
        implicit none
        integer, intent(in) :: ifunc
        integer :: i,j,k
        integer :: ict,nzero,nadd
        integer :: nput1,nput2,rnews,rnews2
        integer, dimension(-1:5) :: ictgh
        integer, allocatable, dimension(:)     :: nlap
        integer, allocatable, dimension(:,:,:) :: olap
        double precision :: dd,factor,facmin
        double precision, dimension(     ndim) :: xchc
        double precision, dimension(   1,   1) :: d1y,d2y
        double precision, dimension(   1,ndim) :: vdx
        double precision, dimension(ndim,   1) :: vg
        double precision, dimension(ndim,ndim) :: vh,xgrd
        double precision, allocatable, dimension(:) :: x,deg
        double precision, allocatable, dimension(:,:) :: xlhs,dlhs
        double precision, allocatable, dimension(:) :: tfun
        double precision, allocatable, dimension(:,:) :: tsampl
        character(len=6), allocatable, dimension(:)   :: tinf

        if(     ngput.eq.-2)then
           nput1 = ndim
        else if(ngput.eq.-1)then
           nput1 = 1
        else
           stop'unknown ngput'
        end if
        if(     nhput.eq.-4)then
           nput2 = nhes + ndim
        else if(nhput.eq.-3)then
           nput2 = ndim*2
        else if(nhput.eq.-2)then
           nput2 = ndim
        else if(nhput.eq.-1)then
           nput2 = 1
        else if(nhput.ge.1)then
           nput2 = nhput
        else
           stop'unknown nhput'
        end if

! estimation
        ictgh = 0
        do i=1,rsample
           if(     inf(i)(3:6).eq.'F   ')then
             ictgh(0) = ictgh(0) + 1
           else if(inf(i)(3:6).eq.'FG  ')then
             ictgh(1) = ictgh(1) + 1
           else if(inf(i)(3:6).eq.'FH  ')then
             ictgh(2) = ictgh(2) + 1
           else if(inf(i)(3:6).eq.'FGH ')then
             ictgh(3) = ictgh(3) + 1
           else
             stop'Hessvec for indirect?'
           end if
           if(inf(i)(3:6).eq.'F   ')tdx(i) = 0.d0
           if(itrust.eq.0)then
             if(inf(i)(3:6).eq.'FG  ')then
               if(tdx(i).ne.tdxinit(1))then !stop'tdx-tdxinit1'
                 tdx(i) = tdxinit(1)
               end if
             else if(inf(i).eq.'FGH ')then
               if(tdx(i).ne.tdxinit(2))stop'tdx-tdxinit2'
             end if
           end if
        end do
        if(ictgh(2).ne.0)stop'indirect for only Hessian?'
        rnews = ictgh(0) * (      1) &
              + ictgh(1) * (nput1+1) &
              + ictgh(3) * (nput2+1)

        if(id_proc.eq.0.and.Cmode(:4).eq.'Make')then
          if(itrust.eq.0)then
            write(*,'(1x,a,2f10.6)')'>> dx is fixed to ',(tdxinit(i),i=1,2)
            write(*,'(6x,a,i5)') & 
            '>> Max Possible Number of Sample Points = ',rnews+1
          else if(itrust.eq.1)then
            write(*,'(1x,a)')'>> Variable dx is used'
            write(*,'(6x,a,i5)') & 
            '>> Max Possible Number of Sample Points = ',rnews+1
          else
            stop'unknown itrust'
          end if
        end if

! overlap
        allocate(olap(rsample,rsample,2))
        allocate(nlap(rsample))
        call check_overlap(nlap,olap)

! arrangement for additional points (gradient)
        xgrd = 0.d0
        if(ngput.eq.-1)then
          do i=1,ndim
            xgrd(1,i) = 1.d0/dsqrt(dble(ndim))
          end do
        else if(ngput.eq.-2)then
          do i=1,ndim
            xgrd(i,i) = 1.d0
          end do
        else
          stop'unknown ngput'
        end if

! arrangement for additional points (hessian)
        allocate(xlhs(nput2,ndim))
        if(nhput.eq.-1)then
          if(nput2.ne.   1)stop'nput2 for -1'
          do i=1,ndim
            xlhs(1,i) = 1.d0/dsqrt(dble(ndim))
          end do
        else if(nhput.eq.-2)then
          if(nput2.ne.ndim)stop'nput2 for -2'
          do i=1,ndim
            xlhs(i,i) = 1.d0
          end do
        else if(nhput.eq.-3)then
          if(nput2.ne.ndim*2)stop'nput2 for -2'
          xlhs = 0.d0
          ict  = 0
          do i=1,ndim
            do j=-1,1,2
               ict = ict + 1
               xlhs(ict,i) = dble(j)
            end do
          end do
          if(ict.ne.ndim*2)stop'ict.ne.ndimx2'
        else
          allocate (dlhs(nput2,ndim-1))
          allocate (deg(ndim-1))
          allocate (x(ndim))
          call lhs_polar(ndim,nput2,dlhs)
          do i=1,nput2
             deg(:) = dlhs(i,:)
             call calc_polar(ndim,1.d0,deg,x,dd)
             if(dabs(dd-1.d0).ge.0.001d0)stop'dd in lhs_polar'
             xlhs(i,:) = x(:)
             do j=1,i-1
                dd = 0.d0
                do k=1,ndim
                   dd = dd + (xlhs(i,k)-xlhs(j,k))**2
                end do
                dd = dsqrt(dd)
                if(ndim.ne.1)then
                  if(dd.le.1.d-10) &
                  stop'additional pts are so close, change ilhs'
                  if(dd.le.1.d0/dble(nput2)) &
                  write(*,'(a,2i5,f15.8)')'*Additionals are so close?',i,j,dd
                end if
             end do
          end do
          deallocate( dlhs,deg,x )
        end if

! increase the size of dimension
        allocate(tsampl(rsample,ndim))
        allocate(tfun(rsample))
        allocate(tinf(rsample))
        do i=1,rsample
          tsampl(i,:)  = sampl(i,:)
          tfun(i)      =  fun(i)
          tinf(i)      = inf(i)
        end do
        deallocate(sampl,fun,nparent,dxadd,inf)
        allocate(nparent(rnews+1)) ! +1 is for special treatment
        allocate(dxadd(rnews+1))
        allocate(sampl(rnews+1,ndim))
        allocate(fun(rnews+1))
        allocate(inf(rnews+1))
        inf = '      '
        do i=1,rsample
          sampl(i,:) = tsampl(i,:)
          fun(i)     = tfun(i)
          inf(i)     = tinf(i)
        end do
        deallocate(tsampl,tfun,tinf)

! make additional points
        ict     = rsample
        nzero   = 0
        nparent = 0
        dxadd   = 0.d0
        do 100 i=1,rsample
          if(tdx(i).le.0.d0)then
            if(inf(i)(3:6).eq.'FG  ')nzero = nzero + nput1
            if(inf(i)(3:6).eq.'FGH ')nzero = nzero + nput2
            go to 100
          end if
          if(inf(i)(3:6).eq.'F   ')then
             go to 100
          else if(inf(i)(3:6).eq.'FG  ')then
             vg(:,1) = gfun(i,:)
             vh      = 0.d0
             nadd    = nput1
          else if(inf(i)(3:6).eq.'FH  ')then
             stop'only Hess for indirect?'
          else if(inf(i)(3:6).eq.'FGH ')then
             vg(:,1) = gfun(i,:)
             vh(:,:) = hfun(i,:,:)
             nadd    = nput2
          else if(inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv')then
             stop'Hessian Vector for indirect?'
          else
             stop'unknown info in indirect'
          end if

          do 150 j=1,nadd
            if(inf(i)(3:6).eq.'FG  ')then
              vdx(1,:) = xgrd(j,:) * tdx(i)
            else if(inf(i)(3:6).eq.'FGH ')then
              vdx(1,:) = xlhs(j,:) * tdx(i)
            end if
            xchc(:) = sampl(i,:) + vdx(1,:)
            facmin  = 1.d0
            do 160 k=1,ndim
              factor = 1.d0
              if(xchc(k).lt.0.d0)then
                factor = (0.d0-sampl(i,k))/(xchc(k)-sampl(i,k))
              else if(xchc(k).gt.1.d0)then
                factor = (1.d0-sampl(i,k))/(xchc(k)-sampl(i,k))
              end if
              facmin = min(facmin,factor)
160         continue
            if(facmin*tdx(i).le.tdxmin.or.facmin.le.0.25d0)then
              nzero = nzero + 1
              go to 150
            end if

            if(inf(i)(3:6).eq.'FG  ')then
              vdx(1,:) = xgrd(j,:) * tdx(i) * facmin
            else if(inf(i)(3:6).eq.'FGH ')then
              vdx(1,:) = xlhs(j,:) * tdx(i) * facmin
            end if
            d1y = matmul(vdx,vg)
            d2y = matmul( matmul(vdx,vh),transpose(vdx) )
            ict = ict + 1
            sampl(ict,:) = sampl(i,:) + vdx(1,:)
            fun(ict)     = fun(i) + d1y(1,1) + 0.5d0*d2y(1,1)
            nparent(ict) = i
            dxadd(ict)   = tdx(i)*facmin
150       continue
100     continue ! real sample loop(i)
        if(ict.ne.rnews-nzero)then
          write(*,*)ict,rnews,nzero
          stop'rnews-nzero'
        end if

! elimination for overlap
        call reduce_additional_pts(rnews,nlap,olap,rnews2)
        rsample = rnews2
!       rsample = rnews

        deallocate(xlhs)
        deallocate(olap,nlap)

        call check_nugget
        end subroutine indirect
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine check_nugget!(ndebug,id_proc,rsample,nparent,dxadd,Cnug,Vnug)
        use dimKrig
        implicit none
!       integer, intent(in) :: ndebug,id_proc,rsample
!       integer, dimension(rsample), intent(in) :: nparent
!       double precision, dimension(rsample), intent(in) :: dxadd
!       character(len=1), intent(in) :: Cnug
!       double precision, intent(in) :: Vnug
        integer :: i,j
        double precision :: value,vmin,vmax

        vmin = 1.d0
        vmax = 0.d0
        do i=1,rsample
        do j=i,rsample
           call calc_nugget(rsample,nparent,dxadd,i,j,value)
           if(ndebug.eq.1)then
             write(*,'(11x,a,f10.5,2i5)')'>> Nugget = ',value,i,j
           end if
           vmin = min(vmin,value)
           vmax = max(vmax,value)
        end do
        end do
        if(id_proc.eq.0) &
        write(*,'(6x,a,2e15.5)')'>> Nugget_min/max = ',vmin,vmax

!       call stop_all
        end subroutine check_nugget
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine calc_nugget(rsample,nparent,dxadd,n1,n2,vnugget)
        implicit none
        integer, intent(in) :: rsample,n1,n2
        integer, dimension(rsample), intent(in) :: nparent
        double precision, dimension(rsample), intent(in) :: dxadd
!       character(len=1), intent(in) :: Cnug
!       double precision, intent(in) :: Vnug
        double precision, intent(out) :: vnugget
        integer :: np1,np2
        double precision :: v1,v2,dx1,dx2

        np1 = nparent(n1)
        np2 = nparent(n2)
        dx1 = dxadd(n1)
        dx2 = dxadd(n2)
        if(np1.eq.0.and.np2.eq.0)then ! both are real pts
           vnugget = 1.d0
        else 
           v1 = 0.d0
           v2 = 0.d0
           if(np1.ne.0)call calc_nugget_cont(np1,dx1,v1)
           if(np2.ne.0)call calc_nugget_cont(np2,dx2,v2)
           vnugget = (1.d0-v1) * (1.d0-v2)
        end if

        end subroutine calc_nugget
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine calc_nugget_cont(np,dx,v)
        use dimKrig ! only for Cnug,Vnug
        implicit none
!       character(len=1), intent(in) :: Cnug
!       double precision, intent(in) :: Vnug
        integer, intent(in) :: np
        double precision, intent(in) :: dx
        double precision, intent(out) :: v

        if(Cnug.eq.'N')then
           v = Vnug
        else if(Cnug.eq.'C')then !.or.Cnug.eq.'M')then
           if(dx.le.0.d0)stop'dx=0 in calc_nugget_cont'
           v = Vnug * (dx**2)
        else
           stop'unknown Cnug'
        end if

        end subroutine calc_nugget_cont
