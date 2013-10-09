      subroutine LUDecomposition(mode,ideb,isize,ALU,Rinv,det,detlog,verr)
! mode=0 Try to LU Decomp Only
! mode=1 Try to Cholesky_Decomp Firstly, then try to LU Decomp if failure
!     elimination of allocations
      implicit none
      integer, intent(in) :: mode,ideb,isize
      double precision, dimension(0:isize-1,0:isize-1), intent(in)  :: ALU !Rij
      double precision, dimension(isize,isize), intent(out) :: Rinv
      double precision, intent(out) :: det,detlog,verr
      double precision, dimension(0:isize-1,0:isize-1) :: DL,DU,TMP
      integer :: i,j,IDIM,istat,ifail

      IDIM = isize
!     allocate(ALU(0:IDIM-1,0:IDIM-1),stat=istat)
!     if(istat.ne.0)stop'allocation miss'
!     allocate(DL(0:IDIM-1,0:IDIM-1),stat=istat)
!     if(istat.ne.0)stop'allocation miss'
!     allocate(DU(0:IDIM-1,0:IDIM-1),stat=istat)
!     if(istat.ne.0)stop'allocation miss'
!     allocate(TMP(0:IDIM-1,0:IDIM-1),stat=istat)
!     if(istat.ne.0)stop'allocation miss'

!     do i=0,IDIM-1
!       do j=0,IDIM-1
!         ALU(i,j) = Rij(i+1,j+1)
!       end do
!     end do

      if(ideb.eq.2)then
         write(*,*)'------ Write Matrix A'
         call WRITE_A(IDIM,ALU)
      end if

! beacuse all diagonals are 1 for correlation matrix
!     call PIBOT(IDIM,ALU)
!     if(ideb.eq.2)then
!        write(*,*)'------ Write Matrix A after Pibot'
!        call WRITE_A(IDIM,ALU)
!     end if

      if(mode.eq.1)then
        ifail = 1
        call CB_DECOMP(ideb,IDIM,ALU,Rinv,DL,verr,ifail)
      else
        ifail = 1
      end if

      if(ifail.eq.0)then
      else
        call LU_DECOMP(ideb,IDIM,ALU,DL,DU,TMP,verr)
        if(ideb.ge.1) &
        write(*,'(a,e15.5)') &
        '           >> Error in the LU Decomposition = ',verr
        if(ideb.eq.2)then
         write(*,*)'------ Write Matrix L*U'
         call WRITE_A(IDIM,TMP)
         write(*,*)'------ Write Matrix L'
         call WRITE_A(IDIM,DL)
         write(*,*)'------ Write Matrix U'
         call WRITE_A(IDIM,DU)
        end if
        call LU_INVERSION(ideb,IDIM,ALU,DL,DU,Rinv,verr)
      end if

      if(ideb.ge.1) &
      write(*,'(a,e15.5)') &
      '           >> Error in the Inversion of A   = ',verr

      det = 1.d0
      detlog = 0.d0
      do i=0,IDIM-1
        if(ifail.eq.0)then
          det    = det * (DL(i,i)**2)
          detlog = detlog + 2.d0*dlog(dabs(DL(i,i)))
        else
          det    = det * DU(i,i)
          detlog = detlog + dlog(dabs(DU(i,i)))
        end if
!       do j=0,IDIM-1
!         Rinv(i+1,j+1) = TMP(i,j)
!       end do
      end do
      if(ideb.ge.1) &
      write(*,'(a,2e15.5)')'           >> det,detlog = ',det,detlog

      if(ideb.eq.2)then
         write(*,*)'------ Write Matrix Ainv'
         call WRITE_A(IDIM,TMP)
         DL = matmul(ALU,TMP)
         write(*,*)'------ Write Matrix A*Ainv'
         call WRITE_A(IDIM,DL)
      end if

!     deallocate(ALU)
!     deallocate(DL)
!     deallocate(DU)
!     deallocate(TMP)
      end subroutine LUDecomposition
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine LU_INVERSION(ideb,IDIM,ALU,DL,DU,TMP,verr)
      implicit none
      integer, intent(in) :: IDIM,ideb
      double precision, dimension(0:IDIM-1,0:IDIM-1), intent(in) &
      :: ALU,DL,DU
      double precision, dimension(0:IDIM-1,0:IDIM-1), intent(out) &
      :: TMP
      double precision, intent(out) :: verr
      integer :: i,j,k,l,istat
      double precision :: vsum
      double precision, dimension(0:IDIM-1,0:IDIM-1) :: DLinv,DUinv

!     allocate(DLinv(0:IDIM-1,0:IDIM-1),stat=istat)
!     if(istat.ne.0)stop'allocation miss'
!     allocate(DUinv(0:IDIM-1,0:IDIM-1),stat=istat)
!     if(istat.ne.0)stop'allocation miss'
      DLinv = 0.d0
      DUinv = 0.d0

      do i=0,IDIM-1
        do j=0,i
           if(i.eq.j)then
              if(DL(i,j).eq.0.d0)stop'DL(i,j)'
              DLinv(i,j) = 1.d0 / DL(i,j)
           else
              if(DL(i,i).eq.0.d0)stop'DL(i,i)'
              vsum = 0.d0
              do k= 0,i-1
                 vsum = vsum - DL(i,k)*DLinv(k,j)
              end do
              DLinv(i,j) = vsum / DL(i,i)
           end if
        end do
      end do

      do i=IDIM-1,0,-1
        do j=IDIM-1,i,-1
           if(i.eq.j)then
              if(DU(i,j).eq.0.d0)then
                 write(*,*)'DU(i,j) = 0.0, IDIM = ',IDIM
                 if(IDIM.le.10)then
                   do k=0,IDIM-1
                    write(*,'(999f8.2)')(ALU(k,l),l=0,IDIM-1)
                   end do
                 end if
                 stop
              end if
              DUinv(i,j) = 1.d0 / DU(i,j)
           else
              if(DU(i,i).eq.0.d0)then
                 write(*,*)'DU(i,j) = 0.0, IDIM = ',IDIM
                 if(IDIM.le.10)then
                   do k=0,IDIM-1
                    write(*,'(999f8.2)')(ALU(k,l),l=0,IDIM-1)
                   end do
                 end if
                 stop
              end if
              vsum = 0.d0
              do k= IDIM-1,i+1,-1
                 vsum = vsum - DU(i,k)*DUinv(k,j)
              end do
              DUinv(i,j) = vsum / DU(i,i)
           end if
        end do
      end do


      if(ideb.eq.1)then
        TMP  = matmul(DL,DLinv)
        call calc_error_inversion(IDIM,TMP,verr)
        write(*,'(a,e15.5)') &
        '           >> Error in the Inversion of L   = ',verr
        TMP  = matmul(DU,DUinv)
        call calc_error_inversion(IDIM,TMP,verr)
        write(*,'(a,e15.5)') &
        '           >> Error in the Inversion of U   = ',verr
      end if

      TMP   = matmul(DUinv,DLinv)
      DLinv = matmul(ALU,TMP)     ! Rij*Rinv
      call calc_error_inversion(IDIM,DLinv,verr)

!     deallocate (DLinv,DUinv)
      end subroutine LU_INVERSION
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine LU_DECOMP(ideb,IDIM,ALU,DL,DU,TMP,verr)
      implicit none
      integer, intent(in) :: IDIM,ideb
      double precision, dimension(0:IDIM-1,0:IDIM-1), intent(in) :: ALU
      double precision, dimension(0:IDIM-1,0:IDIM-1), intent(out) :: DL,DU,TMP
      double precision, intent(out) :: verr
      integer :: i,j,k

      DU = 0.d0
      DL = 0.d0
      do i=0,IDIM-1
         DL(i,i) = 1.d0
      end do

      do i = 0,IDIM-1
         do j = i,IDIM-1
            DU(i,j) = ALU(i,j)
            do k = 0,i-1,1
               DU(i,j) = DU(i,j) - DU(k,j)*DL(i,k)
            end do
         end do
         do j = i+1,IDIM-1
            DL(j,i) = ALU(j,i)
            do k = 0,i-1,1
               DL(j,i) = DL(j,i) - DU(k,i)*DL(j,k)
            end do
            DL(j,i) = DL(j,i)/DU(i,i)
         end do
      end do

      if(ideb.ge.1)then
        TMP  = matmul(DL,DU)
        verr = 0.d0
        do i=0,IDIM-1
        do j=0,IDIM-1
           verr = verr + (ALU(i,j)-TMP(i,j))**2
        end do
        end do
        verr = dsqrt(verr)
      end if

      end subroutine LU_DECOMP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine PIBOT(IDIM,ALU)
      implicit none
      integer, intent(in) :: IDIM
      double precision, dimension(0:IDIM-1,0:IDIM-1), intent(inout) :: ALU
      integer :: ipibot,j,k,jref
      double precision :: Vmax,V,swap

99    continue
      ipibot = 0
      do 100 j=0,IDIM-1
           if(ALU(j,j).ne.0.d0)go to 100
           ipibot = ipibot + 1
           Vmax = 0.d0
           jref = -1
           do 101 k=0,IDIM-1
              if(ALU(j,k).eq.0.d0)go to 101
              V = abs(ALU(k,j))
              if(V.gt.Vmax)then
                 Vmax = V
                 jref = k
              end if
101        continue
           if(jref.eq.-1)then
              stop'no candidate in pibot'
           end if
           do 102 k = 0,IDIM-1
              swap        = ALU(jref,k)
              ALU(jref,k) = ALU(j,k)
              ALU(j   ,k) = swap
102        continue
!          write(*,'(a10,i3,a3,i3)')'piboting',j,' & ',jref
100   continue
      if(ipibot.ne.0)then
         stop'not considered yet'
         go to 99
      end if

      do 103 j=0,IDIM-1
         if(ALU(j,j).eq.0.d0)then
           write(*,*)'PIBOTING MISSED...',j,ALU(j,j)
           stop
         end if
103   continue

      end subroutine PIBOT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine WRITE_A(IDIM,A)
      implicit none
      integer, intent(in) :: IDIM
      double precision, dimension(0:IDIM-1,0:IDIM-1), intent(in) :: A
      integer :: i,j
      do i=0,IDIM-1
         write(*,'(i3,a,99e10.2)')i+1,' : ',(A(i,j),j=0,IDIM-1)
      end do
      end subroutine WRITE_A
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine calc_error_inversion(n,T,verr)
      implicit none
      integer, intent(in) :: n
      double precision, dimension(n,n), intent(in)  :: T
      double precision, intent(out) :: verr
      integer :: i,j

      verr = 0.d0
      do i=1,n
      do j=1,n
         if(i.eq.j)then
           verr = verr + (T(i,j)-1.d0)**2
         else
           verr = verr + (T(i,j)-0.d0)**2
         end if 
      end do
      end do
      verr = dsqrt(verr)

      end subroutine calc_error_inversion
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine CB_DECOMP(ideb,IDIM,ALU,TMP,DL,verr,ifail) 
      implicit none
      integer, intent(in) :: IDIM,ideb
      double precision, dimension(0:IDIM-1,0:IDIM-1), intent(in)  :: ALU
      double precision, dimension(0:IDIM-1,0:IDIM-1), intent(out) :: DL,TMP
      double precision, intent(out) :: verr
      integer, intent(out) :: ifail
      integer :: i,j

      call choldc(ideb,IDIM,ALU,TMP,DL,ifail,verr)
      if(ifail.ne.0)return

      end subroutine CB_DECOMP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine choldc(ideb,n,a,ainv,l,ifail,verr)
      implicit none
      integer, intent(in) :: ideb,n
      double precision, dimension(n,n), intent(in)  :: a
      double precision, dimension(n,n), intent(out) :: ainv,l
      integer, intent(out) :: ifail
      double precision, intent(out) :: verr

      double precision, dimension(n,n) :: linv,tmp
      integer :: i,j,k
      double precision :: sum
      double precision, dimension(n) :: p

      ifail = 0
      l = 0.d0
!     do i=1,n
!     do j=1,i
!        if(a(i,j).lt.0.d0)  stop'Aij < 0  in Cholesky Decomposition'
!        if(a(i,j).ne.a(j,i))stop'Aij!=Aji in Cholesky Decomposition'
!     end do
!     end do

      do 13 i=1,n
        do 12 j=i,n
            sum=a(i,j)
            do 11 k=j-1,1,-1
              sum=sum-l(i,k)*l(j,k)
11          continue
            if(i.eq.j)then
              if(sum.le.0.d0)then
                 write(*,'(11x,a,2i5,f15.8)') &
                 '>> cholesky-banachiewicz failed',i,j,sum
                 ifail = 1
                 return
              end if
              l(i,i)=sqrt(sum)
            else
              l(j,i)=sum/l(i,i)
            endif
!           write(*,'(1x,a,2i5,f15.8)')'>> i,j,l = ',j,i,l(j,i)
12      continue
13    continue
!     write(*,'(6x,a)')'>> cholesky-banachiewicz finish'

      linv = l
      do 23 i=1,n
        linv(i,i)=1.d0/l(i,i)
        do 22 j=i+1,n
          sum=0.d0
          do 21 k=i,j-1
            sum=sum-linv(j,k)*linv(k,i)
21        continue
          linv(j,i)=sum/l(j,j)
22      continue
23    continue

      if(ideb.eq.1)then
        tmp  = matmul(l,transpose(l))
        verr = 0.d0
        do i=1,n
        do j=1,n
           verr = verr + (a(i,j)-tmp(i,j))**2
        end do
        end do
        verr = dsqrt(verr)
        write(*,'(11x,a,e15.5)')'>> Error in Cholesky Decomp      = ',verr
        tmp = matmul(l,linv)
        call calc_error_inversion(n,tmp,verr)
        write(*,'(11x,a,e15.5)')'>> Error in L*Linv               = ',verr
        tmp = matmul(transpose(l),transpose(linv))
        call calc_error_inversion(n,tmp,verr)
        write(*,'(11x,a,e15.5)')'>> Error in Lt*Linvt             = ',verr
      end if

      ainv = matmul(transpose(linv),linv)
      tmp  = matmul(a,ainv)
      call calc_error_inversion(n,tmp,verr)
!     write(*,'(1x,a,e15.5)')'>> Error in A*Ainv   = ',verr

      end subroutine choldc
