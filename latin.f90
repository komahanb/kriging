        subroutine Latin_Hypercube(ndim,nsam,bound,sample)
        implicit none
        integer, intent(in) :: ndim,nsam
        double precision, dimension(2,ndim), intent(in) :: bound
        double precision, dimension(ndim,nsam), intent(out) :: sample

        integer :: i,j,k,l,key
        double precision :: ran,dv01
        double precision, dimension(0:nsam) :: split
        integer, dimension(ndim,nsam) :: iflag_split

        if(nsam.eq.0)return
        if(nsam.lt.0.or.ndim.le.0)stop'LHS dimension/nsample'
        sample = 0.d0

! Set Cartesian
        do 100 i=0,nsam
           split(i) = dble(i)/dble(nsam)
100     continue

! Arrangement
        do 200 i=1,nsam
           do 210 j=1,ndim
250           continue
              call random_number(ran)
              do 220 k=1,nsam
                 if(ran.ge.split(k-1).and.ran.le.split(k))then
                   do 230 l=1,i-1
                      if(iflag_split(j,l).eq.k)go to 250 ! already used, try again
230                continue    ! previous sample loop
                   iflag_split(j,i) = k
                   go to 210
                 end if
220           continue         ! split loop
              stop'out of order at split loop'
210        continue            ! dimension loop
200     continue               ! sample loop
        
! Check
        do 300 j=1,ndim
          do 310 i=1,nsam
             if(iflag_split(j,i).le.0.or.iflag_split(j,i).gt.nsam) &
             stop'iflag_split'
             do 320 k=i+1,nsam
                if(iflag_split(j,i).eq.iflag_split(j,k)) &
                stop'overlap of iflag_split'
320          continue
310       continue
300     continue

! Generation
        do 400 i=1,nsam
          do 410 j=1,ndim
             key  = iflag_split(j,i)
             call random_number(ran)      ! random set
             ran = 0.01d0 + 0.99d0*ran
!            ran  = 0.5d0                 ! center of cartesian
             dv01 =  split(key-1) + ran*(split(key)-split(key-1))
             sample(j,i) = bound(1,j) + (bound(2,j)-bound(1,j))*dv01
410       continue
400     continue

        end subroutine Latin_Hypercube
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine calc_polar(ndv,rr,deg,x,dd)
      implicit none
      integer, intent(in) :: ndv
      double precision, intent(in)  :: rr
      double precision, dimension(ndv-1), intent(in) :: deg
      double precision, dimension(ndv),  intent(out) :: x
      double precision, intent(out) :: dd
      integer :: i,j

        x(1) = rr * cos(deg(1))
!       x(2) = rr * sin(deg(1)) * cos(deg(2))
        do i=2,ndv-1
           x(i) = rr
           do j=1,i-1
              x(i) = x(i) * sin(deg(j))
           end do
           x(i) = x(i) * cos(deg(i))
        end do
        x(ndv) = rr
        do j=1,ndv-1
          x(ndv) = x(ndv) * sin(deg(j))
        end do

        dd = 0.d0
        do i=1,ndv
           dd = dd + x(i)**2
        end do
        dd = dsqrt(dd)
      end subroutine calc_polar
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine lhs_polar(mdv,msample,DV)
      implicit none
      integer, intent(in) :: mdv,msample
      double precision, dimension(msample,mdv-1), intent(out) :: DV
      integer, dimension(msample,mdv-1)          :: isplit
      double precision, dimension(0:msample)     :: split
      double precision, dimension(0:msample,mdv) :: wsplit
      integer :: i,j,k,l,ilhs,flag,ndiv
      double precision :: x,ran1,vtarget,pow,pi
      double precision :: x1,x2,f1,f2,vdiv,vi
      double precision, allocatable, dimension(:,:) :: tmpreg
! # of sampling is  for mdv - 1
! for 1<i<mdv-2:  0<deg< pi
! for   i=mdv-1:  0<deg<2pi

      ilhs   = -130
      pi     = 4.d0*datan(1.d0)

      DV     = 0.d0
      split  = 0.d0
      wsplit = 0.d0
      isplit = 0

      do 100 i=0,msample
         split(i)  = dble(i) / dble(msample)                      ! [0:1]
!        vtarget   = 1.d0 - 2.d0*dble(i)/dble(msample)            ! [1:-1]
!        vtarget2    = vtarget**pow
!        wsplit(i,j) = (dacos(vtarget2))                          ! [0:pi]
100   continue

      ndiv   = max(msample*10,1001)
      allocate(tmpreg(ndiv,2))
      do 110 j=1,mdv-2
        pow    = dble(mdv-1-j)
        vi     = 0.d0
        tmpreg = 0.d0
        do 120 k=1,ndiv-1
           x1 = dble(k-1)/dble(ndiv-1) * pi
           x2 = dble(k  )/dble(ndiv-1) * pi
           f1 = (sin(x1))**pow
           f2 = (sin(x2))**pow
           vi = vi + (x2-x1)*(f1+f2)*0.5d0
           tmpreg(k+1,1) = x2
           tmpreg(k+1,2) = vi
!          write(1000+j,'(2f15.8)')x2,vi
120     continue
        vdiv = vi / dble(msample)
        wsplit(0,j) = 0.d0
        do 130 i=1,msample-1
          vtarget = vdiv * dble(i)
          do 140 k=2,ndiv
            if(tmpreg(k-1,2).le.vtarget.and.tmpreg(k,2).ge.vtarget)then
              wsplit(i,j) = 0.5d0*(tmpreg(k-1,1)+tmpreg(k,1))
            end if
140       continue
!         write(2000+j,'(2f15.8)')wsplit(i,j),vdiv*(dble(i)-0.5d0)
130     continue
        wsplit(msample,j) = pi                   ! [0:pi]
!check
        do 150 i=1,msample-1
           if(wsplit(i,j).ge.wsplit(i+1,j))then
              write(*,*)'ndv=',j
              write(*,*)'partion=',i,i+1
              write(*,*)wsplit(i,j),wsplit(i+1,j)
              stop'wsplit'
           end if
150     continue
110   continue
      deallocate(tmpreg)

      do 200 i=1,msample
        do 210 j=1,mdv-1
250      continue
         x = ran1(ilhs)                           ! [0:1]
!        call random_number(x)
         do 220 k=1,msample
            if(x.ge.split(k-1).and.x.le.split(k))then
              do 230 l=1,i-1
               if(isplit(l,j).eq.k)go to 250
230           continue
              isplit(i,j) = k
              go to 210
            end if
220      continue
         stop'out of order'
210     continue
200   continue
!check
      do 300 j=1,mdv-1
        do 310 i=1,msample
           if(isplit(i,j).le.0.or.isplit(i,j).gt.msample)stop'region'
           do 320 k=i+1,msample
              if(isplit(i,j).eq.isplit(k,j))stop'LHS'
320        continue
310     continue
300   continue
!check
      do 400 i=1,msample
        do 410 j=1,mdv-1
          flag = isplit(i,j)
!         x = ran1(ilhs)
          x = 0.5d0
          if(j.le.mdv-2)then        ! [0: pi]
            DV(i,j)  = wsplit(flag-1,j) + x*(wsplit(flag,j)-wsplit(flag-1,j))
          else if(j.eq.mdv-1)then   ! [0:2pi]
            DV(i,j)  =  split(flag-1)   + x*( split(flag)  - split(flag-1)  )
            DV(i,j)  = DV(i,j) * 2.d0*pi
          end if
410     continue
400   continue

      end subroutine lhs_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision function ran1(idum)
      implicit none
      integer, intent(inout) :: idum
      integer :: IA,IM,IQ,IR,NTAB,NDIV,iy
      double precision :: AM,EPS,RNMX
      parameter (IA=16807,IM=2147483647,IQ=127773,IR=2836, &
                 NTAB=32,NDIV=1+(IM-1)/NTAB)
!     parameter (AM=1.d0/dble(IM))
      parameter (EPS=1.2e-7)
      parameter (RNMX=1.d0-EPS)
!     parameter (AM=1.d0/dble(IM),EPS=1.2e-7,RNMX=1.d0-EPS)
      integer, dimension(NTAB) :: iv
      data iv /NTAB*0/,iy /0/
      integer :: j,k
!     write(*,*)'iy',iy
      AM = 1.d0 / dble(IM)

      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif

      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      end function ran1

