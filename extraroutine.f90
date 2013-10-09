        subroutine combination(n,m,l)
        implicit none
! l = nCm = (n!)/(m!)/((n-m)!)
        integer, intent(in)  :: n,m
        integer, intent(out) :: l
        integer :: i,i1,i2,i3,nbig,npet
        if(n.le.0.or.m.le.0)stop'n,m<1 in combination'
        if(n.lt.m)stop'n<m in combination'

!       i1 = 1
!       i2 = 1
!       i3 = 1
!       do i=1,n
!          i1 = i1 * i
!       end do
!       do i=1,m
!          i2 = i2 * i
!       end do
!       do i=1,n-m
!          i3 = i3 * i
!       end do
!       l = int( dble(i1)/dble(i2*i3) )

        if(m.ne.n-m)then
          nbig = max(m,n-m)
          npet = min(m,n-m)
        else
          nbig = m
          npet = n-m
        end if
! l = nCm = (n!)/(nbig!)/(npet!) = (n*(n-1)*...*(nbig+1))/(npet!)
        i1 = 1
        i2 = 1
        do i=nbig+1,n
          i1 = i1 * i
        end do
        do i=1,npet
          i2 = i2 * i
        end do
        l = int( dble(i1)/dble(i2) )

        end subroutine combination
subroutine knn(SC,sample,knnptr,ndim,nhs,NCP)
  ! Subroutine to find NCP closest neighbours from array sample to point SC
implicit none
integer :: ndim,nhs,NCP,j,k,knnptr(NCP)
real*8 ::SC(ndim),sample(nhs,ndim),norm2(nhs),sn

! Calculate all distances
do j=1,nhs
  norm2(j)=0.0
  do k=1,ndim
     norm2(j)=norm2(j)+(sample(j,k)-SC(k))**2
  end do
end do

! Find NCP closest neighbours

do k=1,NCP
  sn=10000000000.0
  do j=1,nhs
     if (norm2(j).lt.sn) then
        sn=norm2(j)
        knnptr(k)=j
     end if
  end do
  norm2(knnptr(k))=10000000000.0
end do

end subroutine knn
    subroutine mirtunableparams(fct,ndim,nhs,ncp,taylororder,NTOEX)
      implicit none
      integer,INTENT(IN)::fct,ndim,nhs
      INTEGER,INTENT(OUT)::NCP,NTOEX,TAYLORORDER
      
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !                 EXP
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      if (Fct.eq.4) then
         if (nhs.le.11)  then  
           NCP=nhs     
           Taylororder=3!nhs!2!2!INT(NHS/4)

        else if (nhs.gt.11 .and. nhs.le.15)  then  
           NCP=10
           Taylororder=5!ncp

        else  if (nhs.gt.15 .and. nhs.le.25)  then  
           NCP=15
           Taylororder=5!ncp
           
        else  if (nhs.gt.25 .and. nhs.le.35)  then  
           NCP=20
           Taylororder=5!ncp
           
        else if (nhs.gt.35 .and. nhs.le.100)  then  
           NCP=20!20+3*ndim
           Taylororder=5
       ! else if (nhs.gt.70 .and. nhs.le.100)  then  
       !   NCP=30!50+0.1*nhs
       !    Taylororder=17       
        else
           NCP=20!+5*ndim
           tAYLORORDER=5
        end if
        NTOEX=(30-ndim)*NCP
     end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!    
!                COS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
     if (Fct.eq.1) then
        if (nhs.le.11)  then  
           NCP=nhs
           Taylororder=ncp!2!2!INT(NHS/4)
        else if (nhs.gt.11 .and. nhs.le.15)  then  
           NCP=10
           Taylororder=5!ncp

        else if (nhs.gt.15 .and. nhs.le.25)  then  
           NCP=15
           Taylororder=5!ncp

        else if (nhs.gt.25 .and. nhs.le.35)  then  
           NCP=20
           Taylororder=5!ncp

        else if (nhs.gt.35 .and. nhs.le.70)  then  
           NCP=20!50+0.1*nhs
           Taylororder=5
        else if (nhs.gt.70 .and. nhs.le.100)  then  
           NCP=20!50+0.1*nhs
           Taylororder=5       
        else
           NCP=20
           tAYLORORDER=5
        end if
NTOEX=(30-ndim)*NCP
     end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!               RUNGE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!     
     if (fct.eq.2) then
        if (nhs.le.20)  then  
           NCP=nhs
           Taylororder=5!nhs!2!2!INT(NHS/4)

           !   else if (nhs.gt.11 .and. nhs.le.15)  then  
           !      NCP=nhs
           !      Taylororder=7!ncp
           !      
           !   else if (nhs.gt.15 .and. nhs.le.25)  then  
           !      NCP=nhs
           !      Taylororder=13!ncp
           !      
           !   else if (nhs.gt.25 .and. nhs.le.35)  then  
           !      NCP=nhs
           !      Taylororder=17!ncp
           !      
           !   else if (nhs.gt.35 .and. nhs.le.70)  then  
           !      NCP=30!50+0.1*nhs
           !      Taylororder=17
        else  if (nhs.gt.20 .and. nhs.le.100)  then  
           NCP=20!50+0.1*nhs
           Taylororder=5       
        else
           NCP=20
           tAYLORORDER=5
        end if
NTOEX=(30-ndim)*NCP
     end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!           ROSENBROCK    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
     if (Fct.eq.6) then
        if (nhs.le.11)  then  
           NCP=5
           Taylororder=ncp!2!2!INT(NHS/4)
        else if (nhs.gt.11 .and. nhs.le.15)  then  
           NCP=10
           Taylororder=5!ncp

        else if (nhs.gt.15 .and. nhs.le.25)  then  
           NCP=15
           Taylororder=5!ncp

        else if (nhs.gt.25 .and. nhs.le.35)  then  
           NCP=25
           Taylororder=5!ncp

        else if (nhs.gt.35 .and. nhs.le.70)  then  
           NCP=35!50+0.1*nhs
           Taylororder=5
        else if (nhs.gt.70 .and. nhs.le.100)  then  
           NCP=50!50+0.1*nhs
           Taylororder=5       
        else
           NCP=70
           tAYLORORDER=5
        end if
NTOEX=(30-ndim)*NCP
     end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                      CFD   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
     if (fct.eq.20) then
        if (nhs.le.10)  then  
           NCP=nhs
           Taylororder=nhs!nhs!2!2!INT(NHS/4)
        else if (nhs.gt.10 .and. nhs.le.15)  then  
           NCP=nhs
           Taylororder=10!ncp
           
        else  if (nhs.gt.15 .and. nhs.le.25)  then  
           NCP=nhs
           Taylororder=10!ncp
           
        else  if (nhs.gt.25 .and. nhs.le.35)  then  
           NCP=nhs
           Taylororder=7!ncp
           
        else if (nhs.gt.35 .and. nhs.le.45)  then  
           NCP=30
           Taylororder=7
           
        else if (nhs.gt.45 .and. nhs.le.100)  then  
           NCP=30!50+0.1*nhs
           Taylororder=7
        else
           NCP=35
           Taylororder=7
           
        end if
        NTOEX=(30-ndim)*NCP

     else


        if (nhs.le.11)  then  
           NCP=nhs
           Taylororder=ncp!2!2!INT(NHS/4)
        else if (nhs.gt.11 .and. nhs.le.15)  then  
           NCP=10
           Taylororder=5!ncp

        else if (nhs.gt.15 .and. nhs.le.25)  then  
           NCP=15
           Taylororder=5!ncp

        else if (nhs.gt.25 .and. nhs.le.35)  then  
           NCP=20
           Taylororder=5!ncp

        else if (nhs.gt.35 .and. nhs.le.70)  then  
           NCP=20!50+0.1*nhs
           Taylororder=5
        else if (nhs.gt.70 .and. nhs.le.100)  then  
           NCP=20!50+0.1*nhs
           Taylororder=5       
        else
           NCP=20
           tAYLORORDER=5
        end if
        NTOEX=(30-ndim)*NCP


     end if ! end of CFD 



   end subroutine mirtunableparams
  recursive subroutine recurcorn(dimact,ndim,sample,nhs,bounds)

    implicit none

    common/global/counter

    integer :: dimact,ndim,nhs,counter,j,counterin
    real*8 :: sample(ndim,nhs),bounds(2,ndim)

    counterin=counter
    do j=2,1,-1
       sample(dimact,counter)=bounds(j,dimact)
       if (dimact.gt.1) then
          sample(1:dimact-1,counter)=sample(1:dimact-1,counterin)
       end if
       if (dimact.lt.ndim) call recurcorn(dimact+1,ndim,sample,nhs,bounds)
       if (j.ne.1) counter=counter+1
    end do


  end subroutine recurcorn







subroutine latin_random ( dim_num, point_num, seed, x )

!*****************************************************************************80
!
!! LATIN_RANDOM returns points in a Latin Random square.
!
!  Discussion:
!
!    In each spatial dimension, there will be exactly one
!    point whose coordinate value lies between consecutive
!    values in the list:
!
!      ( 0, 1, 2, ..., point_num ) / point_num
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator,
!    needed if the portable D_UNIFORM_01 routine is being used.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) :: base = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) perm(point_num)
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) x(dim_num,point_num)
!
!  Pick DIM_NUM * POINT_NUM random numbers between 0 and 1.
!
!  For fast results, use the FORTRAN90 standard RANDOM_NUMBER routine.
!  For reproductible results, use the UNIFORM routine.
!
  if ( .false. ) then

    call random_number ( harvest = x(1:dim_num,1:point_num) )

  else

    do i = 1, dim_num
      do j = 1, point_num
        x(i,j) = r8_uniform_01 ( seed )
      end do
    end do

  end if
!
!  For spatial dimension I, 
!    pick a random permutation of 1 to POINT_NUM,
!    force the corresponding I-th components of X to lie in the
!    interval ( PERM(J)-1, PERM(J) ) / POINT_NUM.
!
  do i = 1, dim_num

    call perm_uniform ( point_num, base, seed, perm )

    do j = 1, point_num
      x(i,j) = ( real ( perm(j) - 1, kind = 8 ) + x(i,j) ) &
               / real ( point_num, kind = 8 )
    end do

  end do

  return
end subroutine latin_random



subroutine perm_uniform ( n, base, seed, p )

!*****************************************************************************80
!
!! PERM_UNIFORM selects a random permutation of N objects.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input, integer ( kind = 4 ) BASE, is 0 for a 0-based permutation and 1 for 
!    a 1-based permutation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) P(N), the permutation.  P(I) is the "new"
!    location of the object originally at I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  do i = 1, n
    p(i) = ( i - 1 ) + base
  end do

  do i = 1, n
    j = i4_uniform ( i, n, seed )
    k    = p(i)
    p(i) = p(j)
    p(j) = k
  end do

  return
end subroutine perm_uniform


function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer ( kind = 4 ) arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end function r4_uniform_01


function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer ( kind = 4 ) arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer ( kind = 4 ),
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function r8_uniform_01


function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end function i4_uniform


subroutine lhs_normal_dist(mdv,msample,xc,st,DV)
      implicit none
      integer, intent(in) :: mdv,msample
      double precision, dimension(mdv), intent(in) :: xc
      double precision, intent(in) :: st
      double precision, dimension(msample,mdv), intent(out) :: DV

      integer, dimension(msample,mdv)            :: isplit
      double precision, dimension(0:msample)     :: split
      double precision, dimension(0:msample,mdv) :: wsplit
      double precision, dimension(mdv,2) :: xbd
      integer :: i,j,k,l,ilhs,flag,ndiv
      double precision :: x,ran1,vtarget,eps
      double precision :: x1,x2,f1,f2,vdiv,vi
      double precision, allocatable, dimension(:,:) :: tmpreg

      ilhs   = -130
      DV     = 0.d0
      split  = 0.d0 ! equal spacing
      wsplit = 0.d0 ! variable spacing
      isplit = 0    ! location for (sample,dv)

      do 100 i=0,msample
         split(i)  = dble(i) / dble(msample)                      ! [0:1]
100   continue

      eps = 1.d-5 ! determine bound
      write(*,'(6x,a,e15.5)')'>> LHS for Normal Dist',eps
      do k=1,mdv
        call find_x(1,     eps,xc(k),st,xbd(k,1))
        call find_x(1,1.d0-eps,xc(k),st,xbd(k,2))
        write(*,'(11x,a,i3,a,2f12.7)')'>> Bound for',k,' = ',(xbd(k,j),j=1,2)
      end do

      ndiv   = max(msample*10,1001)
      allocate(tmpreg(ndiv,2))
      do 110 j=1,mdv
        vi     = 0.d0
        tmpreg = 0.d0
        tmpreg(1,1) = xbd(j,1)
        do 120 k=1,ndiv-1
           x1 = xbd(j,1) + (xbd(j,2)-xbd(j,1))*dble(k-1)/dble(ndiv-1)
           x2 = xbd(j,1) + (xbd(j,2)-xbd(j,1))*dble(k  )/dble(ndiv-1)
           call normal_dist(x1,xc(j),st,f1)
           call normal_dist(x2,xc(j),st,f2)
           vi = vi + (x2-x1)*(f1+f2)*0.5d0
           tmpreg(k+1,1) = x2
           tmpreg(k+1,2) = vi
!          write(1000+j,'(2f15.8)')x2,vi
120     continue
        vdiv = vi / dble(msample)
        wsplit(0,j) = xbd(j,1)
        do 130 i=1,msample-1
          vtarget = vdiv * dble(i)
          do 140 k=2,ndiv
            if(tmpreg(k-1,2).le.vtarget.and.tmpreg(k,2).ge.vtarget)then
              wsplit(i,j) = 0.5d0*(tmpreg(k-1,1)+tmpreg(k,1))
            end if
140       continue
!         write(2000+j,'(2f15.8)')wsplit(i,j),vdiv*(dble(i)-0.5d0)
130     continue
        wsplit(msample,j) = xbd(j,2)
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
        do 210 j=1,mdv
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
      do 300 j=1,mdv
        do 310 i=1,msample
           if(isplit(i,j).le.0.or.isplit(i,j).gt.msample)stop'region'
           do 320 k=i+1,msample
              if(isplit(i,j).eq.isplit(k,j))stop'LHS'
320        continue
310     continue
300   continue
!check
      do 400 i=1,msample
        do 410 j=1,mdv
          flag = isplit(i,j)
!         x = ran1(ilhs)
          x = 0.5d0
          DV(i,j)  = wsplit(flag-1,j) + x*(wsplit(flag,j)-wsplit(flag-1,j))
410     continue
400   continue

      end subroutine lhs_normal_dist


subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer   ( kind = 4 ) seed
  real      ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 )  today
  integer   ( kind = 4 ) values(8)
  character ( len = 5 )  zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) /  11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) /  30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) /  23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp                                    /   6.0D+00

  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( huge ( 1 ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer ( kind = 4 ).
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == huge ( 1 ) ) then
    seed = seed - 1
  end if

  return
end subroutine get_seed

subroutine i_to_s(intval,s)

  integer idig,intval,ipos,ival,i
  character ( len = * ) s

  s = ' '

  ival = intval

  !  Working from right to left, strip off the digits of the integer
  !  and place them into S(1:len ( s )).
  !
  ipos = len(s) 

  do while ( ival /= 0 )

     idig = mod( ival, 10 )
     ival = ival / 10
     s(ipos:ipos) = char(idig + 48 )
     ipos = ipos - 1

  end do
  !
  !  Fill the empties with zeroes.
  !
  do i = 1, ipos
     s(i:i) = '0'
  end do

  return
end subroutine i_to_s

