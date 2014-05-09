!% Library to do distance operations on input array
!% Input:  Supply a multidimensional array X(ndim,npts)
!% Output: (1) Get avgdist, mindist, maxdist
!%         (2) Find distance between two points
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine findavg(npts,ndim,x,avgdist)
  implicit none

  integer :: npts,ndim
  real*8  :: x(ndim,npts)
  real*8  :: avgdist,dist(npts*(npts-1)/2)
  integer :: cnt,i,j

  cnt=0
  do i=1,npts
     do j=i+1,npts
        cnt=cnt+1
        call finddist(ndim,x(1:ndim,i),x(1:ndim,j),dist(cnt))
     end do
  end do

  avgdist=sum(dist)/dble(npts*(npts-1)/2)

  return
end subroutine findavg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine findmax(npts,ndim,x,maxdist)
  implicit none

  integer :: npts,ndim
  real*8  :: x(ndim,npts)
  real*8  :: maxdist,dist(npts*(npts-1)/2)
  integer :: cnt,i,j

  cnt=0
  do i=1,npts
     do j=i+1,npts
        cnt=cnt+1
        call finddist(ndim,x(1:ndim,i),x(1:ndim,j),dist(cnt))
     end do
  end do

  maxdist=maxval(dist)

  return
end subroutine findmax
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine findmin(npts,ndim,x,mindist)
  implicit none

  integer :: npts,ndim
  real*8  :: x(ndim,npts)
  real*8  :: mindist,dist(npts*(npts-1)/2)
  integer :: cnt,i,j

  cnt=0
  do i=1,npts
     do j=i+1,npts
        cnt=cnt+1
        call finddist(ndim,x(1:ndim,i),x(1:ndim,j),dist(cnt))
     end do
  end do

  mindist=minval(dist)

  return
end subroutine findmin

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine finddist(n,x1,x2,dist)

  implicit none

  integer,intent(in) :: n
  real*8,intent(in)  :: x1(n),x2(n)
  real*8,intent(out) :: dist
  integer::i

  dist=0.0d0
  do i=1,n
     dist=dist+(x1(i)-x2(i))**2
  end do
  dist=sqrt(dist)

  return
end subroutine finddist

!++++++++++++++++++++++++++++++++++++++++++++++++++++++
