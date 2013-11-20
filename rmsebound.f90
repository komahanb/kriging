
!++++++++++++++++++++++++++++++++++++++
subroutine find_mean(array,length,meanout)
implicit none

! Input and output variables
integer,intent(in):: length
real*8,intent(in) :: array(length)
real*8,intent(out):: meanout

! Internal data types
integer:: i
real*8::mean

mean=0.0d0

do i =1,length
mean=mean+array(i)
end do

mean = mean/dble(length)

meanout=mean
return
end subroutine find_mean
!++++++++++++++++++++++++++++++++++++++

subroutine find_std(array,length,stdout)
implicit none

! Input and output variables
integer,intent(in):: length
real*8,intent(in) :: array(length)
real*8,intent(out):: stdout

! Internal data types
integer:: i
real*8::std
real*8::mean


if (length.eq.1) then
   std=0.d0
   return
end if

mean=0.0
call find_mean(array,length,mean)

std=0.0d0
do i =1,length
std=std+(array(i)-mean)**2
end do

std = dsqrt(std/dble(length-1))

stdout=std

return
end subroutine find_std
!++++++++++++++++++++++++++++++++++++++
!!$program test
!!$  implicit none
!!$
!!$
!!$  real*8::x(5),mean,std,bnd(3)
!!$  integer::ks
!!$
!!$  x=(/0.5841,0.1078,0.9063,0.8797,0.8178/)
!!$!x=1.0
!!$
!!$!  call find_mean(x,5,mean)
!!$!  call find_std(x,5,std)
!!$!  ks=1
!!$  call make_bound(5,x,1.0,bnd)
!!$!  print*,bnd
!!$
!!$end program test



subroutine make_bound(length,array,ks,bnd)
implicit none

! Makes Error bar like upper and lower bound when
! Mean, SD, K (number of SD's) are given and returns BND(low,mid,up)

! Input and output variables
integer,intent(in):: length
real*8,intent(in) :: array(length)
real*8,intent(in)::ks
real*8,intent(out):: bnd(2)

real*8::avg,std
real*8::low,up,mid

if (ks.le.0.0) stop'Wrong Ks value'


call find_mean(array,length,avg)
call find_std(array,length,std)

!!$mid=avg
!!$low = mid - dble(ks)*std
!!$up  = mid + dble(ks)*std

bnd(1)=avg
bnd(2)=dble(ks)*std

!!$bnd(3)=up

return
end subroutine make_bound
!!$
!!$subroutine make_bound1(avg,std,ks,val)
!!$implicit none
!!$
!!$! Makes Error bar like upper and lower bound when
!!$! Mean, SD, K (number of SD's) are given and returns VAL(low,mid,up)
!!$
!!$real*8,intent(out):: val(3)
!!$real*8,intent(in)::avg,std
!!$integer,intent(in)::ks
!!$
!!$real*8::low,up,mid
!!$
!!$if (ks.eq.0.0) stop'Wrong Ks value'
!!$
!!$mid=avg
!!$low = mid - dble(ks)*std
!!$up  = mid + dble(ks)*std
!!$
!!$val(1)=low
!!$val(2)=mid
!!$val(3)=up
!!$
!!$return
!!$end subroutine make_bound1

subroutine  matrix_process(nruns)
  use dimKrig,only:rmsemat,loopcounter,outfile
  implicit none

  real*8::vec(nruns)
  real*8::nrows
  real*8::bnd(2)
  integer::i,j,k,nruns

  nrows=loopcounter
  bnd(:)=0.0


  outfile(1:4)='AVKR'
  open(23,file='norm/'//outfile,form='formatted',status='unknown')
  write(23,'(3a)') 'NPOINTS   ','   MeanRMSE   ','   BOUND'
  
  do i=1,nrows !nrows
     vec=rmsemat(1:nruns,i,2)
     call make_bound(nruns,vec,1.0,bnd)
     write(23,'(i8,2e15.8)')int(rmsemat(1,i,1)),bnd(1),bnd(2)
  end do
  close(23)

  return
end subroutine matrix_process

!!$program rmsebound
!!$implicit none
!!$
!!$integer::nruns
!!$real*8,allocatable(:,:):: rmsemat
!!$
!!$
!!$
!!$
!!$! Set the number of runs
!!$nruns=10
!!$! 
!!$
!!$allocate(rmsemat(nruns))
!!$
!!$end program rmsebound
!!$
!!$
!!$subroutine processmatrix(nrows,nruns,nfunc,matin,matout)
!!$implicit none
!!$
!!$integer::nruns
!!$integer::nrows
!!$integer::nfunc
!!$
!!$real*8::matin(nrows,nruns,nfunc)
!!$real*8::matout(nrows,nfunc,3)
!!$
!!$!1=low bound
!!$!2=mean
!!$!3=upper bound
!!$
!!$integer::i,j,k
!!$
!!$! Main program
!!$
!!$! Set some data
!!$
!!$
!!$nrows=5
!!$nruns=3
!!$nfunc=1
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$return
!!$end subroutine processmatrix
