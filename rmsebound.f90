
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

subroutine make_bound(length,array,ks,bnd)
implicit none

! Makes Error bar like upper and lower bound when
! Mean, SD, K (number of SD's) are given and returns BND(low,mid,up)

! Input and output variables
integer,intent(in):: length
real*8,intent(in) :: array(length)
real*8,intent(in)::ks
real*8,intent(out):: bnd(7)

real*8::avg,std
real*8::low,up,mid

if (ks.le.0.0) stop'Wrong Ks value'


call find_mean(array,length,avg)
call find_std(array,length,std)

bnd(1)=avg

bnd(2)=minval(array) !minimum of all observed errors
bnd(3)=maxval(array) !maximum of all observed errors

bnd(4)=abs(bnd(1)-bnd(2)) 
bnd(5)=abs(bnd(3)-bnd(1))

bnd(6)=abs(bnd(2)-bnd(1))/dlog(10.0)
bnd(7)=abs(bnd(3)-bnd(1))/dlog(10.0)

return
end subroutine make_bound

subroutine  matrix_process(nruns)
  use dimKrig,only:rmsemat,loopcounter,outfile
  implicit none

  real*8::vec(nruns)
  real*8::nrows
  real*8::bnd(7)
  integer::i,j,k,nruns

  nrows=loopcounter
  bnd(:)=0.0

  outfile(1:4)='AVKR'
  open(23,file='norm/'//outfile,form='formatted',status='unknown')
  write(23,'(3a)') 'NPOINTS   ','   MeanRMSE   ','   BOUND'
  
  do i=1,nrows !nrows
     vec=rmsemat(1:nruns,i,2)
     call make_bound(nruns,vec,1.0,bnd)
     write(23,'(i8,7e15.8)')int(rmsemat(1,i,1)),bnd(1),bnd(2),bnd(3),bnd(4),bnd(5),bnd(6),bnd(7)
  end do
  close(23)

  return
end subroutine matrix_process

