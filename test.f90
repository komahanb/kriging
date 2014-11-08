program main

  implicit none

!  include 'mpif.h'
  integer :: ierr,i,j,imax,nstattmp,ndimtmp,NCP,lenc,ndimin,NMCin,fctindxin,Casemode
 integer :: nDIMint
  parameter (nDIMint=2)

  double precision :: diffconv,xavgin(ndimint),xstdin(ndimint),fmeanout,fvarout,fmeanprimeout(ndimint),fvarprimeout(ndimint)
  character*60 :: filename
  character*2 :: dimnumber,fctnumber
  character*3 :: lsnumber
  double precision :: freal,df(20),d2f(20,20),v(20),Javg,Jvar
  character*2 :: nptstoaddpercycnum

  integer:: nsamples,counter,numberpointstudy
  integer,parameter::timing=1 !other number if timing results not needed
  double precision :: Initialmach, Finalmach, InitialAOA,FinalAOA

  integer ::evlfnc,fuct
  common/evf/ evlfnc

  call MPI_START

ndimin=2
xavgin(:)=0.0

xstdin(:)=0.15d0
NMCin=100000
fctindxin=0

call Krigingestimate(ndimin,xavgin,ndimint,xstdin,fmeanout,fvarout,fmeanprimeout,fvarprimeout,NMCin,fctindxin)

print *, fmeanout,fvarout
call stop_all
end program main
