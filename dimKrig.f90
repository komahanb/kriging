      Module dimKrig
      implicit none

      integer::loopcounter,runnum
! mpi parallel
      integer :: id_proc, num_proc,OUUflag

      integer::probtype(20) ! 1=fixed interval, 2= varialbe interval

! problem set
      character(len=20) :: Cmode
      integer :: nOPT,nCON,fct,hstat,lstat,nstyle,iterDEL,randomini,randomtestl,dynamics,ifid,reusesamples,maxsample,maxsamplewant
      integer :: nKRI,nDCK,nICK,nDMF,nCOK,nVFM,nptstoaddpercyc,Dutchorderg,readMCsamples,selectedevaluation,evlfnc
      integer :: nMXS,nMXG,nMXH,nMXV,ndiffloc,NMCS,fctindx,filenum,ctr
      integer, dimension(10)            :: nfOPT,nfCON
      integer, dimension(10)            :: nfKRI,nfDCK,nfICK,nfDMF,nfCOK,nfVFM
      double precision fmean,fvar
      double precision, dimension(10)   :: vCON
      double precision, dimension(2,20)   :: DS
      double precision, dimension(20)   :: xavg,xstd,fmeanprime,fvarprime,xavgt,xstdt
      double precision, dimension(1000)   :: difflocar
      character(len=1), dimension(10)   :: cCON1
      character(len=7), dimension(10)   :: cCON2

! direct
      logical:: mainprog
      logical:: lhsdyn
      integer :: randomflag
      integer:: mode_dck
      
! indirect(trust)
      integer :: itrust,ngput,nhput,it4opt
      character(len=1) :: Cnug
      double precision :: Vnug
      double precision, dimension(2) :: tdxinit
      double precision, dimension(4) :: tdd_thres
      double precision :: tdxmin,tdxmax,dpar,pkl,pku,vgfac,mgfac,bdfac

! Kriging
      integer :: ndim,ndimt,nfunc,nhes
      integer :: zsample,nsample,rsample ! L+H, H, H-D12
      integer :: nhs,nls
      integer :: lmax,tdim
      integer :: ndebug,iran,newflag,mEI,mLC
      integer :: nsize,kreg,kord,kreg_orig,iscf
      double precision :: devi,devmin,devmax,diffloc,distloc
      double precision :: llfd_best
      integer, allocatable, dimension(:,:) :: mreg,iRij
      double precision, allocatable, dimension(:) &
              :: yy,deviratio,deviratio_best
      double precision, allocatable, dimension(:,:) &
              :: Rij,Rinv,mean,r,  FB,RYFB,FRFinv,FR, &
                 d_theta,d_power
      integer, dimension(0:10,0:10) :: ict_sample
      integer, dimension(0:10     ) :: ict_dummy
      double precision :: ccrf,Reps
      double precision,dimension(20)::dat
! All Krig
      integer :: i_dim,d_dim,j_dim,c_dim
      integer, allocatable, dimension(:,:) &
              :: ipv,iak,ibl,icm
      double precision, allocatable, dimension(:,:) &
              :: dak
      integer, allocatable, dimension(:,:) &
              :: jbl
      character(len=6), allocatable, dimension(:,:) &
              :: ccm
      integer :: ncr1,ncr2
      double precision, allocatable, dimension(:,:) &
              :: dadd,dmul

! opt
      double precision :: cpena
      integer :: lopt 
      integer, dimension(10) :: IDopt
      double precision, dimension(10) :: OFopt

! basic
      character(len=6), allocatable, dimension(:) &
              :: info, inf
      integer, allocatable, dimension(:) &
              :: icone, &
                 nparent
      double precision, allocatable, dimension(:) &
              :: tdxx, &
                 tdx,dxadd,fun
      double precision, allocatable, dimension(:,:) &
              :: sample,func, &
                 sampl ,gfun
      double precision, allocatable, dimension(:,:,:) &
              :: gfunc, &
                 hfun,hvec
      double precision, allocatable, dimension(:,:,:,:) &
              :: hfunc,hvect

      ! Ensemble
      double precision,allocatable,dimension(:,:,:)::rmsemat
      character*60 :: outfile

      End Module dimKrig
