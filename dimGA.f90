      Module dimGA
      implicit none

! mpi parallel
      integer :: id_proc, num_proc
! evaluation counting
      integer :: ict_eva(10),ictg_eva(10)

! input parameters
      integer :: npop, ngen, ndv, nobj, ndat, neva
      character (len= 3) :: Cini,Cobj(10)
      character (len=30) :: Cinifile,Ctype_sh,Cross,Cmut,Cprob,Cpareto
      integer :: ngen_sh,ngen_mat,ngen_stop,ndebug, &
                 ngen_out,ngen_fin,ngen_grad,ntype_sh,neva_max
      integer :: id_sus,if_crs,iaxs_nec
      double precision :: &
                 alpha_sh,beta_cr,r_mut,eta_mut,fac_prob, &
                 p_mut,p_muts,p_mute, &
                 p_drx,p_drxs,p_drxe, &
                 p_nec,p_necs,p_nece

! major dimensionals
      integer :: ipool
      integer, dimension(10,2) :: elite
      double precision, dimension(10) :: Fmin,Fmax
      double precision, allocatable, dimension(:,:,:) :: d,parent
      double precision, allocatable, dimension(:,:) :: rank,Frank,Rrank
      integer, allocatable, dimension(:,:) :: pool
! tmp
      double precision, allocatable, dimension(:)   :: bd_t
      double precision, allocatable, dimension(:,:) :: dv_t,dv_p
      double precision :: vrange
      integer :: ntmp,ttmp

      End Module dimGA
