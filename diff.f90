        subroutine diff_estimation(xin,dy,d2y,                          &
                                   id_proc,ndim,kreg_orig,              &
                                   iscf,mode_dck,                       &
                                   ccrf,                                &
                                   rsample,nsize,kreg,lmax,tdim,        &
                                   train,mean,devi,deviratio,           &
                                   theta,power,                         &
                                   RYFB,dxadd,hvec,                     &
                                   ict_sample,nparent,                  &
                                   inf)
        implicit none
        integer, intent(in) :: id_proc
        double precision, dimension(ndim),      intent(in)  :: xin
        double precision, dimension(ndim),      intent(out) :: dy
        double precision, dimension(ndim,ndim), intent(out) :: d2y

        integer,          intent(in) :: ndim,kreg,kreg_orig
        integer,          intent(in) :: rsample,nsize,iscf,mode_dck
        integer,          intent(in) :: lmax,tdim
        double precision, intent(in) :: devi,ccrf
        double precision, dimension(rsample,ndim),   intent(in) :: train
        double precision, dimension(kreg,1),         intent(in) :: mean
        double precision, dimension(lmax),           intent(in) :: deviratio
        double precision, dimension(ndim,tdim),      intent(in) :: theta,power
        double precision, dimension(nsize,1),        intent(in) :: RYFB
        double precision, dimension(rsample),        intent(in) :: dxadd
        integer,          dimension(0:10,0:10),      intent(in) :: ict_sample
        integer,          dimension(rsample),        intent(in) :: nparent
        character(len=6), dimension(rsample),        intent(in) :: inf
        double precision, dimension(rsample,ndim,2), intent(in) :: hvec

        integer :: i,k,l
        double precision, dimension(1,1)       :: y1,y2
        double precision, dimension(1,kreg)    :: freg
        double precision, dimension(nsize,1)   :: r
        double precision, dimension(kreg_orig) :: yout

        dy  = 0.d0
        d2y = 0.d0

! for 1st-order gradient
        do 100 k=1,ndim
          call makeMatrixPetitR(1,k,0, &
               ndim,nsize,rsample,iscf,mode_dck, &
               lmax,tdim,ict_sample, &
               train,xin,nparent,dxadd,inf,hvec, &
               theta,power,ccrf,r)
          call get_regression(1,k,0,xin,yout)
          freg = 0.d0
          do i=1,kreg_orig
            freg(1,i) = yout(i)
          end do
          y1    = matmul(freg,mean)               ! (1*k)x(k*1)
          y2    = matmul(transpose(r),RYFB)       ! (1*ns)x(ns*1)
          dy(k) = y1(1,1) + y2(1,1)
100     continue

! for 2nd-order gradient
        do 200 k=1,ndim
        do 200 l=k,ndim
          call makeMatrixPetitR(2,k,l, &
               ndim,nsize,rsample,iscf,mode_dck, &
               lmax,tdim,ict_sample, &
               train,xin,nparent,dxadd,inf,hvec, &
               theta,power,ccrf,r)
          call get_regression(2,k,l,xin,yout)
          freg = 0.d0
          do i=1,kreg_orig
            freg(1,i) = yout(i)
          end do
          y1 = matmul(freg,mean)               ! (1*k)x(k*1)
          y2 = matmul(transpose(r),RYFB)       ! (1*ns)x(ns*1)
          d2y(k,l) = y1(1,1) + y2(1,1)
          d2y(l,k) = d2y(k,l)
200     continue

        end subroutine diff_estimation
