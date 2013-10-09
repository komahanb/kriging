        subroutine eva_bfgs(mode,n,xin,f,g)
        use dimKrig
        implicit none
        include 'mpif.h'
        integer, intent(in)                         :: mode,n
        double precision, dimension(n), intent(in)  :: xin
        double precision, intent(out)               :: f
        double precision, dimension(n), intent(out) :: g
        integer, dimension(0:num_proc-1,1000) :: ispl
        integer, dimension(0:num_proc-1)      :: ispm
        double precision, dimension(n)        :: xtmp
        double precision, dimension(ndim,tdim):: theta,power
        integer :: ierr,id,iall,i,k,l,imode
        double precision :: dx,ftmp,llfd,verr
        double precision :: pena,yhat,RMSE,EI,yhatprime(ndim)
        double precision, dimension(nfunc) :: ycon,econ

        if(id_proc.eq.0) &
        call spliting(n+1,num_proc,0,ispl,ispm)
        call MPI_BCAST(ispl(0,1),num_proc*1000,MPI_INTEGER,0, &
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ispm(0),num_proc,MPI_INTEGER,0, &
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(xin(1),n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        dx = 1.d-4
        do i=1,ispm(id_proc)
           id = ispl(id_proc,i)
           if(id.lt.0.or.id.gt.n)stop'id in eva_bfgs'
           xtmp(:) = xin(:)
           if(id.ne.0)xtmp(id) = xtmp(id) + dx

           if(mode.ge.1.and.mode.le.3)then
             if(mode.eq.1)then
               do l=1,tdim
                do k=1,ndim
                  theta(k,l)   = xtmp(          (l-1)*ndim+k)
                  if(iscf.eq.1)then
                    power(k,l) = xtmp(ndim*tdim+(l-1)*ndim+k)
                  else
                    power(k,l) = 2.d0
                  end if
                end do
               end do
             else if(mode.eq.2)then
               do l=1,tdim
                do k=1,ndim
                  theta(k,l)   = xtmp(     k)
                  if(iscf.eq.1)then
                    power(k,l) = xtmp(ndim+k)
                  else
                    power(k,l) = 2.d0
                  end if
                end do
               end do
             else if(mode.eq.3)then
               do l=1,tdim
                do k=1,ndim
                  theta(k,l)   = xtmp(1)
                  if(iscf.eq.1)then
                    power(k,l) = xtmp(2)
                  else
                    power(k,l) = 2.d0
                  end if
                end do
               end do
             end if
             call likelihood(theta,power,llfd,verr)
             if(verr.lt.1.d-3)then
               ftmp = -1.d0*llfd
             else
               ftmp = -1.d0*llfd + verr*cpena
             end if
           else if(mode.ge.10.and.mode.le.30)then
             if(n.ne.ndim)stop'eva_bfgs, ndim'
             if(mode.le.20)then
                imode = 0
             else
                imode = 1
             end if
             do k=1,nfunc-1
                call meta_call(k,imode,xtmp,yhat,yhatprime,RMSE,EI)
                ycon(k) = yhat
                econ(k) = EI
             end do
             call constraint(ycon,pena)
             if(mode.le.20)then
               if(pena.eq.0.d0)then
                 ftmp = ycon( nfOPT(mode-10) )
               else
                 ftmp = ycon( nfOPT(mode-10) ) + pena*10.d0
               end if
             else
               if(pena.eq.0.d0)then
                 ftmp = -1.d0*econ( nfOPT(mode-20) )
               else
                 ftmp = -1.d0*econ( nfOPT(mode-20) ) + pena*10.d0
               end if
             end if
           else
             stop'unknown mode in eva_bfgs'
           end if

           if(id.eq.0)then
             if(id_proc.ne.0)stop'why original is not node(0)?'
             f = ftmp
           else
             g(id) = ftmp
           end if
        end do

        call MPI_BCAST(f,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        do i=0,num_proc-1
           id   = ispl(i,1)
           iall = ispm(i)
           if(id.eq.0)then
             if(i.ne.0)stop'@ eva_bfgs1'
             id   = ispl(i,2)
             if(abs(id).ne.1)stop'@ eva_bfgs2'
             iall = iall -1
           end if
           if(iall.ne.0) &
           call MPI_BCAST(g(id),iall,MPI_DOUBLE_PRECISION,i, &
                          MPI_COMM_WORLD,ierr)
        end do

        do i=1,n
           g(i) = (g(i)-f)/dx
        end do

        end subroutine eva_bfgs
