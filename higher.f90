        subroutine Post_Higher(ifac)
        use dimKrig
        implicit none
        include 'mpif.h'
        integer, intent(in) :: ifac
        integer :: i,j,k,istat,ierr
        integer :: idec,is,ie,imode,id
        integer :: ict,ictglb,idy,idyy,ide,idee
        integer, dimension(ndim) :: itp
        double precision :: yhat,RMSE,EI,maxerror,maxerrorglb
        double precision :: EImax,EIglb,yhatmin,yhatglb,fsum,fsumglb,f
        double precision, dimension(ndim) :: x,xe,xy,df,v,yhatprime
        double precision, dimension(ndim,ndim) :: d2f
        double precision, allocatable, dimension(:) :: ANV
        double precision, allocatable, dimension(:) :: ANL

        if(Cmode(:12).eq.'Post_HigherD')then
           imode = 1
        else
           imode = 0
!          if(id_proc.eq.0)then
!            allocate(ANV(ifac**ndim),stat=istat)
!            if(istat.ne.0)stop'allocation ANV'
!            ANV = 0.d0
!          end if
           allocate(ANL(ifac**ndim),stat=istat)
           if(istat.ne.0)stop'allocation ANL'
           ANL = 0.d0
        end if

        if(Cmode(:16).eq.'Post_ANOVA_Local')then
           if(nOPT.ne.1)stop'only for Mono-Objective'
           if(IDopt(1).le.0.or.IDopt(1).gt.nsample)stop'no optimal'
        end if

        do 100 k=1,1!nfunc-1

           idec = dble(ifac**ndim)/dble(num_proc)
           is   = idec*id_proc + 1
           ie   = idec*(id_proc+1)
           if(id_proc.eq.num_proc-1)ie = ifac**ndim
           write(*,'(11x,a,2i10,a,i3)') '>> [',is,ie,' ] for Processor',id_proc
           ict     = 0
           fsum    = 0.d0
           EImax   = 0.d0
           yhatmin = 1.d10
           maxerror= 0.d0

           do 200 i=is,ie

              call make_cartesian_higherD(i,ndim,ifac,itp,x)

              if(Cmode(:16).eq.'Post_ANOVA_Local')then
                 do 310 j=1,ndim
                    x(j) = sample(IDopt(1),j) + 0.01d0*(x(j)-0.5d0)
310              continue
              end if
!!$              x(1)=1.0d0
!!$              x(2)=1.0d0
!!$              x(3)=1.0d0
              call meta_call(k,imode,x,yhat,yhatprime,RMSE,EI)
              if(yhat.le.yhatmin)then
                 yhatmin = yhat
                 xy(:)   = x(:)
              end if
              if(Cmode(:10).eq.'Post_ANOVA')then
                 ANL(i) = yhat
              else if(Cmode(:12).eq.'Post_HigherD')then
                 if (fct.lt.20) call evalfunc(x,ndim,fct,0,0,f,df,d2f,v)

!!                 print *,'x:',X,'Exact:',f,'PC:',yhat
                 fsum = fsum + (yhat-f)**2
                 if (abs(f-yhat).gt.maxerror) maxerror=abs(f-yhat)
                 if(rmse.lt.0.d0)ict = ict + 1
                 if(EI.ge.EImax)then
                    EImax   = EI
                    xe(:)   = x(:)
                 end if
              end if
200        continue ! main loop for Cartesian (i)
! MPI
           if(Cmode(:10).eq.'Post_ANOVA')then
!            call MPI_GATHER(ANL(is),ie-is+1,MPI_DOUBLE_PRECISION,   &
!                            ANV( 1),ie-is+1,MPI_DOUBLE_PRECISION,0, &
!                            MPI_COMM_WORLD,ierr)
              do id=0,num_proc-1
                 is   = idec*id + 1
                 ie   = idec*(id+1)
                 if(id.eq.num_proc-1)ie = ifac**ndim
                 call MPI_BCAST(ANL(is),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
              end do
           end if
           call MPI_ALLREDUCE(ict,ictglb,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
           call MPI_ALLREDUCE(fsum,fsumglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
           fsumglb = dsqrt(fsumglb/dble(ifac**ndim)) 
           call MPI_ALLREDUCE(maxerror,maxerrorglb,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
           call MPI_ALLREDUCE(yhatmin,yhatglb,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
           call MPI_ALLREDUCE(EImax,EIglb,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
           idy = -1
           ide = -1
           if(yhatmin.eq.yhatglb) idy = id_proc
           if(  EImax.eq.  EIglb) ide = id_proc
           call MPI_ALLREDUCE(idy,idyy,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
           call MPI_ALLREDUCE(ide,idee,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
           if(idyy.lt.0.or.idyy.ge.num_proc) stop'idyy'
           if(idee.lt.0.or.idee.ge.num_proc) stop'idee'
           call MPI_BCAST(xy(1),ndim,MPI_DOUBLE_PRECISION,idyy,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(xe(1),ndim,MPI_DOUBLE_PRECISION,idee,MPI_COMM_WORLD,ierr)

           if(id_proc.eq.0)then
              write(*,'(1x,a,i3)')'>> Post Process for Function -',k
              write(*,'(6x,a,e15.5,a,99f8.3)') '>> Yhat Min = ',yhatglb,' @',(xy(i),i=1,ndim)
              if(Cmode(:12).eq.'Post_HigherD')then
                 write(*,'(6x,a,e15.5,a,99f8.3)') '>> EI   Max = ',  EIglb,' @',(xe(i),i=1,ndim)
                 write(*,'(6x,a,i8,a)')'>> Negative RMSE on ',ict,' Pts'
                 write(*,*)
                 write(*,'(6x,a,e20.10)') '>> RMSE compared to Analytical function = ',fsumglb
                 write(*,*)
                 write(93,'(i8,4e15.8)') nhs,fsumglb,maxerrorglb,diffloc,distloc
                 rmsemat(runnum,loopcounter,1)=nhs
                 rmsemat(runnum,loopcounter,2)=fsumglb
              end if


              if(Cmode(:10).eq.'Post_ANOVA')then
                 call ANOVA(ifac,ndim,ANL)
              end if
           end if
100     continue ! function loop (k)

        if(Cmode(:10).eq.'Post_ANOVA' .and. id_proc.eq.0) deallocate(ANL)
 
        end subroutine Post_Higher
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine make_cartesian_higherD(inp,ndim,ifac,itp,x)
        implicit none
        integer, intent(in) :: inp,ndim,ifac
        integer, dimension(ndim), intent(out) :: itp
        double precision, dimension(ndim), intent(out) :: x
        integer :: j,l,indx

        j  = inp

        do 100 l=1,ndim

           indx = int(j/ifac**(ndim-l))
           if(mod(j,ifac**(ndim-l)).eq.0) indx = indx - 1
           if(indx.eq.0.or.indx.eq.ifac-1)then
           else if(indx.ge.1.and.indx.le.ifac-2)then
           else
              stop'invalid indx'
           end if
           itp(l) = indx
           x(l)   = dble(indx)/dble(ifac-1)
           j      = j - indx*(ifac**(ndim-l))
           if (x(l).lt.0.d0.or.x(l).gt.1.d0) stop'x(l)>1 or x(l)<0'

100     continue

        end subroutine make_cartesian_higherD
