subroutine MonteCarlo
  use dimKrig
  implicit none
  include 'mpif.h'
  integer :: i,j,k,ii,kk
  integer :: idec,is,ie,id
  integer :: ierr,istat,seed
  integer :: ntgt,ftgt,npdf
  double precision :: st,pst,dx,p1,p2,xmintmp,xmaxtmp,dinvnorm
  double precision :: yhat,RMSE,EI,ran
  double precision, dimension(ndim)      :: x,df,x0,Ddiff,sigma,xbar,v,yhatprime
  double precision, dimension(ndim,ndim) :: d2f
  integer :: ict,ictglb,ihst
  double precision :: ymin,ymax,yminglb,ymaxglb,yhmin,yhmax
  double precision :: MCm,MCmglb,MCd,MCdglb,hglb,width,pdf
  double precision, dimension(ndim)      :: xmin,xmax,xminglb,xmaxglb,MCmprime,MCmprimeglb,MCdprime,MCdprimeglb
  double precision, allocatable, dimension(:)   :: MNCf
  double precision, allocatable, dimension(:,:) :: MNCx
  double precision, allocatable, dimension(:,:) :: HST,HSTglb
  character(len=6) :: cdum
  real*8::average(ndim)
  double precision :: f,muy1,muy2,sigmay1,sigmay2,fobjlin,fobjquad,Javg(3),Jvar(3),freal
  real*8::fvtemp
  character*60 :: histname
  character*3  :: fctindxnumber
  integer :: expensive
  logical writeMCSamples

  writeMCSamples=.false.

  open(10,file='MC.inp',form='formatted',status='unknown')
  read(10,*)! (xavg(i),i=1,ndim)
  read(10,*)! (xstd(i),i=1,ndim)     
  read(10,*)
  read(10,*)
  read(10,*) NMCS!,ndimtmp
  read(10,*) npdf
  read(10,*) readmcsamples
  read(10,*) evlfnc
  read(10,*) expensive
  close(10)

  if (readmcsamples.eq.1) writeMCSamples=.true.

  call find_Optimal
  if(id_proc.eq.0) write(filenum,'(1x,a)')'>> MonteCarlo Simulation on Kriging Model'

  ihst = 0

  if (readMcsamples.eq.0 .and. NMCS.eq.0) then
     open(10,file='MC.inp',form='formatted',status='unknown')
     read(10,*)(xavg(i),i=1,ndim)
     read(10,*,err=11) pst,cdum,dx
     call find_st(pst,dx,st)
     call MPI_BCAST(st,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     go to 12
11   continue
     !       write(*,*)pst,cdum,dx
     backspace 10
     read(10,*) (xstd(i),i=1,ndim)
     pst = 0.9    ! dummy
     dx  = 0.03d0 ! dummy
12   continue
     read(10,*) ntgt
     read(10,*) ftgt
     read(10,*) NMCS
     read(10,*,err=13) npdf, yhmin, yhmax
     go to 14
13   continue
     backspace 10
     read(10,*) npdf
     ihst = 1
14   continue
     close(10)
     !       write(*,*)st,ntgt,ftgt,NMCS
     if(pst.le.0.d0.or.pst.ge.1.d0)stop'pst in MC.inp'
     if(dx.le.0.d0.or.dx.ge.1.d0)stop'dx in MC.inp'
     if(ntgt.lt.0.or.ntgt.gt.ndim)stop'ntgt in MC.inp'
     if(ftgt.le.0.or.ftgt.gt.nfunc-1)stop'ftgt in MC.inp'

  else

     ntgt=0
     ftgt=1
!     npdf=30

  end if

  allocate( MNCf(NMCS), stat=istat)
  if(istat.ne.0)stop'allocation MNCf in MonteCarlo'
  if(ntgt.eq.0 .or. readMcsamples.eq.1)then
     allocate( MNCx(ndim,NMCS), stat=istat)
  else
     allocate( MNCx(1,NMCS), stat=istat)
  end if
  if(istat.ne.0)stop'allocation MNCx in MonteCarlo'
  allocate( HST(npdf,3), stat=istat)
  if(istat.ne.0)stop'allocation HST in MonteCarlo'
  allocate( HSTglb(npdf,3), stat=istat)
  if(istat.ne.0)stop'allocation HSTglb in MonteCarlo'
  MNCf   = 0.d0
  MNCx   = 0.d0
  HST    = 0.d0
  HSTglb = 0.d0

  if(id_proc.eq.0)then
!!$           write(*,'(6x,a,e15.5)')'>> Standard Deviation = ',st
!!$           call CDF(          -1.d0*dx,0.d0,st,p1)
!!$           call CDF_Numerical(-1.d0*dx,0.d0,st,p2)
!!$           write(*,'(6x,a,2e15.5,a,f12.7)') &
!!$           '>> CDF = ',p1,p2,' @ dx =  ',-1.d0*dx
     write(filenum,'(6x,a,i3)')'>> Target Function = ',ftgt
     if(ntgt.ne.0)then
        write(filenum,'(6x,a,i3)')'>> Target Design Variable = ',ntgt
     else
        write(filenum,'(6x,a)')   '>> All Design Variables Fluctuated'
     end if
     write(filenum,'(6x,a,i8)')'>> # of Monte-Carlo Samples = ',NMCS

     if (readMcsamples.eq.1) then

        if (fct.lt.20) then
           open(10,file='MCsamp.dat',form='formatted',status='unknown')
           read(10,'(2i8)') NMCS,ndim
           read(10,*) (xavg(i),i=1,ndim)
           read(10,*) (xstd(i),i=1,ndim)
           do j=1,NMCS
              read(10,*) (MNCx(k,j),k=1,ndim)
              !                    print *, j, MNCX(:,j)
           end do
           close(10)

        else ! fct =10
!!$
!!$           open(10,file='MC.inp',form='formatted',status='unknown')
!!$           read(10,*)! (xavg(i),i=1,ndim)
!!$           read(10,*)! (xstd(i),i=1,ndim)     
!!$           read(10,*)
!!$           read(10,*)
!!$           read(10,*) !NMCS!,ndimtmp
!!$           close(10)

           open(17,file='MCvalues.dat',form='formatted',status='unknown')
           do j=1,nmcs
              read(17,*) (MNCx(kk,j),kk=1,ndim),fvtemp  !Just read the cost function
           end do
           close(17)

        end if

     else

        call get_seed(seed)
        call latin_random(ndim,NMCS,seed,MNCx)
        
        if (writeMCSamples) then

           open(10,file='MCsamp.dat',form='formatted',status='unknown')
           write(10,'(2i8)') NMCS,ndim
           write(10,*) (xavg(i),i=1,ndim)
           write(10,*) (xstd(i),i=1,ndim)

        end if

        do j = 1, NMCS
           do k=1,ndim 
!!$              
!!$              if (fct.eq.24) then !wing optimization
!!$
!!$                 if (k.eq.3) then
!!$
!!$                    MNCx(k,j)=xavg(k)+dinvnorm(MNCx(k,j))*xstd(k)
!!$
!!$                 else
!!$
!!$                    MNCx(k,j)=xavg(k)+MNCx(k,j)*xstd(k)
!!$
!!$                 end if
!!$
!!$              else
!!$
              MNCx(k,j)=xavg(k)+dinvnorm(MNCx(k,j))*xstd(k)
!!$
!!$              end if

           end do
           if (writeMCSamples) write(10,*) (MNCx(kk,j),kk=1,ndim)
        end do
        if (writeMCSamples) close(10)


     end if
     
     do j = 1, NMCS
        do k=1,ndim
           MNCx(k,j) = (MNCx(k,j)-DS(1,k))/(DS(2,k)-DS(1,k))
        end do
     end do

  end if !id_proc.eq.0

  if(ntgt.eq.0)then
     call MPI_BCAST(MNCx(1,1),NMCS*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     stop
     call MPI_BCAST(MNCx(1,1),NMCS*1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if

  do k=ftgt,ftgt

     idec = dble(NMCS)/dble(num_proc)
     is   = idec*id_proc + 1
     ie   = idec*(id_proc+1)
     if(id_proc.eq.num_proc-1)ie = NMCS
     write(filenum,'(11x,a,2i10,a,i3)') '>> [',is,ie,' ] for Processor',id_proc
     ymin = 1.d10
     ymax =-1.d10
     xmin = 1.d10
     xmax =-1.d10
     MCm  = 0.d0
     MCd  = 0.d0
     MCmprime(:)=0.d0
     MCdprime(:)=0.d0
     ict  = 0

     do i=is,ie ! Main Loop for MC

        x(:)    = MNCx(:,i)                           
        !     print *,i,x(:)
        do j=1,ndim
           xmin(j) = min(xmin(j),x(j))
           xmax(j) = max(xmax(j),x(j))
        end do

!!$              xmintmp=0.0
!!$              xmaxtmp=1.0
!!$              do j=1,ndim
!!$                 xmintmp=min(xmintmp,x(j))
!!$                 xmaxtmp=max(xmaxtmp,x(j))
!!$              end do
!!$
!!$              if ( xmaxtmp.gt.1.0 .or. xmintmp.lt.0.0  ) then
!!$                 call meta_call(2,0,x,yhat,yhatprime,RMSE,EI) 
!!$              else
!!$                 call meta_call(1,0,x,yhat,yhatprime,RMSE,EI)
!!$              end if

        !call meta_call(1,0,x,yhat,yhatprime,RMSE,EI)
        call meta_call(1,2,x,yhat,yhatprime,RMSE,EI)

        !              call evalfunc(x,ndim,fct,0,1,freal,df,d2f,v)
        !!              print *, 'x:',x!,'Kr:',yhat,'Ex:',freal ,id_proc 

!!$              print *,'Derivative at',x
!!$              print *,yhatprime
!!$              print *,df
!!$              print *
!!$              ict    = ict + 1
!!$              MCm    = MCm + freal    ! for mean
!!$              MCd    = MCd + freal**2 ! for variance              
!!$              do j=1,ndim
!!$                 MCmprime(j) = MCmprime(j) + df(j)  !for derivative of mean
!!$                 MCdprime(j) = MCdprime(j) + freal*df(j)  !for derivative of variance
!!$              end do

        ict    = ict + 1
        MCm    = MCm + yhat    ! for mean
        MCd    = MCd + yhat**2 ! for variance              
        do j=1,ndim
           MCmprime(j) = MCmprime(j) + yhatprime(j)  !for derivative of mean
           MCdprime(j) = MCdprime(j) + yhat*yhatprime(j)  !for derivative of variance
        end do

        MNCf(i) = yhat
        if(yhat.ge.ymax)then
           ymax = yhat
        end if
        if(yhat.le.ymin)then
           ymin = yhat
        end if

     end do! main loop for MonteCarlo (i)     

     ! Information Sharing
     do id=0,num_proc-1
        is   = idec*id + 1
        ie   = idec*(id+1)
        if(id.eq.num_proc-1)ie = NMCS
        call MPI_BCAST(MNCf(is),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
     end do
     call MPI_ALLREDUCE(ict,ictglb,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MCm,MCmglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MCd,MCdglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(ymin,yminglb,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(ymax,ymaxglb,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
     do i=1,ndim
        call MPI_ALLREDUCE(xmin(i),xminglb(i),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(xmax(i),xmaxglb(i),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MCmprime(i),MCmprimeglb(i),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MCdprime(i),MCdprimeglb(i),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     end do

     if(ictglb.ne.NMCS)stop'ictglb.ne.NMCS'

     fmean = MCmglb / dble(ictglb)
     fvar = MCdglb / dble(ictglb) - fmean**2  
     do i=1,ndim
        fmeanprime(i) = MCmprimeglb(i) / dble(ictglb)
        fvarprime(i) = 2.0 * MCdprimeglb(i) / dble(ictglb) - 2.0 * fmean * fmeanprime(i)
     end do

     ! Histogram

     yhmin = yminglb
     yhmax = ymaxglb
     width = (yhmax-yhmin)/dble(npdf)
     do i=1,npdf
        HST(i,1) = yhmin + dble(i-1)*width
        HST(i,2) = yhmin + dble(i  )*width
     end do
     HST(npdf,2) = yhmax
     HST(:,3)    = 0.d0
     MCd = 0.d0
     is   = idec*id_proc + 1
     ie   = idec*(id_proc+1)
     if(id_proc.eq.num_proc-1)ie = NMCS
     do 250 i=is,ie
        if(MNCf(i).le.yhmin)then
           HST(1,   3) = HST(1,   3) + 1.d0
        else if(MNCf(i).ge.yhmax)then
           HST(npdf,3) = HST(npdf,3) + 1.d0
        else
           do 260 j=1,npdf
              if(MNCf(i).ge.HST(j,1).and.MNCf(i).lt.HST(j,2))then
                 HST(j,3) = HST(j,3) + 1.d0
                 go to 250
              end if
260           continue
              write(*,'(2i8,3e15.5)')id_proc,i,yhmin,MNCf(i),yhmax
              stop'without the border of histogram?'
           end if

250        continue

           do i=1,npdf
              call MPI_ALLREDUCE(HST(i,3),hglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
              HSTglb(i,1) = HST(i,1)
              HSTglb(i,2) = HST(i,2)
              HSTglb(i,3) = hglb
           end do

           ! output
           if(id_proc.eq.0)then
              write(filenum,*)
              write(filenum,'(6x,a,6e15.5)')'>> xavg ',xavg(1:ndim)
              write(filenum,'(6x,a,6e15.5)')'>> xstd ',xstd(1:ndim)   
              write(filenum,'(6x,a,i4)')'>> Fctindx ',fctindx
              write(filenum,*)

              write(filenum,'(6x,a,2e15.5)')'>> Mean and Variance = ',fmean,fvar
              write(filenum,'(6x,a,20e15.5)')'>> Meanprime = ',fmeanprime(1:ndim)/(DS(2,1:ndim)-DS(1,1:ndim))
              write(filenum,'(6x,a,20e15.5)')'>> Varprime = ',fvarprime(1:ndim)/(DS(2,1:ndim)-DS(1,1:ndim))
              write(filenum,'(6x,a,99f8.3)')'>> Range of X(min) = ',(xminglb(i),i=1,ndim)
              write(filenum,'(6x,a,99f8.3)')'>> Range of X(max) = ',(xmaxglb(i),i=1,ndim)
              write(filenum,'(6x,a,2e15.5)')'>> Ymax/min = ',ymaxglb,yminglb
              write(filenum,'(6x,a,2(e15.5,a))')'>> Output Histogram  [',yhmin,' :',yhmax,' ]'

              histname(1:8)='HISTGidx'
              call i_to_s(fctindx,fctindxnumber)
              histname(9:11)=fctindxnumber
              histname(12:15)='.dat'

              pdf = 0.d0

              open(10,file=histname,form='formatted',status='unknown')
              do i=1,npdf
                 pdf = pdf + HSTglb(i,3)/dble(NMCS)
                 write(10,'(4e15.5)')(HSTglb(i,1)+HSTglb(i,2))*0.5d0,HSTglb(i,3),HSTglb(i,3)/dble(NMCS),pdf
              end do
              close(10)
           end if

        end do! dummy loop for function index (k)


        ! Plotting hack

!!$        if (ndim.eq.2) then
!!$
!!$           if (id_proc.eq.0) then
!!$              open(10,file='MCout.dat',form='formatted',status='unknown')
!!$              do j=1,NMCS
!!$                 if (fct.lt.20) call evalfunc(MNCx(:,j),ndim,fct,0,0,freal,df,d2f,v)
!!$                 write(10,*) (MNCx(k,j)*(DS(2,k)-DS(1,k))+DS(1,k),k=1,ndim),freal
!!$              end do
!!$              close(10)              
!!$           end if
!!$ 
!!$           Cmode='Post_2D'
!!$           call Post_1or2D(101)
!!$           Cmode='Post_MonteCarlo'
!!$           
!!$        end if

        ! Moment methods, extrapolation and real function evaluation comparisons if desired

        if (id_proc.eq.0) then

           if (expensive.ne.1) then ! if not expensive then execute the rest

           write(filenum,*)
           write(filenum,*)

           do k=1,ndim
              sigma(k) = xstd(k)/(DS(2,k)-DS(1,k))
              xbar(k) = (xavg(k)-DS(1,k))/(DS(2,k)-DS(1,k))
           end do

           if (hstat.gt.0) then

              f=func(1,1)
              df(:)=gfunc(1,nfCOK(1),:)
              d2f(:,:)=hfunc(1,nfCOK(1),:,:)

              !MM1
              muy1=f
              sigmay1=0.0
              do k=1,ndim
                 sigmay1=sigmay1+(sigma(k)*df(k))**2
              end do

              if (hstat.eq.3) then

                 ! MM2
                 muy2=muy1
                 do k=1,ndim
                    muy2=muy2+0.5*sigma(k)**2*d2f(k,k)
                 end do

                 sigmay2=sigmay1
                 do j=1,ndim
                    do k=1,ndim
                       sigmay2=sigmay2+0.5*(sigma(j)*sigma(k)*d2f(j,k))**2
                    end do
                 end do

              end if

           end if



           Javg(:)=0.0
           Jvar(:)=0.0

           if (fct.lt.20.or.fct.eq.22) then

              open(10,file='MCvalues.dat',form='formatted',status='unknown')


           else if (fct.eq.20.or.fct.eq.22)then

              if (fctindx.eq.0) then

                 open(10,file='MCCFDvalues00.dat',form='formatted',status='unknown')

              else if (fctindx.eq.4) then

                 open(10,file='MCCFDvalues04.dat',form='formatted',status='unknown')

              end if
           end if


          
           do i=1,NMCS

              ! Real function evaluation

              if (evlfnc.eq.1) then

                 call evalfunc(MNCx(:,i),ndim,fct,0,0,freal,df,d2f,v) 
!                 print*,MNCx(:,i)*(DS(2,1)-DS(1,1))+DS(1,1),freal

!                call evalfunc(MNCx(:,i),ndim,6,0,0,freal,df,d2f,v) 
!               print*,MNCx(:,i),freal
!              stop
                 write(10,*) (MNCx(kk,i),kk=1,ndim),freal 

              else

                 read(10,*) (MNCx(kk,i),kk=1,ndim),freal               

              end if

              Javg(1)=Javg(1)+freal
              Jvar(1)=Jvar(1)+freal**2

              if (hstat.gt.0) then

                 Ddiff(:)=MNCx(:,i)-xbar(:)

                 ! Linear extrapolation        
                 fobjlin=f
                 do k=1,ndim
                    fobjlin=fobjlin+df(k)*Ddiff(k)
                 end do

                 Javg(2)=Javg(2)+fobjlin
                 Jvar(2)=Jvar(2)+fobjlin**2

                 if (hstat.eq.3) then

                    ! Quadratic extrapolation
                    fobjquad=fobjlin
                    do j=1,ndim
                       do k=1,ndim
                          fobjquad=fobjquad+0.5*d2f(j,k)*Ddiff(j)*Ddiff(k)
                       end do
                    end do

                    Javg(3)=Javg(3)+fobjquad
                    Jvar(3)=Jvar(3)+fobjquad**2

                 end if

              end if

           end do
           close(10)


           Javg(:)=Javg(:)/real(NMCS)      
           Jvar(:)=Jvar(:)/real(NMCS)-Javg(:)**2
!!$           if (fct.eq.20) then
!!$              if (xstd(1).eq.0.01) then
!!$                 Javg(1)= 0.25906E+00    
!!$                 Jvar(1)= 0.31653E-01
!!$              else if (xstd(1).eq.0.0075) then
!!$                 Javg(1)= 0.26361E+00    
!!$                 Jvar(1)= 0.18123E-01
!!$              else if (xstd(1).eq.0.005) then
!!$                 Javg(1)= 0.26671E+00    
!!$                 Jvar(1)= 0.86688E-02
!!$              else if (xstd(1).eq.0.0025) then
!!$                 Javg(1)= 0.26819E+00    
!!$                 Jvar(1)= 0.20390E-02
!!$              end if
!!$           end if

           write(filenum,'(6x,a,2e15.5)')'Real: Mean and Variance',Javg(1),Jvar(1)
           write(filenum,*)

           write(filenum,'(6x,a,2e15.5)')'Krig: Mean and Variance',fmean,fvar
           write(filenum,*)

           write(filenum,'(6x,a,2e15.5)')'Error: Mean and Variance',abs(fmean-Javg(1)),abs(fvar-Jvar(1))
           write(filenum,*)

           !           if (hstat.eq.0)
           write(94,'(i4,6e16.8)') nhs,Javg(1),Jvar(1),fmean,fvar,abs(fmean-Javg(1)),abs(fvar-Jvar(1))

!!$           if (hstat.gt.0) then
!!$           
!!$              write(filenum,'(6x,a,2e15.5)')'MM1: Mean and Variance',muy1,sigmay1
!!$              write(filenum,*)
!!$           
!!$              write(filenum,'(6x,a,2e15.5)')'Lin: Mean and Variance',Javg(2),Jvar(2)
!!$              write(filenum,*)
!!$              
!!$              if (hstat.ne.3) write(94,'(i4,8e16.8)') nhs,Javg(1),Jvar(1),fmean,fvar,abs(fmean-Javg(1)),abs(fvar-Jvar(1)),muy1,sigmay1,Javg(2),Jvar(2)
!!$           
!!$              if (hstat.eq.3) then
!!$                 write(filenum,'(6x,a,2e15.5)')'MM2: Mean and Variance',muy2,sigmay2
!!$                 write(filenum,*)
!!$                 write(filenum,'(6x,a,2e15.5)')'Quad: Mean and Variance',Javg(3),Jvar(3)
!!$                 write(filenum,*)
!!$                 write(94,'(i4,12e16.8)') nhs,Javg(1),Jvar(1),fmean,fvar,abs(fmean-Javg(1)),abs(fvar-Jvar(1)),muy1,sigmay1,Javg(2),Jvar(2),muy2,sigmay2,Javg(3),Jvar(3)
!!$              end if
!!$
        end if
     end if
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     !!     end if  !readMCsamples.ne.0 ?
     
     deallocate(MNCf,MNCx) 
     deallocate(HST,HSTglb) 
     
   end subroutine MonteCarlo
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine find_st(pst,dx,st)
        implicit none
! find st which ensures the probability within [-dx:dx] is pst
        double precision, intent(in)  :: pst,dx
        double precision, intent(out) :: st
        integer :: iter
        double precision :: pout,s,ds,vv,vp,dv,sini

        if(pst.le.0.d0.or.pst.ge.1.d0)stop'pst in find_st'
        pout = (1.d0-pst)/2.d0
        sini = 1.d0
190     continue
        s    = sini
        ds   = sini*1.d-3
        iter = 0
200     continue
        iter = iter + 1
        call CDF(-1.d0*dx,0.d0,s,   vv)
        if(dabs(vv-pout).le.1.d-10)go to 210
        call CDF(-1.d0*dx,0.d0,s+ds,vp)
        dv = (vp-vv)/ds
!       write(*,'(5e15.5)')s,vv,pout,vv-pout,dv
        if(dv.eq.0.d0)stop'dv = 0.d0 in find_st'
        if(iter.ge.100)stop'iter>100 in find_st'
        s = s - (vv-pout)/dv
        if(s.le.0.d0)then
          sini = sini * 0.1d0
          go to 190
        end if
        go to 200
210     continue
!       write(*,'(4e15.5)')s,vv,pout,vv-pout
        st = s

        end subroutine find_st
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine find_x(mode,ran,xc,st,xout)
        implicit none
! find xout which ensures the CDF at xout is ran
! mode=0 : analytical CDF (fast but less robust?)
! mode=1 : numerical  CDF (time comsuming, but robust)
        integer, intent(in) :: mode
        double precision, intent(in)  :: ran,xc,st
        double precision, intent(out) :: xout
        integer :: iter 
        double precision :: x,vv,dv,vp,dx

        x    = xc
        dx   = 1.d-4
        iter = 0
200     continue
        iter = iter + 1
        if(mode.eq.0)then
          call CDF(x,xc,st,vv)
        else
          call CDF_Numerical(x,xc,st,vv)
        end if
        if(dabs(vv-ran).le.1.d-10)go to 210
        if(mode.eq.0)then
          call CDF(x+dx,xc,st,vp)
        else
          call CDF_Numerical(x+dx,xc,st,vp)
        end if
        dv = (vp-vv)/dx
!       write(*,'(4e15.5)')x,vv,vv-ran,dv
        if(dv.eq.0.d0)stop'dv=0 in find_x'
        if(iter.ge.100)stop'iter>100 in find_x'
        x = x - (vv-ran)/dv
        go to 200
210     continue
!       write(*,'(3e15.5)')x,vv,vv-ran
        xout = x

        end subroutine find_x
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine CDF(xin,xc,st,vout)
        implicit none
        double precision, intent(in)  :: xin,xc,st
        double precision, intent(out) :: vout
        double precision :: vtmp
!       vout = 0.5d0 * (1.d0 + erf( (xin-xc)/(st*sqrt(2.d0)) ))
        call ERF_MINE1( (xin-xc)/(st*sqrt(2.d0)), vtmp )
        vout = 0.5d0 * (1.d0 + vtmp)
        end subroutine CDF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine DCDF(xin,xc,st,dvout)
        implicit none
        double precision, intent(in)  :: xin,xc,st
        double precision, intent(out) :: dvout
        double precision :: dvtmp
        call DERF_MINE( (xin-xc)/(st*sqrt(2.d0)), dvtmp )
        dvout = 0.5d0*dvtmp/(st*sqrt(2.d0))
        end subroutine DCDF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine ERF_MINE1(xin,yout)
        implicit none
        double precision, intent(in)  :: xin
        double precision, intent(out) :: yout
        integer :: i,k,n
        double precision :: vsum,kai
! n is the order of Taylor
! Maybe accurate within the range of [-4:4] with n=100
        n = 100
        vsum = 0.d0
        do i=0,n
          kai = 1.d0
          do k=1,i
            kai = kai * (-1.d0) * xin**2 / dble(k)
          end do
          vsum = vsum + kai*xin/dble(2*i+1)
        end do
        yout = vsum*2.d0/(sqrt(3.141592653589793238d0))

        if(yout.gt.1.d0)write(*,'(a,2e15.5)')'*ERF>1 ',xin,yout-1.d0
        end subroutine ERF_MINE1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine DERF_MINE(xin,dyout)
        implicit none
        double precision, intent(in)  :: xin
        double precision, intent(out) :: dyout
        double precision :: vsum

        vsum  = exp(-1.d0*xin**2)
        dyout = vsum*2.d0/(sqrt(3.141592653589793238d0))

        end subroutine DERF_MINE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine CDF_Numerical(xin,xc,st,cdf)
        implicit none
        double precision, intent(in)  :: xin,xc,st
        double precision, intent(out) :: cdf
        double precision :: vtmp
        integer :: i,num
        double precision :: xs,xe,dx,x1,x2,pdf1,pdf2
        if(xin.lt.xc)then
          cdf = 0.d0
          xs  = xin -2.d0
          xe  = xin
        else if(xin.ge.xc)then
          cdf = 0.5d0
          xs  = xc
          xe  = xin
        end if
        num = 1001
        dx  = (xe-xs)/dble(num-1)
        do i=1,num-1
           x1 = xs + dble(i-1)*dx
           x2 = xs + dble(i  )*dx
           call normal_dist(x1,xc,st,pdf1)
           call normal_dist(x2,xc,st,pdf2)
           cdf = cdf + (pdf1+pdf2)*dx*0.5d0
        end do
        end subroutine CDF_Numerical
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine normal_dist(xin,xc,st,y)
        implicit none
        double precision, intent(in)  :: xin,xc,st
        double precision, intent(out) :: y
        double precision :: pi

        pi = 4.d0*datan(1.d0)
        y = exp(-1.d0*(xin-xc)**2/2.d0/st/st)/sqrt(2.d0*pi*st**2)

        end subroutine normal_dist

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function dinvnorm(p)

    implicit none

      real*8 :: dinvnorm,p,p_low,p_high
      real*8 :: a1,a2,a3,a4,a5,a6
      real*8 :: b1,b2,b3,b4,b5
      real*8 :: c1,c2,c3,c4,c5,c6
      real*8 :: d1,d2,d3,d4
      real*8 :: z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=sqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=sqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   dinvnorm=z

      return

    end function dinvnorm
