subroutine Make_Sample
  use dimKrig
  implicit none
  include 'mpif.h'
  common/global/counter
  integer :: idx,is,ie,idec,id,ierr
  integer :: nseed,mode,counter,nhstmp
  integer :: i,j,k
  double precision :: f,ran,prd,factor,x(ndim)
  double precision, allocatable, dimension(:)   :: df,hv,v,ftmp,fdata
  double precision, allocatable, dimension(:,:)   :: xtmp,dftmp,hvtmp,vtmp,dfdata
  double precision, allocatable, dimension(:,:) :: d2f
  double precision, allocatable, dimension(:,:,:) :: d2ftmp,d2fdata
  double precision, allocatable, dimension(:,:) :: bound
  character(len=17) :: Cfile
  character(len=6)  :: Cstat,Cstatl
  character(len=6), allocatable,dimension(:)  :: Cstattmp
  character*2 :: fctindxnumber
  character*60 :: filename
  logical :: pointinbox

  allocate(sample(ndim,nhs))
  allocate(sampl(ndim,nls))
  allocate(df(ndim))
  allocate(d2f(ndim,ndim))    
  allocate(hv(ndim))
  allocate(v(ndim))
  allocate(bound(2,ndim))

  allocate(xtmp(ndim,300))
  allocate(ftmp(300))
  allocate(dftmp(ndim,300))
  allocate(d2ftmp(ndim,ndim,300))    
  allocate(hvtmp(ndim,300))
  allocate(vtmp(ndim,300))
  allocate(Cstattmp(300))

  ! for parallel function calls
  allocate(fdata(nhs))
  allocate(dfdata(nhs,ndim))
  allocate(d2fdata(nhs,ndim,ndim))

  nhstmp=0

  ! Cstatus
  if(     hstat.eq.0)then
     Cstat = '1_F   '
  else if(hstat.eq.1)then
     Cstat = '1_FG  '
  else if(hstat.eq.3)then
     Cstat = '1_FGH '
  else if(hstat.eq.5)then
     Cstat = '1_FGHv'
  else
     stop'unknown hstat'
  end if

  if (nls.gt.0) then

     ! Cstatusl
     if(     lstat.eq.0)then
        Cstatl = '2_F   '
     else if(lstat.eq.1)then
        Cstatl = '2_FG  '
     else if(lstat.eq.3)then
        Cstatl = '2_FGH '
     else if(lstat.eq.5)then
        Cstatl = '2_FGHv'
     else
        stop'unknown lstat'
     end if

  end if


  if (reusesamples.eq.1 .and. randomini.eq.1) then

     filename='Function'
     call i_to_s(fctindx,fctindxnumber)
     filename(9:10)=fctindxnumber
     filename(11:14)='.dat'

     if(id_proc.eq.0)  open(97,file=filename,form='formatted',status='unknown')

     do while (nhstmp.lt.2*maxsamplewant)

        if(id_proc.eq.0)  read(97,100,end=19) Cstat,(x(j),j=1,ndim),f,(df(j),j=1,ndim),&
             ((d2f(k,j),j=1,ndim),k=1,ndim)

        call MPI_BCAST(Cstat,6,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(x,ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(f,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(df,ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(d2f,ndim*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        print*,cstat,x,f,df,d2f,id_proc

        pointinbox=.true.
        do j=1,ndim
           if (x(j).gt.DS(2,j) .or. x(j).lt.DS(1,j)) pointinbox=.false.
        end do

        if (pointinbox) then
           nhstmp=nhstmp+1
           Cstattmp(nhstmp)=Cstat
           ftmp(nhstmp)=f
           do j=1,ndim
              xtmp(j,nhstmp)=(x(j)-DS(1,j))/(DS(2,j)-DS(1,j))
              dftmp(j,nhstmp)=df(j)*(DS(2,j)-DS(1,j))
           end do
           if (hstat.le.3) then
              do j=1,ndim
                 do k=1,ndim
                    d2ftmp(k,j,nhstmp)=d2f(k,j)*(DS(2,j)-DS(1,j))*(DS(2,k)-DS(1,k))
                 end do
              end do
           else if (hstat.gt.3) then
              do j=1,ndim
                 vtmp(j,nhstmp)=v(j)
                 d2ftmp(j,1,nhstmp)=d2f(j,1)
              end do
           end if
        end if

     end do

19   if(id_proc.eq.0) close(97)

     ! Make sure to do at least one cycle of adaption
     maxsample=max(maxsamplewant,nhstmp+nptstoaddpercyc)

     if(id_proc.eq.0) write (filenum,'(a,i3,a,i3,a,i3)') 'Could reuse',nhstmp,&
          ' samples and still need',maxsample-nhstmp,&
          ' samples for function',fctindx

  end if ! reuse samples






  if (nhstmp.lt.nhs) then

     if (id_proc.eq.0) then

        if (randomini.eq.0) then

           if(id_proc.eq.0) write(filenum,*)'>> Initial Sample Points in corners'

           if (randomtestl.eq.1) then
              bound(1,:)=0.0d0!0.05
              bound(2,:)=1.0d0!0.95
           else
              bound(1,:)=0.0d0
              bound(2,:)=1.0d0
           end if

           sample(:,1)=0.5d0
           counter=2
           call recurcorn(1,ndim,sample,nhs,bound)

        else if (randomini.eq.1) then



           if (randomflag.eq.1) then

              write(filenum,*)'>> Initial Sample Points by Latin Hypercube'
              sample(:,1)=0.5
              call get_seed(nseed)
              call latin_random(ndim,nhs-1-nhstmp,nseed,sample(:,2:nhs-nhstmp))       

           else if (randomflag.eq.2) then

              write(filenum,*)'>> Initial Sample Points by NIEDER Sequence'
              call get_seed(nseed)
              call nieder(nseed,ndim,nhs-nhstmp,sample(:,1:nhs-nhstmp))
              !                 print*,sample(:,1:nhs-nhstmp)
              !                 stop    
           else if (randomflag.eq.3) then

              write(filenum,*)'>> Initial Sample Points by Halton Sequence'
              call halton_real(ndim,nhs-nhstmp,sample(:,1:nhs-nhstmp))
              !                print*,sample(:,1:nhs-nhstmp)
              !                stop
           else if (randomflag.eq.4) then

              if(id_proc.eq.0)        write(filenum,*)'>> Initial Sample Points by Hammersley Sequence'
              call hammersley_real(ndim,nhs-nhstmp,sample(:,1:nhs-nhstmp))
              !                 print*,sample(:,1:nhs-nhstmp)
              !                stop

           else if (randomflag.eq.5) then

              write(filenum,*)'>> Initial Sample Points by Sobol Sequence'
              call get_seed(nseed)
              call sobol_real(nseed,ndim,nhs-nhstmp,sample(:,1:nhs-nhstmp))
              !                 print*,sample(:,1:nhs-nhstmp)
              !                stop

           else if (randomflag.eq.6) then

              write(filenum,*)'>> Initial Sample Points by Faure Sequence'
              call get_seed(nseed)
              call faure_real(nseed,ndim,nhs-nhstmp,sample(:,1:nhs-nhstmp))
              !                print*,sample(:,1:nhs-nhstmp)
              !                stop

           else !if (randomflag.eq.6) then

              if(id_proc.eq.0)     write(filenum,*)'>> Initial Sample Points by Which Sequence'
              print*,"Cool.. Please go ahead and implement"
              stop

           end if !randomflag

        end if !random ini

     end if !MASTER

     call MPI_BCAST(sample(1:ndim,1:nhs-nhstmp),ndim*(nhs-nhstmp),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     if (nls.gt.0) then

        if(id_proc.eq.0) then
           write(filenum,*)'>> Initial Low-fidelity Sample Points by Latin Hypercube'
           call get_seed(nseed)
           call latin_random(ndim,nls,nseed,sampl)       
        end if

        call MPI_BCAST(sampl(1:ndim,1:nls),ndim*nls,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     end if

     if(nstyle.eq.0)then

        if (id_proc.eq.0) then
           open(10,file='sample.dat',form='formatted',status='unknown')
           !$$ write(10,'(3i8)') ndim,nhs+nls,3
           write(10,'(3i8)') ndim,nhs+nls,2
        end if !master
     end if

     do i=1,nhstmp ! write the available points

        if (hstat.le.3.and.id_proc.eq.0) write(10,100) Cstattmp(i),(xtmp(j,i),j=1,ndim),&
             ftmp(i),0.d0,(dftmp(j,i),j=1,ndim),&
             ((d2ftmp(k,j,i),j=1,ndim),k=1,ndim)  
        if (hstat.gt.3.and.id_proc.eq.0) write(10,100) Cstattmp(i),(xtmp(j,i),j=1,ndim),&
             ftmp(i),0.d0,(dftmp(j,i),j=1,ndim),&
             ((d2ftmp(k,j,i),j=1,ndim),k=1,ndim),(vtmp(j,i),j=1,ndim),&
             (d2ftmp(k,1,i),k=1,ndim)  
     end do

     !====================================================!
     !================ WORK SHARE ========================!
     !====================================================!

     idec = dble(nhs-nhstmp+nls)/dble(num_proc)
     is   = idec*id_proc + 1
     ie   = idec*(id_proc+1)
     if(id_proc.eq.num_proc-1)ie = nhs-nhstmp+nls

     do i=is,ie

        if(i.le.nhs-nhstmp)then
           x(:) = sample(:,i)
           call evalfunc(x,ndim,fct,0,hstat,fdata(i),dfdata(i,:),d2fdata(i,:,:),v)
           write(filenum,*) "Training data #",i,"by proc",id_proc,"f: ",fdata(i)
        else
           x(:) = sampl(:,i-nhs-nhstmp)   
           !ifid=2
           call evalfunc(x,ndim,fct,ifid,lstat,fdata(i),dfdata(i,:),d2fdata(i,:,:),v)
           write(filenum,*) "Training data #",i,"by proc",id_proc,"f: ",fdata(i)
        end if

     end do !is,ie

     ! Information Sharing
     do id=0,num_proc-1
        is   = idec*id + 1
        ie   = idec*(id+1)
        if(id.eq.num_proc-1)ie = nhs-nhstmp+nls

        call MPI_BCAST(fdata(is),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)

        do i=1,ndim
           call MPI_BCAST(dfdata(is,i),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
        end do

        do i=1,ndim
           do j =1,ndim
              call MPI_BCAST(d2fdata(is,i,j),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
           end do
        end do

     end do ! all the procs

     !====================================================!
     !===========  END WORK SHARE ========================!
     !====================================================!

     f=0.0d0
     df=0.0d0
     d2f=0.0d0

        do i=1,nhs-nhstmp+nls


           f=fdata(i)
           df(:)=dfdata(i,:)
           d2f(:,:)=d2fdata(1,:,:)

           if(nstyle.eq.0)then ! with func


              if (i.le.nhs) then ! high-fid

                 x(:)=sample(i,:)

                 !$$ if (hstat.le.3) write(10,100) Cstat,(x(j),j=1,ndim),f,f,0.d0,(df(j),j=1,ndim),((d2f(k,j),j=1,ndim),k=1,ndim)   !0.d0: everything is fine     other: point corrupted
                 !$$ if (hstat.gt.3) write(10,100) Cstat,(x(j),j=1,ndim),f,f,0.d0,(df(j),j=1,ndim),(v(j),j=1,ndim),(d2f(k,1),k=1,ndim)   !0.d0: everything is fine     other: point corrupted
                 if (hstat.le.3.and.id_proc.eq.0) write(10,100) Cstat,(x(j),j=1,ndim),&
                      f,0.d0,(df(j),j=1,ndim),((d2f(k,j),j=1,ndim),k=1,ndim)   !0.d0: everything is fine     other: point corrupted
                 if (hstat.gt.3.and.id_proc.eq.0) write(10,100) Cstat,(x(j),j=1,ndim),&
                      f,0.d0,(df(j),j=1,ndim),(v(j),j=1,ndim),(d2f(k,1),k=1,ndim)   !0.d0: everything is fine     other: point corrupted
              else ! low-fid
                 
                 x(:) = sampl(:,i-nhs-nhstmp)   

                 !$$ if (lstat.le.3) write(10,100) Cstatl,(x(j),j=1,ndim),f,f,0.d0,(df(j),j=1,ndim),((d2f(k,j),j=1,ndim),k=1,ndim)   !0.d0: everything is fine     other: point corrupted
                 !$$ if (lstat.gt.3) write(10,100) Cstatl,(x(j),j=1,ndim),f,f,0.d0,(df(j),j=1,ndim),(v(j),j=1,ndim),(d2f(k,1),k=1,ndim)   !0.d0: everything is fine     other: point corrupted
                 if (lstat.le.3.and.id_proc.eq.0) write(10,100) Cstatl,(x(j),j=1,ndim),&
                      f,0.d0,(df(j),j=1,ndim),((d2f(k,j),j=1,ndim),k=1,ndim)   !0.d0: everything is fine     other: point corrupted
                 if (lstat.gt.3.and.id_proc.eq.0) write(10,100) Cstatl,(x(j),j=1,ndim),&
                      f,0.d0,(df(j),j=1,ndim),(v(j),j=1,ndim),(d2f(k,1),k=1,ndim)   !0.d0: everything is fine     other: point corrupted
              end if

           else

              if(i.ge.1.and.i.le.9)then
                 write(Cfile,111)i
              else if(i.ge.10.and.i.le.99)then
                 write(Cfile,112)i
              else
                 write(Cfile,113)i
              end if
              open(10,file=Cfile,form='formatted',status='unknown')
              write(10,101)(x(j),j=1,ndim),dble(hstat)
              close(10)

           end if !nstyle.eq.0

        end do ! loop over points

     else

        nhs=nhstmp

        if(nstyle.eq.0.and.id_proc.eq.0)then
           open(10,file='sample.dat',form='formatted',status='unknown')
           !$$ write(10,'(3i8)') ndim,nhs+nls,3
           write(10,'(3i8)') ndim,nhs+nls,2
        end if

        do i=1,nhs
           if (hstat.le.3.and.id_proc.eq.0) write(10,100) Cstattmp(i),(xtmp(j,i),j=1,ndim),&
                ftmp(i),0.d0,(dftmp(j,i),j=1,ndim),((d2ftmp(k,j,i),j=1,ndim),k=1,ndim)  
           if (hstat.gt.3.and.id_proc.eq.0) write(10,100) Cstattmp(i),(xtmp(j,i),j=1,ndim),&
                ftmp(i),0.d0,(dftmp(j,i),j=1,ndim),((d2ftmp(k,j,i),j=1,ndim),k=1,ndim)&
                ,(vtmp(j,i),j=1,ndim),(d2ftmp(k,1,i),k=1,ndim)  
        end do

  end if  ! nhstmp.lt.nhs ?

  if (id_proc.eq.0) then

     if(nstyle.eq.0)then   
        close(10)
        write(filenum,'(a,i5)')'      >> Output to sample.dat',nhs
     end if

  end if

  deallocate(bound)
  deallocate(sample)
  deallocate(sampl)
!  deallocate(x)
  deallocate(df,d2f,hv,v)
  deallocate(xtmp,ftmp,dftmp,d2ftmp,hvtmp,vtmp,Cstattmp)
  deallocate(fdata,dfdata,d2fdata)

100 format(a,1x,10000e20.10)
101 format(10000e20.10)
111 format('newsample00',i1,'.dat')
112 format('newsample0',i2,'.dat')
113 format('newsample',i3,'.dat')


end subroutine Make_Sample



  recursive subroutine recurcorn(dimact,ndim,sample,nhs,bounds)

    implicit none

    common/global/counter

    integer :: dimact,ndim,nhs,counter,j,counterin
    real*8 :: sample(ndim,nhs),bounds(2,ndim)

    counterin=counter
    do j=2,1,-1
       sample(dimact,counter)=bounds(j,dimact)
       if (dimact.gt.1) then
          sample(1:dimact-1,counter)=sample(1:dimact-1,counterin)
       end if
       if (dimact.lt.ndim) call recurcorn(dimact+1,ndim,sample,nhs,bounds)
       if (j.ne.1) counter=counter+1
    end do


  end subroutine recurcorn







