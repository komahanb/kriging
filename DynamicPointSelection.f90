subroutine DynamicPointSelection
  use dimKrig
  implicit none
  include 'mpif.h'

  common/global/counter
  integer :: counter,i,ii,j,jj,jjj,k,kk,kp,l,NTOEX,NTOEXtmp,triangle_node(ndim+1,100000),triangle_num,triangle_coor_num,NCP,node,knnptr(20000),orderextmp(0:20000),Dutchorder(100000),nseed,hstatad(nptstoaddpercyc)
  integer :: mode
  double precision :: triangle_coor(ndim,100000),Dtoextmp(ndim),Dtoex(ndim,100000),dist(100000),minftoex(100000),maxftoex(100000),distmean,ftoextry(2,100000),derivdummy(ndim)
  double precision :: diff2,RMSE(100000),RMSEmean,EI,Ddibtmp(ndim,0:1000),Dgdibtmp(ndim,0:1000),fdibtmp(0:1000),gdibtmp(ndim,0:1000),hdibtmp(ndim,ndim,0:1000),diffloctmp,difflocmin,difflocavg, SIGMAmean,distcomp

  double precision, dimension(nptstoaddpercyc) :: f
  double precision, dimension(ndim,nptstoaddpercyc) :: df,Dad,v
  double precision, dimension(ndim,ndim,nptstoaddpercyc) :: d2f
  double precision, DIMENSION(200) :: SIGV
  double precision, DIMENSION(200) :: SIGG
  integer :: Taylororder, IERR, NCPG,idec,is,ie,id,point,kpc,nptstoaddpercyctmp,  nptstoaddpercycorig
  double precision :: BETA, GAMM, SIGMA(100000),distmeanloc,RMSEmeanloc,minftoexloc(10000),maxftoexloc(10000),phi,invphi
  character*60::export


  integer,parameter::makesamples=1 			!1=Make random samples, 0=read from Kriging Samples

  integer::npass

  call find_Optimal
  call read_all_krig

  if (lhsdyn) then

     if (id_proc.eq.0) then

        write(filenum,*) '>> [Picking points dynamically for LHS]'        
        call get_seed(nseed)
        !        call nieder(nseed,ndim,nptstoaddpercyc,Dtoex)
        call latin_random(ndim,nptstoaddpercyc,nseed,Dtoex) 
        !   end if

        !     call MPI_Barrier(MPI_COMM_WORLD,ierr)
        !    call MPI_BCAST(Dtoex(:,:),100000,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

        do ii=1,nptstoaddpercyc
           hstatad(ii)=hstat
           if(nstyle.eq.0) then
              Dad(:,ii)=Dtoex(:,ii)
              call evalfunc(Dtoex(:,ii),ndim,fct,0,hstatad(ii),f(ii),df(:,ii),d2f(:,:,ii),v(:,ii))
           end if
           !           print*, DTOEX(:,ii),f(ii)
        end do !! ii loop
        ! Add successful test candidates to sample points 
        if(nstyle.eq.0)then
           open(10,file='sample.dat',form='formatted',status='unknown')
           !$$ write(10,'(3i8)')ndim,nhs+nls+nptstoaddpercyc,3
           write(10,'(3i8)')ndim,nhs+nls+nptstoaddpercyc,2
           do i=1,nhs+nls
              if (info(i)(3:6).ne.'FGHv') then
                 write(10,101) info(i),(sample(i,j),j=1,ndim),(func(i,j),j=1,nfunc),((gfunc(i,nfCOK(j),k),k=1,ndim),j=1,nCOK),(((hfunc(i,nfCOK(j),k,l),l=1,ndim),k=1,ndim),j=1,nCOK)
              else
                 write(10,101) info(i),(sample(i,j),j=1,ndim),(func(i,j),j=1,nfunc),((gfunc(i,nfCOK(j),k),k=1,ndim),j=1,nCOK),((hvect(i,nfCOK(j),k,1),k=1,ndim),j=1,nCOK),((hvect(i,nfCOK(j),k,2),k=1,ndim),j=1,nCOK)
              end if
           end do
           do ii=1,nptstoaddpercyc
              if (hstatad(ii).le.3) then
                 !$$ write(10,102) hstatad(ii),(Dad(j,ii),j=1,ndim),f(ii),f(ii),0.d0,(df(j,ii),j=1,ndim),((d2f(k,j,ii),j=1,ndim),k=1,ndim)
                 write(10,102) hstatad(ii),(Dad(j,ii),j=1,ndim),f(ii),0.d0,(df(j,ii),j=1,ndim),((d2f(k,j,ii),j=1,ndim),k=1,ndim)
              else
                 !$$ write(10,102) hstatad(ii),(Dad(j,ii),j=1,ndim),f(ii),f(ii),0.d0,(df(j,ii),j=1,ndim),(v(j,ii),j=1,ndim),(d2f(k,1,ii),k=1,ndim)
                 write(10,102) hstatad(ii),(Dad(j,ii),j=1,ndim),f(ii),0.d0,(df(j,ii),j=1,ndim),(v(j,ii),j=1,ndim),(d2f(k,1,ii),k=1,ndim)
              end if
           end do
           close(10)

        end if !nstyle

     end if

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call deallocate_all_krig

     return ! returns the routine
  end if

  phi=(1+sqrt(5.0d0))/2.0d0 !1.618
  invphi=1.0d0/phi

  if (randomtestl.eq.0) then
     if(id_proc.eq.0) then !master thread alone

        ! Make grid via Delaunay triangulation       
        triangle_coor_num=nhs
        do ii=1,nhs
           triangle_coor(1:ndim,ii)=sample(ii,1:ndim)
        end do


        open(unit=44,file='qdelaunayin')
        write(44,*) ndim
        write(44,*) triangle_coor_num
        do ii=1,triangle_coor_num
           write(44,*) triangle_coor(1:ndim,ii)
        end do
        close(44)

        call system('../qdelaunay QJ1e-3 i < qdelaunayin > qdelaunayout')

        open(unit=44,file='qdelaunayout')
        read(44,*) triangle_num
        do ii=1,triangle_num
           read(44,*) triangle_node(1:ndim+1,ii)
           triangle_node(1:ndim+1,ii)=triangle_node(1:ndim+1,ii)+1
        end do
        close(44)

        call fixcolindelaunay(triangle_num,triangle_coor_num,ndim,triangle_coor,triangle_node)

        if (ndim.eq.2) call triangulation_order3_plot('triangulation_plot.eps',triangle_coor_num,triangle_coor,triangle_num,triangle_node,2,2)

        ! Figure out locations of test candidates (midpoints of Delaunay sides and centres of Delaunay triangles),  Calculate local Dutch intrapolations and compare to Kriging values
        NCP=ndim+1
        NTOEX=0
        do ii=1,triangle_num

           NTOEXtmp=NTOEX+1

           !print *, 'Work on triangle',ii,'out of',triangle_num

           do kk=1,ndim+1

              kp=kk+1
              if (kp.gt.ndim+1) kp=1

              do jj=1,ndim             
                 Dtoextmp(jj)=(triangle_coor(jj,triangle_node(kk,ii))+triangle_coor(jj,triangle_node(kp,ii)))/2.0
              end do

              ! Don't consider doubles
              do jj=1,NTOEX
                 do jjj=1,ndim
                    if (Dtoextmp(jjj).ne.Dtoex(jjj,jj)) GOTO 111
                    if (jjj.eq.ndim) GOTO 112
                 end do
111              continue
              end do
              NTOEX=NTOEX+1
              Dtoex(1:ndim,NTOEX)=Dtoextmp(1:ndim)
112           continue

           end do

           NTOEX=NTOEX+1
           do jj=1,ndim
              Dtoex(jj,NTOEX)=0.0
              do kk=1,ndim+1
                 Dtoex(jj,NTOEX)=Dtoex(jj,NTOEX)+triangle_coor(jj,triangle_node(kk,ii))
              end do
              Dtoex(jj,NTOEX)=Dtoex(jj,NTOEX)/real(ndim+1)
           end do

           knnptr(1:ndim+1)=triangle_node(1:ndim+1,ii)

           do j=0,NCP-1  
              node=knnptr(j+1)
              Ddibtmp(1:ndim,j)=sample(node,1:ndim)
              fdibtmp(j)=func(node,1)
              orderextmp(j)=0
              if (info(node).eq.'FG ' .or. info(node).eq.'FGH ') then
                 gdibtmp(1:ndim,j)=gfunc(node,nfCOK(1),1:ndim)
                 orderextmp(j)=1
              end if
              if (info(node).eq.'FGH ' .or. info(node).eq.'FH ') then
                 hdibtmp(:,:,j)=hfunc(node,nfCOK(1),:,:)
                 orderextmp(j)=2
              end if
           end do

           Dutchorder(NTOEXtmp:NTOEX)=1

           call Dutch(Ddibtmp,fdibtmp,gdibtmp,hdibtmp,orderextmp,Dutchorder(NTOEXtmp:NTOEX),Dtoex(:,NTOEXtmp:NTOEX),ftoextry(1,NTOEXtmp:NTOEX),NCP,ndim,NTOEX-NTOEXtmp+1)

           do k=NTOEXtmp,NTOEX

              call meta_call(1,0,Dtoex(:,k),ftoextry(2,k),derivdummy,RMSE(k),EI)

              minftoex(k)=ftoextry(1,k)
              maxftoex(k)=ftoextry(1,k)
              do j=2,2
                 if (ftoextry(j,k).gt.maxftoex(k)) then
                    maxftoex(k)=ftoextry(j,k)
                 else if (ftoextry(j,k).lt.minftoex(k)) then
                    minftoex(k)=ftoextry(j,k)
                 end if
              end do

              !call evalfunc(Dtoex(:,k),ndim,fct,0,0,f(1),df(:,1),d2f(:,:,1),v(:,1))
              !print *, k,ftoextry(:,k),f(1)

           end do

        end do

     end if ! master thread only







     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     !PARALLEL REGION

     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




  else if (randomtestl.eq.1) then


     if(id_proc.eq.0)  write(filenum,*) '>> [Dutch Intrapolation is being used as local surrogate]'

     ! Figure out locations of test candidates (randomly). Calculate local Dutch intrapolations and compare to Kriging values

     call combination(ndim+Dutchorderg,ndim,NCP)
     !     NTOEX=(30-ndim)*NCP

     NTOEX=5000*NDIM

     if (id_proc.eq.0) then
        call get_seed(nseed)
        call latin_random(ndim,NTOEX,nseed,Dtoex) 
     end if

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_BCAST(Dtoex(:,:),100000,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

     !Information sharing by master with slaves
     idec = dble(NTOEX)/dble(num_proc)
     is   = idec*id_proc + 1
     ie   = idec*(id_proc+1)
     if(id_proc.eq.num_proc-1)ie =NTOEX  

     write(filenum,'(11x,a,2i10,a,i3)') '>> [',is,ie,' ] for Processor',id_proc           

     do k=is,ie !Main loop for test candidates

        ! if(id_proc.eq.2)   print *, Dtoex(:,k), id_proc

        ! Still need to make sure points are not collinear in higher dimensions!
        call knn(Dtoex(:,k),sample,knnptr,ndim,nhs,NCP)

        do j=0,NCP-1  
           node=knnptr(j+1)
           Ddibtmp(1:ndim,j)=sample(node,1:ndim)
           fdibtmp(j)=func(node,1)
           orderextmp(j)=0
           if (info(node).eq.'FG ' .or. info(node).eq.'FGH ') then
              gdibtmp(1:ndim,j)=gfunc(node,nfCOK(1),1:ndim)
              orderextmp(j)=1
           end if
           if (info(node).eq.'FGH ' .or. info(node).eq.'FH ') then
              hdibtmp(:,:,j)=hfunc(node,nfCOK(1),:,:)
              orderextmp(j)=2
           end if
        end do

        Dutchorder(k)=Dutchorderg

        !call DutchRBF(Ddibtmp,fdibtmp,gdibtmp,hdibtmp,orderextmp,Dutchorder(k),Dtoex(:,k),ftoextry(1,k),NCP,ndim,1)

        call Dutchgeninterp(Ddibtmp,fdibtmp,gdibtmp,hdibtmp,orderextmp,Dutchorder(k),Dtoex(:,k),ftoextry(1,k),NCP,ndim,1)

        !mode=0 ! return function value only
        mode=1 ! return function, RMSE, EI

        call meta_call(1,mode,Dtoex(:,k),ftoextry(2,k),derivdummy,RMSE(k),EI)
        !call evalfunc(Dtoex(:,k),ndim,fct,0,0,f(1),df(:,1),d2f(:,:,1),v(:,1))


     end do ! loop over NTOEX test candidates


  else if (randomtestl.eq.2) then     

     if(id_proc.eq.0)  write(filenum,*) '>> [MIR is being used as local surrogate]'

     call mirtunableparams(fct,ndim,nhs,ncp,taylororder)

     NTOEX=27000 !*NDIM

     !      NTOEX=int((1000*num_proc)/ndim)

     !     NTOEX=1000

     !     print*,taylororder

     if(id_proc.eq.0)  write(filenum,*) '     >> Number of test candidates',NTOEX
     if (id_proc.eq.0) then
        call get_seed(nseed)
        call latin_random(ndim,NTOEX,nseed,Dtoex) 
        !        call hammersley_real(ndim,NTOEX,DTOEX)
        !        write(export, '(a,i3.3,a)')'testcand.dat'
        !        open(10,file='./KrigSamples/'//export,form='formatted',status='unknown')
        !        read(10,*)(Dtoex(:,i),i=1,NToex)
        !        close(10)
     end if

     ! Information sharing by master with slaves        
     call MPI_Barrier(MPI_COMM_WORLD,ierr)           
     call MPI_BCAST(Dtoex(:,:),100000,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
     idec = dble(NTOEX)/dble(num_proc)
     is   = idec*id_proc + 1
     ie   = idec*(id_proc+1)
     if(id_proc.eq.num_proc-1)ie =NTOEX 
     write(filenum,'(11x,a,2i10,a,i3)') '>> [',is,ie,' ] for Processor',id_proc            

     do k=is,ie !Main loop for test candidates

        !if(id_proc.eq.2)   print *, Dtoex(:,k), id_proc
        ! Still need to make sure points are not collinear in higher dimensions!
        call knn(Dtoex(:,k),sample,knnptr,ndim,nhs,NCP)

        NCPG=0
        do j=0,NCP-1  
           node=knnptr(j+1)
           Ddibtmp(1:ndim,j)=sample(node,1:ndim)
           fdibtmp(j)=func(node,1)
           orderextmp(j)=0
           if (info(node).eq.'FG ' .or. info(node).eq.'FGH ') then 
              Dgdibtmp(1:ndim,NCPG)=sample(node,1:ndim) 
              gdibtmp(1:ndim,NCPG)=gfunc(node,nfCOK(1),1:ndim)  
              NCPG=NCPG+1
           end if
        end do

        ! Calculate the best parameters beta and gamma
        sigv=0.d0
        sigg=0.d0


        CALL MIR_BETA_GAMMA(nfunc-1, ndim, NCP, Ddibtmp(:,0:NCP-1), fdibtmp(0:NCP-1), SIGV, NCPG , Dgdibtmp(:,0:NCPG-1), gdibtmp(:,0:NCPG-1), SIGG, Taylororder, 1, dble(1.0), BETA, GAMM, IERR)
        if (ierr.ne.0) stop'MIR BETA gamma error'
        CALL MIR_EVALUATE(nfunc-1, ndim, 1, Dtoex(:,k), NCP, Ddibtmp(:,0:NCP-1), fdibtmp(0:NCP-1), SIGV, NCPG , Dgdibtmp(:,0:NCPG-1), gdibtmp(:,0:NCPG-1), SIGG, BETA, GAMM, Taylororder, 1, ftoextry(1,k), SIGMA(k), IERR)
        if (ierr.ne.0) stop'MIR evaluate error'

        !mode=0 ! return function value only
        mode=1 ! return function, RMSE, EI

        call meta_call(1,mode,Dtoex(:,k),ftoextry(2,k),derivdummy,RMSE(k),EI)

        !        print *, k,ftoextry(2,k),ftoextry(1,k)

     end do
  end if ! randomtestl


  ! Information Sharing --Exchange the function values
  do id=0,num_proc-1
     is   = idec*id + 1
     ie   = idec*(id+1)
     if(id.eq.num_proc-1)ie = NTOEX
     call MPI_BCAST(ftoextry(1,is:ie),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(ftoextry(2,is:ie),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(RMSE(is:ie),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
     if(randomtestl.eq.2)   call MPI_BCAST(SIGMA(is:ie),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
  end do
  !print *, ftoextry(1,1:NTOEX), ftoextry(2,1:NTOEX)
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !END PARALLEL REGION

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  if (id_proc.eq.0) then

     do k =1,NTOEX
        minftoex(k)=ftoextry(1,k)
        maxftoex(k)=ftoextry(1,k)
        if (ftoextry(2,k).gt.maxftoex(k)) then
           maxftoex(k)=ftoextry(2,k)
        else 
           minftoex(k)=ftoextry(2,k)
        end if
     end do

     ! Figure out distance to closest real sample point and the mean of all these distances
     distmean=0.0
     RMSEmean=0.0
     if(randomtestl.eq.2)  SIGMAmean=0.0
     do k=1,NTOEX
        dist(k)=1000000000000.0
        do kk=1,nhs
           diff2=0.0
           do jj=1,ndim
              diff2=diff2+(Dtoex(jj,k)-sample(kk,jj))**2
           end do
           if ( diff2.lt.dist(k) ) then !current sample closer than previous samples update the smallest distance
              dist(k)=diff2
           end if
        end do
        dist(k)=SQRT(dist(k))
        distmean= distmean+dist(k)
        RMSEmean=RMSEmean+RMSE(k)
        if(randomtestl.eq.2)      SIGMAmean=SIGMAmean+SIGMA(k)
     end do

     distmean=distmean/real(NTOEX)
     RMSEmean=RMSEmean/real(NTOEX)
     if(randomtestl.eq.2)   SIGMAmean=SIGMAmean/real(NTOEX)

     ! Pick test candidate with largest difference in values, but above distcomp distance to nearest neighbours

     !     if (iterDel.eq.1) then
     !        distcomp=1.5*distmean
     !     else
     distcomp=distmean !0.618 ! 
     !     end if


     diffloctmp=0.0d0
     diffloc2=0.0d0
     do k=1,NTOEX
        diffloctmp=(maxftoex(k)-minftoex(k))**2
        diffloc2=diffloc2+diffloctmp
     end do
     diffloc2=diffloc2/dble(ntoex)
     diffloc2=sqrt(diffloc2)

     do ii=1,nptstoaddpercyc

        npass=0
        do while (npass.ne.1)
           kp=0

           diffloctmp=0.0d0

           do k=1,NTOEX
              !! Successful Training point passes the follwing tests:

              !! 1. The local difference between the local and global surrogate should > 
              !! 2. The next training point should be atleast distcomp away from the closest existing sample
              !! 3. RMSE should be greater than RMSE mean of Kriging
              !! 4. SIGMA should be greater than SIGMA mean of MIR 
              !
              if ((maxftoex(k)-minftoex(k)).gt.diffloctmp .and. dist(k).ge.distcomp) then !.and. SIGMA(k).ge.SIGMAmean .and. RMSE(k).ge.RMSEmean
                 !                 print*,RMSE(k),RMSEmean
                 !     if (RMSE(k).ge.RMSEmean) then
                 diffloctmp=maxftoex(k)-minftoex(k)
                 kp=k
              end if
           end do

           diffloc=max(diffloc,diffloctmp)

           if (kp.eq.0) then ! if no successful candidate

              write (filenum,*) '  >>No passing candidate found . . .'
              write (filenum,*) '  >>Relaxing geometric constraint by 2 percent . . .'
              distcomp=0.98*distcomp

!!$           write (filenum,*) 'Could not find suitable test candidate just take the one with largest difference'
!!$           diffloctmp=0.0
!!$           do k=1,NTOEX
!!$              if ((maxftoex(k)-minftoex(k)).gt.diffloctmp) then
!!$                 diffloctmp=maxftoex(k)-minftoex(k)
!!$                 kp=k
!!$              end if
!!$           end do

           else 

              npass=npass+1 ! one successful candidate

              write(filenum,*)
              write(filenum,*) '>>Loc diff is',diffloctmp,' for candidate',ii,' at iteration',iterDEL,kp
              write(filenum,*)

           end if

        end do! while loop execute until a passing candidate is found


        ! Trick to not consider this point again
        maxftoex(kp)=minftoex(kp)
        RMSE(kp)=-10000.0

        ! Update other minimum distances
        do k=1,NTOEX
           if (k.ne.kp) then
              diff2=0.0
              do jj=1,ndim
                 diff2=diff2+(Dtoex(jj,k)-Dtoex(jj,kp))**2
              end do
              diff2=SQRT(diff2)
              if ( diff2.lt.dist(k) ) then
                 dist(k)=diff2
              end if
           end if
        end do

        !   call stop_all

        hstatad(ii)=hstat
        if (selectedevaluation.eq.1) then
           ! If local difference of point to add is smaller than average of past local differences we only want to calculate function values (and maybe gradient values)
           difflocmin=1000000.0
           difflocavg=0.0
           do i=1,ndiffloc
              difflocavg=difflocavg+difflocar(i)
              if (difflocar(i).lt.difflocmin) difflocmin=difflocar(i)
           end do
           difflocavg=difflocavg/ndiffloc

           print *,'Local difference, Avg, Min',diffloctmp,difflocavg,difflocmin

           if (ndiffloc.ge.2 .and. hstat.gt.0) then
              if (diffloctmp.le.difflocavg) hstatad(ii)=1
              if (diffloctmp.le.difflocmin) hstatad(ii)=0
           end if

           ndiffloc=ndiffloc+1
           difflocar(ndiffloc)=diffloctmp

        end if !selected evaluation

        ! Evaluate desired quantities
        if(nstyle.eq.0) then
           Dad(:,ii)=Dtoex(:,kp)
           call evalfunc(Dtoex(:,kp),ndim,fct,0,hstatad(ii),f(ii),df(:,ii),d2f(:,:,ii),v(:,ii))
        end if

     end do !! ii loop

!!$  
!!$  if(kp.eq.0) then
!!$     nptstoaddpercyc=kpnonzeroiter
!!$     print *,kp, ii,nptstoaddpercyc,  kpnonzeroiter
!!$  else
!!$
!!$     open(110,file='npoints.dat',form='formatted',status='unknown')
!!$     read(110,*) nptstoaddpercyc
!!$     close(110)
!!$  end if
!!$

     !     diffloc=diffloc2 ! trick to write the mean diff instead of maxdiff


     ! Add successful test candidates to sample points 
     if(nstyle.eq.0)then
        open(10,file='sample.dat',form='formatted',status='unknown')
        !$$ write(10,'(3i8)')ndim,nhs+nls+nptstoaddpercyc,3
        write(10,'(3i8)')ndim,nhs+nls+nptstoaddpercyc,2
        do i=1,nhs+nls
           if (info(i)(3:6).ne.'FGHv') then
              write(10,101) info(i),(sample(i,j),j=1,ndim),(func(i,j),j=1,nfunc),((gfunc(i,nfCOK(j),k),k=1,ndim),j=1,nCOK),(((hfunc(i,nfCOK(j),k,l),l=1,ndim),k=1,ndim),j=1,nCOK)
           else
              write(10,101) info(i),(sample(i,j),j=1,ndim),(func(i,j),j=1,nfunc),((gfunc(i,nfCOK(j),k),k=1,ndim),j=1,nCOK),((hvect(i,nfCOK(j),k,1),k=1,ndim),j=1,nCOK),((hvect(i,nfCOK(j),k,2),k=1,ndim),j=1,nCOK)
           end if
        end do
        do ii=1,nptstoaddpercyc
           if (hstatad(ii).le.3) then
              !$$ write(10,102) hstatad(ii),(Dad(j,ii),j=1,ndim),f(ii),f(ii),0.d0,(df(j,ii),j=1,ndim),((d2f(k,j,ii),j=1,ndim),k=1,ndim)
              write(10,102) hstatad(ii),(Dad(j,ii),j=1,ndim),f(ii),0.d0,(df(j,ii),j=1,ndim),((d2f(k,j,ii),j=1,ndim),k=1,ndim)
           else
              !$$ write(10,102) hstatad(ii),(Dad(j,ii),j=1,ndim),f(ii),f(ii),0.d0,(df(j,ii),j=1,ndim),(v(j,ii),j=1,ndim),(d2f(k,1,ii),k=1,ndim)
              write(10,102) hstatad(ii),(Dad(j,ii),j=1,ndim),f(ii),0.d0,(df(j,ii),j=1,ndim),(v(j,ii),j=1,ndim),(d2f(k,1,ii),k=1,ndim)
           end if
        end do
        close(10)

     end if !nstyle

  END IF! master thread 


  call MPI_Barrier(MPI_COMM_WORLD,ierr)!slaves wait until master arrives

  call deallocate_all_krig

101 format(a,10000e20.10)
102 format(i1,10000e20.10)

end subroutine DynamicPointSelection
! call knn(Dtoex(:,k),sample,knnptr,ndim,nhs,NCP)

subroutine mirtunableparams(fct,ndim,nhs,ncp,taylororder)
  use dimKrig,only:ndimt
  implicit none
  integer,INTENT(IN)::fct,ndim,nhs
  INTEGER,INTENT(OUT)::NCP,TAYLORORDER

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !                 EXP
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  if (Fct.eq.4) then
     if (nhs.le.11)  then  
        NCP=nhs     
        Taylororder=3!nhs!2!2!INT(NHS/4)

     else if (nhs.gt.11 .and. nhs.le.15)  then  
        NCP=10
        Taylororder=5!ncp

     else  if (nhs.gt.15 .and. nhs.le.25)  then  
        NCP=15
        Taylororder=5!ncp

     else  if (nhs.gt.25 .and. nhs.le.35)  then  
        NCP=20
        Taylororder=5!ncp

     else if (nhs.gt.35 .and. nhs.le.100)  then  
        NCP=20!20+3*ndim
        Taylororder=5
        ! else if (nhs.gt.70 .and. nhs.le.100)  then  
        !   NCP=30!50+0.1*nhs
        !    Taylororder=17       
     else
        NCP=20!+5*ndim
        tAYLORORDER=5
     end if

  end if
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!    
  !                COS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  if (Fct.eq.1) then
     if (nhs.le.11)  then  
        NCP=nhs
        Taylororder=ncp!2!2!INT(NHS/4)
     else if (nhs.gt.11 .and. nhs.le.15)  then  
        NCP=10
        Taylororder=5!ncp

     else if (nhs.gt.15 .and. nhs.le.25)  then  
        NCP=15
        Taylororder=5!ncp

     else if (nhs.gt.25 .and. nhs.le.35)  then  
        NCP=20
        Taylororder=5!ncp

     else if (nhs.gt.35 .and. nhs.le.70)  then  
        NCP=20!50+0.1*nhs
        Taylororder=5
     else if (nhs.gt.70 .and. nhs.le.100)  then  
        NCP=20!50+0.1*nhs
        Taylororder=5       
     else
        NCP=20
        tAYLORORDER=5
     end if

  end if
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !               RUNGE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!     
  if (fct.eq.2) then
     if (nhs.le.20)  then  
        NCP=nhs
        Taylororder=5!nhs!2!2!INT(NHS/4)

        !   else if (nhs.gt.11 .and. nhs.le.15)  then  
        !      NCP=nhs
        !      Taylororder=7!ncp
        !      
        !   else if (nhs.gt.15 .and. nhs.le.25)  then  
        !      NCP=nhs
        !      Taylororder=13!ncp
        !      
        !   else if (nhs.gt.25 .and. nhs.le.35)  then  
        !      NCP=nhs
        !      Taylororder=17!ncp
        !      
        !   else if (nhs.gt.35 .and. nhs.le.70)  then  
        !      NCP=30!50+0.1*nhs
        !      Taylororder=17
     else  if (nhs.gt.20 .and. nhs.le.100)  then  
        NCP=20!50+0.1*nhs
        Taylororder=5       
     else
        NCP=20
        tAYLORORDER=5
     end if

  end if
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !           ROSENBROCK    
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  if (Fct.eq.6) then
     if (nhs.le.11)  then  
        NCP=5
        Taylororder=ncp!2!2!INT(NHS/4)
     else if (nhs.gt.11 .and. nhs.le.15)  then  
        NCP=10
        Taylororder=5!ncp

     else if (nhs.gt.15 .and. nhs.le.25)  then  
        NCP=15
        Taylororder=5!ncp

     else if (nhs.gt.25 .and. nhs.le.35)  then  
        NCP=25
        Taylororder=5!ncp

     else if (nhs.gt.35 .and. nhs.le.70)  then  
        NCP=35!50+0.1*nhs
        Taylororder=5
     else if (nhs.gt.70 .and. nhs.le.100)  then  
        NCP=50!50+0.1*nhs
        Taylororder=5       
     else
        NCP=70
        tAYLORORDER=5
     end if

  end if
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !                      CFD   
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  if (fct.eq.20) then
     if (nhs.le.10)  then  
        NCP=nhs
        Taylororder=nhs!nhs!2!2!INT(NHS/4)
     else if (nhs.gt.10 .and. nhs.le.15)  then  
        NCP=nhs
        Taylororder=10!ncp

     else  if (nhs.gt.15 .and. nhs.le.25)  then  
        NCP=nhs
        Taylororder=10!ncp

     else  if (nhs.gt.25 .and. nhs.le.35)  then  
        NCP=nhs
        Taylororder=7!ncp

     else if (nhs.gt.35 .and. nhs.le.45)  then  
        NCP=30
        Taylororder=7

     else if (nhs.gt.45 .and. nhs.le.100)  then  
        NCP=30!50+0.1*nhs
        Taylororder=7
     else
        NCP=35
        Taylororder=7

     end if


  else


     if (nhs.le.11)  then  
        NCP=nhs
        Taylororder=ncp!2!2!INT(NHS/4)
     else if (nhs.gt.11 .and. nhs.le.15)  then  
        NCP=10
        Taylororder=5!ncp

     else if (nhs.gt.15 .and. nhs.le.25)  then  
        NCP=15
        Taylororder=5!ncp

     else if (nhs.gt.25 .and. nhs.le.35)  then  
        NCP=20
        Taylororder=5!ncp

     else if (nhs.gt.35 .and. nhs.le.70)  then  
        NCP=20!50+0.1*nhs
        Taylororder=5
     else if (nhs.gt.70 .and. nhs.le.100)  then  
        NCP=20!50+0.1*nhs
        Taylororder=5       
     else
        NCP=20
        tAYLORORDER=5
     end if

  end if ! end of CFD 

  if(ndimt.gt.2) then ! use 25% of the existing data points

     !NCP= ceiling(0.25*dble(nhs))

     if (Nhs.lt.50) then
        NCP=nhs
     else
        ncp=50
     end if
     !        print *,'NCP:',ncp
     tAYLORORDER=5
  end if

  return
end subroutine mirtunableparams
subroutine knn(SC,sample,knnptr,ndim,nhs,NCP)
  ! Subroutine to find NCP closest neighbours from array sample to point SC
  implicit none
  integer :: ndim,nhs,NCP,j,k,knnptr(NCP)
  real*8 ::SC(ndim),sample(nhs,ndim),norm2(nhs),sn

  ! Calculate all distances
  do j=1,nhs
     norm2(j)=0.0
     do k=1,ndim
        norm2(j)=norm2(j)+(sample(j,k)-SC(k))**2
     end do
  end do

  ! Find NCP closest neighbours

  do k=1,NCP
     sn=10000000000.0
     do j=1,nhs
        if (norm2(j).lt.sn) then
           sn=norm2(j)
           knnptr(k)=j
        end if
     end do
     norm2(knnptr(k))=10000000000.0
  end do

end subroutine knn


subroutine fixcolindelaunay(triangle_num,triangle_coor_num,ndim,triangle_coor,triangle_node)
  ! Subroutine to check whether delaunay triangles are colinear 
  implicit none
  integer :: triangle_num,triangle_coor_num,ndim,triangle_node(ndim+1,triangle_num),ii,jj,kk,info,ipvt(ndim),Toremove(100000),Numtoremove
  real*8 :: triangle_coor(ndim,triangle_coor_num),x0(ndim),Css(ndim,ndim)


  Numtoremove=0
  do ii=1,triangle_num

     x0(1:ndim)=triangle_coor(1:ndim,triangle_node(1,ii))

     Css(:,:)=0.0
     do jj=2,ndim+1
        do kk=1,ndim
           Css(kk,jj-1)=triangle_coor(kk,triangle_node(jj,ii))-x0(kk)
        end do
     end do

     call dgefa(Css,ndim,ndim,ipvt,info)

     if (info.ne.0) then
        ! Flag triangle
        Numtoremove=Numtoremove+1
        Toremove(Numtoremove)=ii
     end if

  end do

  ! Remove flagged triangles from list
  do ii=1,Numtoremove
     triangle_node(:,Toremove(ii):triangle_num-1)=triangle_node(:,Toremove(ii)+1:triangle_num)
     triangle_num=triangle_num-1
     Toremove(ii+1:Numtoremove)=Toremove(ii+1:Numtoremove)-1
  end do

end subroutine fixcolindelaunay





subroutine triangulation_order3_plot ( file_name, node_num, node_xy, &
     triangle_num, triangle_node, node_show, triangle_show )

  !*****************************************************************************80
  !
  !! TRIANGULATION_ORDER3_PLOT plots a 3-node triangulation of a set of nodes.
  !
  !  Discussion:
  !
  !    The triangulation is most usually a Delaunay triangulation,
  !    but this is not necessary.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 March 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) FILE_NAME, the name of the output file.
  !
  !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
  !
  !    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
  !
  !    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
  !
  !    Input, integer ( kind = 4 ) TRIANGLE_NODE(3,TRIANGLE_NUM), lists, for each 
  !    triangle, the indices of the nodes that form the vertices of the triangle.
  !
  !    Input, integer ( kind = 4 ) NODE_SHOW,
  !    0, do not show nodes;
  !    1, show nodes;
  !    2, show nodes and label them.
  !
  !    Input, integer ( kind = 4 ) TRIANGLE_SHOW,
  !    0, do not show triangles;
  !    1, show triangles;
  !    2, show triangles and label them.
  !
  !  Local parameters:
  !
  !    Local, integer ( kind = 4 ) CIRCLE_SIZE, controls the size of the circles 
  !    depicting the nodes.  Currently set to 5.  3 is pretty small, and 1 is
  !    barely visible.
  !
  implicit none

  integer   ( kind = 4 ) node_num
  integer   ( kind = 4 ) triangle_num

  real      ( kind = 8 ) ave_x
  real      ( kind = 8 ) ave_y
  character ( len = 40 ) date_time
  integer   ( kind = 4 ), parameter :: circle_size = 5
  integer   ( kind = 4 ) delta
  integer   ( kind = 4 ) e
  character ( len = * ) file_name
  integer   ( kind = 4 ) file_unit
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i4_wrap
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) node
  integer   ( kind = 4 ) node_show
  real      ( kind = 8 ) node_xy(2,node_num)
  character ( len = 40 ) string
  integer   ( kind = 4 ) triangle
  integer   ( kind = 4 ) triangle_node(3,triangle_num)
  integer   ( kind = 4 ) triangle_show
  real      ( kind = 8 ) x_max
  real      ( kind = 8 ) x_min
  integer   ( kind = 4 ) x_ps
  integer   ( kind = 4 ) :: x_ps_max = 576
  integer   ( kind = 4 ) :: x_ps_max_clip = 594
  integer   ( kind = 4 ) :: x_ps_min = 36
  integer   ( kind = 4 ) :: x_ps_min_clip = 18
  real      ( kind = 8 ) x_scale
  real      ( kind = 8 ) y_max
  real      ( kind = 8 ) y_min
  integer   ( kind = 4 ) y_ps
  integer   ( kind = 4 ) :: y_ps_max = 666
  integer   ( kind = 4 ) :: y_ps_max_clip = 684
  integer   ( kind = 4 ) :: y_ps_min = 126
  integer   ( kind = 4 ) :: y_ps_min_clip = 108
  real      ( kind = 8 ) y_scale
  !
  !  We need to do some figuring here, so that we can determine
  !  the range of the data, and hence the height and width
  !  of the piece of paper.
  !
  x_max = maxval ( node_xy(1,1:node_num) )
  x_min = minval ( node_xy(1,1:node_num) )
  x_scale = x_max - x_min

  x_max = x_max + 0.05D+00 * x_scale
  x_min = x_min - 0.05D+00 * x_scale
  x_scale = x_max - x_min

  y_max = maxval ( node_xy(2,1:node_num) )
  y_min = minval ( node_xy(2,1:node_num) )
  y_scale = y_max - y_min

  y_max = y_max + 0.05D+00 * y_scale
  y_min = y_min - 0.05D+00 * y_scale
  y_scale = y_max - y_min

  if ( x_scale < y_scale ) then

     delta = nint ( real ( x_ps_max - x_ps_min, kind = 8 ) &
          * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

     x_ps_max = x_ps_max - delta
     x_ps_min = x_ps_min + delta

     x_ps_max_clip = x_ps_max_clip - delta
     x_ps_min_clip = x_ps_min_clip + delta

     x_scale = y_scale

  else if ( y_scale < x_scale ) then

     delta = nint ( real ( y_ps_max - y_ps_min, kind = 8 ) &
          * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

     y_ps_max      = y_ps_max - delta
     y_ps_min      = y_ps_min + delta

     y_ps_max_clip = y_ps_max_clip - delta
     y_ps_min_clip = y_ps_min_clip + delta

     y_scale = x_scale

  end if

  file_unit=99

  open ( unit = file_unit, file = file_name, status = 'replace', &
       iostat = ios )

  if ( ios /= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'TRIANGULATION_ORDER3_PLOT - Fatal error!'
     write ( *, '(a)' ) '  Can not open output file "', trim ( file_name ), '".'
     return
  end if

  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: triangulation_order3_plot.f90'
  write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( file_unit, '(a)' ) '%%CreationDate: ' // trim ( date_time )
  write ( file_unit, '(a)' ) '%%Pages: 1'
  write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
       x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( file_unit, '(a)' ) '%%EndComments'
  write ( file_unit, '(a)' ) '%%BeginProlog'
  write ( file_unit, '(a)' ) '/inch {72 mul} def'
  write ( file_unit, '(a)' ) '%%EndProlog'
  write ( file_unit, '(a)' ) '%%Page: 1 1'
  write ( file_unit, '(a)' ) 'save'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
  write ( file_unit, '(a)' ) 'stroke'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB color to black.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the font and its size.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '/Times-Roman findfont'
  write ( file_unit, '(a)' ) '0.50 inch scalefont'
  write ( file_unit, '(a)' ) 'setfont'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Print a title.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  210  702  moveto'
  write ( file_unit, '(a)' ) '%  (Triangulation)  show'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
       x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
       x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
       x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
       x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
       x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a)' ) 'clip newpath'
  !
  !  Draw the nodes.
  !
  if ( 1 <= node_show ) then
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
     write ( file_unit, '(a)' ) '%'

     do node = 1, node_num

        x_ps = int ( &
             ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
             + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
             / ( x_max                   - x_min ) )

        y_ps = int ( &
             ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
             + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
             / ( y_max                   - y_min ) )

        write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
             circle_size, '0 360 arc closepath fill'

     end do

  end if
  !
  !  Label the nodes.
  !
  if ( 2 <= node_show ) then

     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '%  Label the nodes:'
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
     write ( file_unit, '(a)' ) '/Times-Roman findfont'
     write ( file_unit, '(a)' ) '0.20 inch scalefont'
     write ( file_unit, '(a)' ) 'setfont'
     write ( file_unit, '(a)' ) '%'

     do node = 1, node_num

        x_ps = int ( &
             ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
             + (       + node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
             / ( x_max                   - x_min ) )

        y_ps = int ( &
             ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
             + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
             / ( y_max                   - y_min ) )

        write ( string, '(i4)' ) node
        string = adjustl ( string )

        write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
             ' moveto (' // trim ( string ) // ') show'

     end do

  end if
  !
  !  Draw the triangles.
  !
  if ( 1 <= triangle_show ) then
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '%  Draw the triangles.'
     write ( file_unit, '(a)' ) '%'

     do triangle = 1, triangle_num

        write ( file_unit, '(a)' ) 'newpath'

        do i = 1, 4

           e = i4_wrap ( i, 1, 3 )

           node = triangle_node(e,triangle)

           x_ps = int ( &
                ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
                + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
                / ( x_max                   - x_min ) )

           y_ps = int ( &
                ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
                + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
                / ( y_max                   - y_min ) )

           if ( i == 1 ) then
              write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' moveto'
           else
              write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'
           end if

        end do

        write ( file_unit, '(a)' ) 'stroke'

     end do

  end if
  !
  !  Label the triangles.
  !
  if ( 2 <= triangle_show ) then

     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '%  Label the triangles:'
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '%  Set the RGB color to darker red.'
     write ( file_unit, '(a)' ) '%'
     write ( file_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
     write ( file_unit, '(a)' ) '/Times-Roman findfont'
     write ( file_unit, '(a)' ) '0.20 inch scalefont'
     write ( file_unit, '(a)' ) 'setfont'
     write ( file_unit, '(a)' ) '%'

     do triangle = 1, triangle_num

        ave_x = 0.0D+00
        ave_y = 0.0D+00

        do i = 1, 3

           node = triangle_node(i,triangle)

           ave_x = ave_x + node_xy(1,node)
           ave_y = ave_y + node_xy(2,node)

        end do

        ave_x = ave_x / 3.0D+00
        ave_y = ave_y / 3.0D+00

        x_ps = int ( &
             ( ( x_max - ave_x         ) * real ( x_ps_min, kind = 8 )   &
             + (       + ave_x - x_min ) * real ( x_ps_max, kind = 8 ) ) &
             / ( x_max         - x_min ) )

        y_ps = int ( &
             ( ( y_max - ave_y         ) * real ( y_ps_min, kind = 8 )   &
             + (         ave_y - y_min ) * real ( y_ps_max, kind = 8 ) ) &
             / ( y_max         - y_min ) )

        write ( string, '(i4)' ) triangle
        string = adjustl ( string )

        write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps, ' moveto (' &
             // trim ( string ) // ') show'

     end do

  end if

  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'restore  showpage'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  End of page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%%Trailer'
  write ( file_unit, '(a)' ) '%%EOF'
  close ( unit = file_unit )

  return
end subroutine triangulation_order3_plot


function i4_wrap ( ival, ilo, ihi )

  !*****************************************************************************80
  !
  !! I4_WRAP forces an integer to lie between given limits by wrapping.
  !
  !  Example:
  !
  !    ILO = 4, IHI = 8
  !
  !    I  I4_WRAP
  !
  !    -2     8
  !    -1     4
  !     0     5
  !     1     6
  !     2     7
  !     3     8
  !     4     4
  !     5     5
  !     6     6
  !     7     7
  !     8     8
  !     9     4
  !    10     5
  !    11     6
  !    12     7
  !    13     8
  !    14     4
  !
  !  Modified:
  !
  !    15 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer IVAL, an integer value.
  !
  !    Input, integer ILO, IHI, the desired bounds for the integer value.
  !
  !    Output, integer I4_WRAP, a "wrapped" version of IVAL.
  !
  integer i4_modp
  integer i4_wrap
  integer ihi
  integer ilo
  integer ival
  integer wide

  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
     i4_wrap = ilo
  else
     i4_wrap = ilo + i4_modp ( ival-ilo, wide )
  end if

  return
end function i4_wrap

function i4_modp ( i, j )

  !*****************************************************************************80
  !
  !! I4_MODP returns the nonnegative remainder of integer division.
  !
  !  Formula:
  !
  !    If
  !      NREM = I4_MODP ( I, J )
  !      NMULT = ( I - NREM ) / J
  !    then
  !      I = J * NMULT + NREM
  !    where NREM is always nonnegative.
  !
  !  Comments:
  !
  !    The MOD function computes a result with the same sign as the
  !    quantity being divided.  Thus, suppose you had an angle A,
  !    and you wanted to ensure that it was between 0 and 360.
  !    Then mod(A,360) would do, if A was positive, but if A
  !    was negative, your result would be between -360 and 0.
  !
  !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
  !
  !  Examples:
  !
  !        I     J     MOD  I4_MODP    Factorization
  !
  !      107    50       7       7    107 =  2 *  50 + 7
  !      107   -50       7       7    107 = -2 * -50 + 7
  !     -107    50      -7      43   -107 = -3 *  50 + 43
  !     -107   -50      -7      43   -107 =  3 * -50 + 43
  !
  !  Modified:
  !
  !    02 March 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer I, the number to be divided.
  !
  !    Input, integer J, the number that divides I.
  !
  !    Output, integer I4_MODP, the nonnegative remainder when I is
  !    divided by J.
  !
  integer i
  integer i4_modp
  integer j

  if ( j == 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'I4_MODP - Fatal error!'
     write ( *, '(a,i6)' ) '  I4_MODP ( I, J ) called with J = ', j
     stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
     i4_modp = i4_modp + abs ( j )
  end if

  return
end function i4_modp

