program Kriging
  use dimKrig
  use timer_mod
  implicit none
  include 'mpif.h'
  integer :: ierr,i,j,imax,nstattmp,ndimtmp,NCP,lenc,NMCStmp,Casemode,mirdutchhybrid
  double precision :: diffconv
  character*60 :: filename
  character*2 :: dimnumber,fctnumber,nptstoaddpercycnum
  character*3 :: lsnumber
  double precision :: freal,df(20),d2f(20,20),v(20),Javg,Jvar!,distloc
  integer nsamples,counter
  integer,parameter::timing=1,numberpointstudy=0 !other number if timing results not needed
  double precision :: time, Initialmach, Finalmach, InitialAOA,FinalAOA
  integer::nmcstemp
  integer::ctest,fuct

  call MPI_START

  ! Dimension of problem

  ndim=2
  ndimt=ndim

  randomini=1      ! How intial samples are chosen? 0: Corners of cube 1: latin hypercube. If not dynamic it should be set to 1

  mainprog =.true.
  filenum=6


  do randomtestl=2,2   

     ! 0: Delaunay triangulation with Dutchintrapolation 
     ! 1: latin hypercube with Dutchintrapolation 
     ! 2: latin hypercube with Multivariate Interpolation and Regression (MIR)

     !     maxsamplewant= 90 !

     !     nptstoaddpercyc=100

     ! Low fidelity data

     nls=0      ! number of low-fidelity samples!0
     lstat=0     ! 0: f only  1: f+g  3: f+g+h   
     ifid=0     ! fidelity level of function!2

     ! Dynamic sample point location parameters

     diffconv=1.e-6      ! 

     selectedevaluation=0 ! 0: not selected 1: selected 

     filenum=6 ! 6 for screen

     NMCS=0

     ! CFD solve input data

     fctindx=0        !0 for drag, 4 for lift
     InitialMach=0.5d0
     FinalMach=1.5d0
     InitialAOA=0.0d0  !degrees
     FinalAOA=5.0d0    !degrees      
     ! ---------------------------------------------------------------------------
     ! General parameter
     ! ---------------------------------------------------------------------------

     do  Casemode= 0,0      ! 0: Calculate RMS  1: Output Statistics  2: Optimize using EI
        reusesamples=0   ! 0: no  1: yes 
        Dutchorderg=2    ! Order of Interpolation for Dutch Intrapolation (1 or 2)

        ! How to get MC samples
        readMCsamples=0     ! 0: Make random samples   1: Read from file MCsamp.dat

!!$
!!$        if (randomtestl.ne.0) then
!!$           dynamics=1
!!$        else 
!!$           dynamics=0
!!$        end if     ! 0: not dynamic  1: dynamic 

        do ctest=1,1
           
           if (ctest.eq.1) then

              dynamics=0
                           
              !1=LHS,
              !2=Nieder
              !3=Halton
              !4=Hammersley
              !5=Sobol
              !6=Faure

              randomflag=6
              
           end if
           if (ctest.eq.2) then
              dynamics=1
              lhsdyn=.true.
           end if
           if (ctest.eq.3) then
              dynamics=1
              lhsdyn=.false. !mirdyn is true
           end if

           ! ---------------------------------------------------------------------------

!!$
!!$         open(10,file='MCsamp.dat',form='formatted',status='unknown')
!!$         read(10,*) NMCStmp,ndimtmp
!!$         read(10,*) (xavg(i),i=1,ndim)
!!$         read(10,*) (xstd(i),i=1,ndim) 
!!$         close(1


        do fuct=2,2 !0:exp 1: cos(lin sum) 2: Runge fct 3: Rosenbrock fct 4: Rastrigin 5: Lin (cos plus noise)  6: Trustdesign 7: Quadratic 8: Cubic 9: Short Column, 10:  Cantilever, 11: Three Bar ,20: CFD, 21,22: Optimization

           if (fuct.eq.1) fct=0
           if (fuct.eq.2) fct=2
           if (fuct.eq.3) fct=3

           if (id_proc.eq.0) write(filenum,'(4x,a,i8)')">> Test case number",ctest
           if (id_proc.eq.0) write(filenum,'(4x,a,i8)')">> Test function number",fct

           evlfnc=1  ! CFD case exact evaluation for MC needed or not

           do nstattmp=0,0                ! 0: f only  1: f+g  2: f+g+h  3: f+g+hv
              if (nstattmp.eq.0) then

                 maxsamplewant= 150
                 nptstoaddpercyc=5!160

              else if (nstattmp.eq.1) then

                 maxsamplewant= 150 !1000+5
                 nptstoaddpercyc=5 

              else if (nstattmp.eq.2) then

                 maxsamplewant= 150!15!20+33
                 nptstoaddpercyc=5

              end if

              if (nstattmp.eq.0) hstat=0
              if (nstattmp.eq.1) hstat=1
              if (nstattmp.eq.2) hstat=3
              if (nstattmp.eq.3) hstat=5
              if (id_proc.eq.0) call TimerInit()

              counter=0
              do nsamples=5,5 !2*NDIM+1,2*NDIM+1!40,40!101,101!2*NDIM+1,2*NDIM+1!5,75,5!50,50!50,500,50!500,500!3,47,4 !Makes this many nhs samples per cycle
                 counter=counter+1

                 nhs=nsamples

!!$              if (randomini.eq.0) then
!!$                 nhs=1+2**ndim
!!$                 else
!!$                 call combination(ndim+Dutchorderg,ndim,nhs)
!!$              end if

                 ! Get standard deviation in function space

                 if (Casemode.eq.1) then
                    if (readMcsamples.eq.1) then

                       if (fct.eq.20) then
                          open(10,file='MC.inp',form='formatted',status='unknown')
                          read(10,*) (xavg(i),i=1,ndim)
                          read(10,*) (xstd(i),i=1,ndim)     
                          read(10,*)
                          read(10,*)
                          read(10,*) NMCS!,ndimtmp
                          close(10)
                       else 

                          open(10,file='MCsamp.dat',form='formatted',status='unknown')
                          read(10,*) NMCS,ndimtmp
                          if (ndimtmp.ne.ndim) STOP 'Dimension does not match in MCsamp.dat!'
                          read(10,*) (xavg(i),i=1,ndim)
                          read(10,*) (xstd(i),i=1,ndim)         
                          close(10)

                       end if
!!$                    ! Calculate mean and variance of real function
!!$                    Javg=0.0
!!$                    Jvar=0.0
!!$                    DS(1,:)=0.0
!!$                    DS(2,:)=1.0              
!!$                    do j=1,NMCS
!!$                       read(10,*) (xavg(i),i=1,ndim)
!!$                  !     call evalfunc(xavg,ndim,fct,0,0,freal,df,d2f,v)
!!$                       Javg=Javg+freal
!!$                       Jvar=Jvar+freal**2
!!$                    end do
!!$                    Javg=Javg/real(NMCS)      
!!$                    Jvar=Jvar/real(NMCS)-Javg**2
!!$                    write(*,*)
!!$                    write(*,'(6x,a,2e15.5)')'Real: Mean and Variance',Javg,Jvar
!!$                    write(*,*)
!!$                    !              stop
                       close(10)
                    else if (readMcsamples.eq.0 .and. NMCS.eq.0) then
                       open(10,file='MC.inp',form='formatted',status='unknown')
                       read(10,*) (xavg(i),i=1,ndim)
                       read(10,*) (xstd(i),i=1,ndim)
                       close(10)
                    end if
                 end if


                 ! Domain size in function space
                 do i=1,ndim
                    if (Casemode.eq.1) then
                       DS(1,i)=xavg(i)-2.0*xstd(i)
                       DS(2,i)=xavg(i)+2.0*xstd(i)
                    else
                       if (fct.eq.4) then
                          DS(1,i)=-5.12
                          DS(2,i)=5.12
                       else if (fct.eq.6) then
                          DS(1,i)=0.0d0
                          DS(2,i)=5.0d0
                       else if (fct.ge.20) then
                          !               DS(1,i)=-2.0*0.01
                          !               DS(2,i)=2.0*0.01
                          if (i.eq.1) then
                             DS(1,i)=InitialAOA*4.0d0*atan(1.0)/180.0   ! in radian
                             DS(2,i)=FinalAOA*4.0d0*atan(1.0)/180.0   ! in radian
                          else
                             DS(1,i)=InitialMach !initial mach number
                             DS(2,i)=FinalMach    !final mach number
                          end if
!!$                 if (i.eq.1) then
!!$                    DS(1,i)=-3.517039000000000E-002
!!$                    DS(2,i)=4.062709000000000E-002
!!$                 else
!!$                    DS(1,i)=-4.033711600000000E-002
!!$                    DS(2,i)=3.797685900000000E-002
!!$                 end if

                       else if (fct.eq.12) then

                          DS(1,i)=1.0
                          DS(2,i)=3.0
                       else
                          DS(1,i)=-2.0
                          DS(2,i)=2.0
                       end if
                    end if
                 end do

                 nstyle=0    				! With Function(0) or Separately(1)
                 if (dynamics.ne.0) then
                    if ((randomtestl.eq.0 .and. Dutchorderg.eq.2) .or. Dutchorderg.gt.2 .or. Dutchorderg.lt.1) then
                       write(*,*)
                       if (id_proc.eq.0) write(*,*) 'This global order is not supported (yet)'
                       call stop_all
                    end if
                    if (randomtestl.eq.0 .and. randomini.eq.1)  then
                       if (id_proc.eq.0) write(*,*) 'If dynamic with Delaunay triangulation randomini should be set to 0!'
                       call stop_all
                    end if
                 else if (dynamics.eq.0) then
                    if (randomini.ne.1)  then
                       if (id_proc.eq.0) write(*,*) 'If not dynamic randomini should be set to 1!'
                       call stop_all
                    end if
                 end if

                 maxsample=maxsamplewant

                 if (id_proc.eq.0.and.dynamics.ne.0) call Make_Sample 

                 if (id_proc.eq.0) then

                    call i_to_s(ndim,dimnumber)
                    call i_to_s(fct,fctnumber)
                    call i_to_s(nls,lsnumber)

                    if (numberpointstudy.eq.1) then
                       call i_to_s(nptstoaddpercyc,nptstoaddpercycnum)
                    end if

                    if (Casemode.eq.0) then
                       filename='KRIGerrornormdim'
                       lenc=16
                       filename(lenc+1:lenc+2)=dimnumber
                       lenc=lenc+2
                       filename(lenc+1:lenc+3)='fct'
                       lenc=lenc+3

                       filename(lenc+1:lenc+2)=fctnumber
                       lenc=lenc+2

                       if (numberpointstudy.eq.1) then
                          filename(lenc+1:lenc+3)='TPS'
                          lenc=lenc+3
                          filename(lenc+1:lenc+2)=nptstoaddpercycnum
                          lenc=lenc+2
                       end if

                       if (hstat.eq.5) then

                          filename(lenc+1:lenc+4)='FGHv'
                          lenc=lenc+4

                       else if (hstat.eq.3) then 

                          filename(lenc+1:lenc+3)='FGH'
                          lenc=lenc+3    

                       else if (hstat.eq.1) then

                          filename(lenc+1:lenc+2)='FG'
                          lenc=lenc+2 

                       else if (hstat.eq.0) then

                          filename(lenc+1:lenc+1)='F'
                          lenc=lenc+1 

                       else
                          print *, 'Wrong value in hstat'
                          call stop_all
                       end if

                       if (randomflag.eq.2) then

                          filename(lenc+1:lenc+1)='2'
                          lenc=lenc+1 

                       else if (randomflag.eq.3) then

                          filename(lenc+1:lenc+1)='3'
                          lenc=lenc+1 

                       else  if (randomflag.eq.4) then

                          filename(lenc+1:lenc+1)='4'
                          lenc=lenc+1 

                       else if (randomflag.eq.5) then

                          filename(lenc+1:lenc+1)='5'
                          lenc=lenc+1 

                       else if (randomflag.eq.6) then

                          filename(lenc+1:lenc+1)='6'
                          lenc=lenc+1 
                       end if


                       if (dynamics.eq.1) then
                          if (randomtestl.eq.1) filename(lenc+1:lenc+6)='dutdyn'
                          if (randomtestl.eq.2) filename(lenc+1:lenc+6)='mirdyn'
                          if (lhsdyn)  filename(lenc+1:lenc+6)='lhsdyn'
                          lenc=lenc+6
                       end if

                       if (nls.gt.0) then              
                          filename(lenc+1:lenc+3)='low'
                          lenc=lenc+3
                          filename(lenc+1:lenc+3)=lsnumber
                          lenc=lenc+3
                          if (lstat.eq.5) then
                             filename(lenc+1:lenc+4)='FGHv'
                             lenc=lenc+4
                          else if (lstat.eq.3) then 
                             filename(lenc+1:lenc+3)='FGH'
                             lenc=lenc+3    
                          else if (lstat.eq.1) then
                             filename(lenc+1:lenc+2)='FG'
                             lenc=lenc+2 
                          else if (lstat.eq.0) then
                             filename(lenc+1:lenc+1)='F'
                             lenc=lenc+1 
                          else
                             print *, 'Wrong value in lstat'
                             call stop_all
                          end if

                       end if

                       open(unit=93,file='norm/'//filename)!,position='append')
                       if (counter.eq.1) write(93,*) 'Nhs','    L2diff','          Maxdiff','       Dutchdiff' 

                    else if (Casemode.eq.1) then

                       filename='statsdim'
                       lenc=8
                       filename(lenc+1:lenc+2)=dimnumber
                       lenc=lenc+2
                       filename(lenc+1:lenc+3)='fct'
                       lenc=lenc+3
                       filename(lenc+1:lenc+2)=fctnumber
                       lenc=lenc+2

                       if (hstat.eq.5) then
                          filename(lenc+1:lenc+4)='FGHv'
                          lenc=lenc+4
                       else if (hstat.eq.3) then 
                          filename(lenc+1:lenc+3)='FGH'
                          lenc=lenc+3    
                       else if (hstat.eq.1) then
                          filename(lenc+1:lenc+2)='FG'
                          lenc=lenc+2 
                       else if (hstat.eq.0) then
                          filename(lenc+1:lenc+1)='F'
                          lenc=lenc+1 
                       else
                          print *, 'Wrong value in hstat'
                          call stop_all
                       end if

                       if (dynamics.eq.1) then
                          filename(lenc+1:lenc+6)='mirdyn'
                          lenc=lenc+3
                       end if

                       if (nls.gt.0) then              
                          filename(lenc+1:lenc+3)='low'
                          lenc=lenc+3
                          filename(lenc+1:lenc+3)=lsnumber
                          lenc=lenc+3
                          if (lstat.eq.5) then
                             filename(lenc+1:lenc+4)='FGHv'
                             lenc=lenc+4
                          else if (lstat.eq.3) then 
                             filename(lenc+1:lenc+3)='FGH'
                             lenc=lenc+3    
                          else if (lstat.eq.1) then
                             filename(lenc+1:lenc+2)='FG'
                             lenc=lenc+2 
                          else if (lstat.eq.0) then
                             filename(lenc+1:lenc+1)='F'
                             lenc=lenc+1 
                          else
                             print *, 'Wrong value in lstat'
                             call stop_all
                          end if

                       end if


                       open(unit=94,file='stats/'//filename)
                       if (hstat.eq.0) write(94,*) 'Nhs','      RealAVG','        RealVAR', '         KrigAVG', '    KrigVAR', '     ErrorAVG','       ErrorVAR'
                       if (hstat.eq.1) write(94,'(a131)') 'Nhs      RealAVG        RealVAR         KrigAVG         KrigVAR        ErrorAVG       ErrorVAR      MM1AVG         MM1VAR          LinAVG           LinVAR'
                       if (hstat.eq.3) write(94,'(a196)') 'Nhs      RealAVG        RealVAR         KrigAVG         KrigVAR        ErrorAVG       ErrorVAR   MM1AVG         MM1VAR          LinAVG           LinVAR          MM2AVG         MM2VAR          QuadAVG        QuadVAR'

                    end if

                 end if ! master thread

                 call MPI_Barrier(MPI_COMM_WORLD,ierr)
                 call MPI_BCAST(nhs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                 call MPI_BCAST(maxsample,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)       


                 diffloc=1000.0
                 ndiffloc=0

                 iterDEL=0

                 distloc=1.0

                 do while (nhs.le.maxsample .and. diffloc.gt.diffconv)

                    if (dynamics.eq.0 .and. id_proc.eq.0) call Make_Sample 
                    call MPI_Barrier(MPI_COMM_WORLD,ierr)
                    call MPI_BCAST(nhs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                    call MPI_BCAST(maxsample,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


                    iterDEL=iterDEL+1
                    call Read_Set
                    call Read_Sample
                    call Check_Sample
                    call Make_Mreg


                    if ((dynamics.eq.1 .and. nhs+nptstoaddpercyc.le.maxsample) .or. Casemode.eq.2) then
                       imax=3
                    else
                       imax=2
                    end if

                    do i=1,imax

                       if (i.eq.1) then
                          Cmode='Make_by_GA'
                       else if (i.eq.2) then
                          if (Casemode.eq.1) then
                             Cmode='Post_MonteCarlo'
                          else if (Casemode.eq.2) then
                             Cmode='Search_by_GA2_Local'
                          else if (Casemode.eq.0) then
                             if (ndim.eq.2) then
                                Cmode='Post_2D'
                             else if (ndim.gt.2) then
                                Cmode='Post_HigherD'
                             end if
                          end if
                       else if (i.eq.3) then
                          !Cmode='Search_by_MaxVar'
                          if (Casemode.le.1) then
                             Cmode='DynamicPointSelection'
                          else if (Casemode.eq.2) then
                             Cmode='Update'
                          end if
                       end if

                       call MPI_BCAST(Cmode,20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)          

                       if(Cmode(:15).eq.'Make_and_Search')then
                          Cmode='Make_by_GA'
                          call Make_Krig
                          Cmode='Search_by_GA2_Local'
                          call Search_Krig

                       else if(Cmode(:5).eq.'Dynam')then !.and. id_proc.eq.0
                          if (id_proc.eq.0) call TimerStart('Dyn Training Pt')        
                          call DynamicPointSelection
                          if (id_proc.eq.0) call TimerStop('Dyn Training Pt')

                       else if(Cmode(:5).eq.'Make_')then
                          call Make_Krig

                       else if(Cmode(:7).eq.'Search_')then
                          call Search_Krig

                       else if(Cmode(:5).eq.'Post_')then
                          if (id_proc.eq.0) call TimerStart('Post process')
                          call Post_Process
                          if (id_proc.eq.0) call TimerStop('Post process')

                       else if(Cmode(:5).eq.'Rank ')then
                          call Rank_New_One

                       else if(Cmode(:7).eq.'Vrange ')then
                          call Variable_Range

                       else if(Cmode(:12).eq.'Trust_Region')then
                          call Check_Trust_Region

                       else if(Cmode(:7).eq.'Update')then

                          if (id_proc.eq.0) call Update

                          if (fct.lt.20 .and. ndim.eq.2) then
                             Cmode='Post_2D'
                             if (id_proc.eq.0) call TimerStart('Post process')
                             call Post_Process
                             if (id_proc.eq.0) call TimerStop('Post process')


                             !open(90,file='KrigingSamples'nsamples, form='formatted',status='unknown')
                             !write(90,'(a)')'TITLE = " "'
                             !write(90,'(a)')'VARIABLES = "x" "y" "f"'
                             !write(90,'(3(a,i5),a)')'ZONE T="hsample", I=',ict,', J=',1,', K=',1,', F=BLOCK'
                             !write(90,'(9999f15.8)')(xvis(1,i),i=1,ict)
                             !write(90,'(9999f15.8)')(xvis(2,i),i=1,ict)
                             ! write(90,'(9999f15.8)')(fvis(i,1),i=1,ict)
                             !close(90)



                          end if

                       else if(Cmode(:16).eq.'Trust_and_Update')then
                          call Check_Trust_Region
                          call Update
                       end if

                       if (i.eq.3) call MPI_BCAST(diffloc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                       if (i.eq.3) call MPI_BCAST(distloc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

                    end do

                    if (Casemode.ne.2) then
                       nhs=nhs+nptstoaddpercyc
                    else
                       nhs=nhs+1
                    end if

                    ! Use Dutchorderg=2 once we exceed a certain amount of sample points
                    call combination(ndim+2,ndim,NCP)
                    if (randomtestl.eq.1 .and. nhs.gt.NCP) Dutchorderg=2

                    call Deallocating      

                 end do

                 if (id_proc.eq.0) then
                    ! if (Casemode.eq.0) close(93)
                    if (Casemode.eq.1) close(94)
                 end if

              end do

              if (id_proc.eq.0) then
                 if (Casemode.eq.0) close(93)
                 !           if (Casemode.eq.1) close(94)
              end if

              if (id_proc.eq.0.and.timing.eq.1) call TimerReport()

           end do ! F or G loop

        end do ! Function number loop
        
     end do ! test case loop
     
  end do !casemode loop
  
end do !random testl DI or MIR loop

!!$
!!$if (id_proc.eq.0) then
!!$print *, fmean,fmeanprime
!!$print *,fvar,fvarprime
!!$end if

  call stop_all

end program Kriging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine i_to_s(intval,s)

  integer idig,intval,ipos,ival,i
  character ( len = * ) s

  s = ' '

  ival = intval

  !  Working from right to left, strip off the digits of the integer
  !  and place them into S(1:len ( s )).
  !
  ipos = len(s) 

  do while ( ival /= 0 )

     idig = mod( ival, 10 )
     ival = ival / 10
     s(ipos:ipos) = char(idig + 48 )
     ipos = ipos - 1

  end do
  !
  !  Fill the empties with zeroes.
  !
  do i = 1, ipos
     s(i:i) = '0'
  end do

  return
end subroutine i_to_s

