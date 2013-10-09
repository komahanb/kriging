        subroutine Krigingestimate(ndimin,xavgin,ndimint,xstdin,fmeanout,fvarout,fmeanprimeout,fvarprimeout,NMCin,fctindxin)
        use dimKrig
        implicit none
        include 'mpif.h'
        integer :: ierr,i,j,imax,nstattmp,ndimtmp,NCP,lenc,ndimin,NMCin,fctindxin,ndimint,Casemode
        double precision :: diffconv,xavgin(ndimint),xstdin(ndimint),fmeanout,fvarout,fmeanprimeout(ndimint),fvarprimeout(ndimint)
        character*60 :: filename
        character*2 :: dimnumber,fctnumber
        character*3 :: lsnumber
        double precision :: freal,df(20),d2f(20,20),v(20),Javg,Jvar

        ! Dimension of problem
     
        ndim=ndimin
        ndimt=ndimint
        
        xavgt(:)=xavgin(:)
        xstdt(:)=xstdin(:)
       
        NMCS=NMCin

        fctindx=fctindxin

        filenum=77 ! 6 for screen

        ! Function choice and extrapolation order
        
        do fct=12,12     ! 1: cos(lin sum) 2: Runge fct 3: Rosenbrock fct 4: Rastrigin 5: Lin (cos plus noise)  6: Trustdesign 7: Quadratic 8: Cubic 9: Heaviside (Franke) 10: Eulersolve 11: Trustdesign maximization 12: Airfoil maximization/minimization
    do nstattmp=0,0    ! 0: f only  1: f+g  2: f+g+h  3: f+g+hv

        if (nstattmp.eq.0) hstat=0
        if (nstattmp.eq.1) hstat=1
        if (nstattmp.eq.2) hstat=3
        if (nstattmp.eq.3) hstat=5

        ! ---------------------------------------------------------------------------
        ! General parameter
        ! ---------------------------------------------------------------------------

        Casemode=1    ! 0: Calculate RMS  1: output stats  2: Optimize using EI

        dynamics=1       ! 0: not dynamic  1: dynamic 

        reusesamples=0   ! 0: no  1: yes

        Dutchorderg=1    ! Order of Interpolation for Dutch Intrapolation (1 or 2)

        randomini=1      ! 0: Corners of cube 1: latin hypercube

        ! ---------------------------------------------------------------------------

        ! Number of initial points

        if (randomini.eq.0) then
           nhs=1+2**ndim
        else
           call combination(ndim+Dutchorderg,ndim,nhs)
        end if 

        nls=0       ! number of low-fidelity samples
        lstat=0     ! 0: f only  1: f+g  3: f+g+h   
        ifid=2      ! fidelity level of function


        ! Maximum number of high-fidelity points and how many points to add per cycle

        maxsamplewant=13
        nptstoaddpercyc=ndim
        

        ! Dynamic sample point location parameters

        diffconv=1.e-6        

        randomtestl=1    ! 0: Delaunay triangulation with Dutchintrapolation 1: latin hypercube with Dutchintrapolation 2: latin hypercube with Multivariate Interpolation and Regression (MIR)

        selectedevaluation=0  ! 0: not selected 1: selected 

        ! How to get MC samples
               
        readMCsamples=0     ! 0: Make random samples   1: Read from file MCsamp.dat
       
        ! Get standard deviation in function space

        if (Casemode.eq.1) then
           if (readMcsamples.eq.1) then
              open(10,file='MCsamp.dat',form='formatted',status='unknown')
              read(10,*) NMCS,ndimtmp
              if (ndimtmp.ne.ndim) STOP 'Dimension does not match in MCsamp.dat!'
              read(10,*) (xavg(i),i=1,ndim)
              read(10,*) (xstd(i),i=1,ndim)              
!!$              ! Calculate mean and variance of real function
!!$              Javg=0.0
!!$              Jvar=0.0
!!$              DS(1,:)=0.0
!!$              DS(2,:)=1.0              
!!$              do j=1,NMCS
!!$                 read(10,*) (xavg(i),i=1,ndim)
!!$                 call evalfunc(xavg,ndim,fct,0,0,freal,df,d2f,v)
!!$                 Javg=Javg+freal
!!$                 Jvar=Jvar+freal**2
!!$              end do
!!$              Javg=Javg/real(NMCS)      
!!$              Jvar=Jvar/real(NMCS)-Javg**2
!!$              write(*,*)
!!$              write(*,'(6x,a,2e15.5)')'Real: Mean and Variance',Javg,Jvar
!!$              write(*,*)
!!$              stop
              close(10)
           else if (readMcsamples.eq.0 .and. NMCS.eq.0) then
              open(10,file='MC.inp',form='formatted',status='unknown')
              read(10,*) (xavg(i),i=1,ndim)
              read(10,*) (xstd(i),i=1,ndim)
              close(10)
           else if (readMcsamples.eq.0 .and. NMCS.ne.0) then
              xavg(1:ndim)=xavgin(ndimt-ndim+1:ndimt)
              xstd(1:ndim)=xstdin(ndimt-ndim+1:ndimt)
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
              else if (fct.ge.10) then
                 DS(1,i)=-2.0*0.01
                 DS(2,i)=2.0*0.01
!!$                 if (i.eq.1) then
!!$                    DS(1,i)=0.0*4.0*atan(1.0)/180.0   ! in radian
!!$                    DS(2,i)=4.0*4.0*atan(1.0)/180.0   ! in radian
!!$                 else
!!$                    DS(1,i)=0.5
!!$                    DS(2,i)=1.5
!!$                 end if                 
!!$                 if (i.eq.1) then
!!$                    DS(1,i)=-3.517039000000000E-002
!!$                    DS(2,i)=4.062709000000000E-002
!!$                 else
!!$                    DS(1,i)=-4.033711600000000E-002
!!$                    DS(2,i)=3.797685900000000E-002
!!$                 end if
              else
                 DS(1,i)=-2.0
                 DS(2,i)=2.0
              end if
           end if
        end do
          
        nstyle=0    ! With Function(0) or Separately(1)
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

        if (id_proc.eq.0) then

           if (dynamics.ne.0) call Make_Sample 
           call i_to_s(ndim,dimnumber)
           call i_to_s(fct,fctnumber)
           call i_to_s(nls,lsnumber)

           if (Casemode.eq.0) then

              filename='errornormdim'
              lenc=12
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
                 filename(lenc+1:lenc+3)='dyn'
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
 
              open(unit=93,file=filename)
              write(93,*) 'Nhs','    L2diff','          Maxdiff','       Dutchdiff'

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
                 filename(lenc+1:lenc+3)='dyn'
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

              
              open(unit=94,file=filename)
              if (hstat.eq.0) write(94,*) 'Nhs','      RealAVG','        RealVAR', '         KrigAVG', '         KrigVAR'
              if (hstat.eq.1) write(94,'(a131)') 'Nhs      RealAVG        RealVAR         KrigAVG         KrigVAR           MM1AVG         MM1VAR          LinAVG           LinVAR'
              if (hstat.eq.3) write(94,'(a196)') 'Nhs      RealAVG        RealVAR         KrigAVG         KrigVAR           MM1AVG         MM1VAR          LinAVG           LinVAR          MM2AVG         MM2VAR          QuadAVG        QuadVAR'

           end if  !Casemode.eq.0?

        end if

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nhs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(maxsample,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)       

        diffloc=1000.0
        ndiffloc=0

        iterDEL=0

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

              call MPI_Barrier(MPI_COMM_WORLD,ierr)
              call MPI_BCAST(Cmode,20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)                        
              if(Cmode(:15).eq.'Make_and_Search')then
                 Cmode='Make_by_GA'
                 call Make_Krig
                 Cmode='Search_by_GA2_Local'
                 call Search_Krig

              else if(Cmode(:5).eq.'Dynam' .and. id_proc.eq.0)then
                 call DynamicPointSelection

              else if(Cmode(:5).eq.'Make_')then
                 call Make_Krig
                 
              else if(Cmode(:7).eq.'Search_')then
                 call Search_Krig
                 
              else if(Cmode(:5).eq.'Post_')then
                 call Post_Process
                 
              else if(Cmode(:5).eq.'Rank ')then
                 call Rank_New_One
                 
              else if(Cmode(:7).eq.'Vrange ')then
                 call Variable_Range
              
              else if(Cmode(:12).eq.'Trust_Region')then
                 call Check_Trust_Region
                 
              else if(Cmode(:7).eq.'Update ')then

	         if (id_proc.eq.0) call Update
                 
                 if (fct.lt.20 .and. ndim.eq.2) then
                    Cmode='Post_2D'
                    call Post_Process
		 else 
		    call find_Optimal
                 end if
                                  
              else if(Cmode(:16).eq.'Trust_and_Update')then
                 call Check_Trust_Region
                 call Update
              end if

              if (i.eq.3) then
                 call MPI_Barrier(MPI_COMM_WORLD,ierr)
                 call MPI_BCAST(diffloc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
              end if

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
           if (Casemode.eq.0) close(93)
           if (Casemode.eq.1) close(94)
        end if


   end do
end do

      fmeanout=fmean
      fvarout=fvar

      fmeanprimeout(ndimt-ndim+1:ndimt)=fmeanprime(1:ndim)/(DS(2,1:ndim)-DS(1,1:ndim))
      fvarprimeout(ndimt-ndim+1:ndimt)=fvarprime(1:ndim)/(DS(2,1:ndim)-DS(1,1:ndim))

    end subroutine Krigingestimate




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
