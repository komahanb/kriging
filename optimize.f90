      subroutine optimize(ndvar,D,ndvart,fobj,dfdD,low,up,gtol,maximize,outputscreen,fct)

        implicit none

        integer   ::   mmax,ndvar,ndvart,Nouter,maxfev,fct
        parameter      (mmax=80)  !mmax is the maximum number of limited memory corrections.

        double precision :: D(ndvart),fobj,dfdD(ndvart),dfdDtmp(ndvar),dfdDD(ndvart,ndvart),v(ndvart),gtol,low(ndvar),up(ndvar)
 
        character*60 :: csave,task
        logical ::      lsave(4),maximize,outputscreen
        integer ::      n, mlim, iprint, iwa(3*ndvar), isave(44), nbd(ndvar)
        double precision :: factr,pgtol,dsave(29),normgrad,dnrm2
        double precision :: wa(2*mmax*ndvar+4*ndvar+12*mmax*mmax+12*mmax)
        

        maxfev=100
        Nouter=100

!     We now provide nbd which defines the bounds on the variables:
!       nbd(i)=0 if x(i) is unbounded,
!              1 if x(i) has only a lower bound,
!              2 if x(i) has both lower and upper bounds, 
!              3 if x(i) has only an upper bound. 

        do n=1,ndvar
           nbd(n)=2
        end do

!     We suppress the default output.

        iprint = -1
     
!     We suppress both code-supplied stopping tests because the
!        user is providing his own stopping criteria.

        factr=0.0d0
        pgtol=0.0d0

!     We specify the number m of limited memory corrections stored.  
        
        mlim=mmax

!     We start the iteration by initializing task.
! 
        task = 'START'
        
!        ------- the beginning of the loop ----------
 
111     continue

      
!     This is the call to the L-BFGS-B code.
 
        call setulb(ndvar,mlim,D,low,up,nbd,fobj,dfdDtmp,factr,pgtol,wa,iwa,task,iprint,csave,lsave,isave,dsave)
 

        if (task(1:2) .eq. 'FG') then
!        the minimization routine has returned to request the
!        function f and gradient dfdDtmp values at the current D.


           if(fct.ge.1 .and. fct.le.4) then
         
              call CalcstuffBFGS(D,ndvart,fobj,dfdD,fct)

           else

              call Eulersolve(D,ndvart,0,fobj,dfdD,dfdDD,1,v,fct-10)
              
           end if

           dfdDtmp(1:ndvar)=dfdD(1:ndvar)

           if (maximize) then
              fobj=-fobj
              dfdDtmp(:)=-dfdDtmp(:)
           end if

         
           if (isave(34).eq.0) then
              
              normgrad = dnrm2(ndvar, dfdDtmp, 1)
              
              if (normgrad .le. gtol) task='STOP: THE GRADIENT IS SUFFICIENTLY SMALL'
              
              if(outputscreen) write (*,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iter',isave(30),'nfg =',1,'f =',fobj,'|grad| =',normgrad
              
           end if
           
           !        go back to the minimization routine.
           goto 111
        endif
        

        if (task(1:5) .eq. 'NEW_X') then   
!     
!        the minimization routine has returned with a new iterate.
!        At this point have the opportunity of stopping the iteration 
!        or observing the values of certain parameters

!        Note: task(1:4) must be assigned the value 'STOP' to terminate  
!          the iteration and ensure that the final results are
!          printed in the default format. The rest of the character
!          string TASK may be used to store other information.

           if (isave(34) .gt. maxfev) task='STOP: TOTAL NO. of f EVALUATIONS EXCEEDS LIMIT'
           
           if (isave(30) .ge. Nouter) task='STOP: TOTAL NO. OF OUTER ITERATIONS EXCEEDS LIMIT'

           normgrad = dnrm2(ndvar, dfdDtmp, 1)
           
           if (normgrad .le. gtol) task='STOP: THE GRADIENT IS SUFFICIENTLY SMALL'

!        We now wish to print the following information at each
!        iteration:
!        
!          1) the current iteration number, isave(30),
!          2) the total number of f and g evaluations, isave(34),
!          3) the value of the objective function fobj,
!          4) the norm of the projected gradient,  dsave(13)
         
           if(outputscreen) write (*,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iter',isave(30),'nfg =',isave(34),'f =',fobj,'|grad| =',normgrad

 
           !        go back to the minimization routine.
           goto 111
           
        endif

        if (maximize) then
           fobj=-fobj
        end if

        if(outputscreen) write(*,*) task

        if (task(1:5) .ne. 'CONVE') then
           write(*,*)
           write(*,*) 'Problem in epistemic optimization!'
           write(*,*)
           write(*,*) task
           write(*,*)
           stop
        end if

      end subroutine optimize



