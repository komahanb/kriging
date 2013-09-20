      subroutine BFGS(mode,id_proc,num_proc, &
                      ndim,tdim,iscf,ndebug,nevamax, &
                      bound, &
                      theta,power )
      implicit none
      include 'mpif.h'
! mode= 1  : Likelihood Maximization for tdim*ndim
! mode= 2  : Likelihood Maximization for      ndim
! mode= 3  : Likelihood Maximization for      1
! mode=10~ : Yhat       Minimization
! mode=20~ : EI         Maximization
      integer, intent(in) :: mode,id_proc,num_proc
      integer, intent(in) :: ndim,tdim,iscf,nevamax,ndebug
      double precision, dimension(ndim+1), intent(in)  :: bound
      double precision, dimension(ndim,tdim), intent(inout) :: theta,power

      integer :: ierr
      integer :: nmax, mmax
      parameter (nmax=50, mmax=50)

!     nmax is the dimension of the largest problem to be solved.
!     mmax is the maximum number of limited memory corrections.
!     Declare the variables needed by the code.
!       A description of all these variables is given at the end of 
!       driver1.
 
      character(len=60) :: task, csave
      logical  ::      lsave(4)
      integer  ::      n, m, iprint,  &
                       nbd(nmax), iwa(3*nmax), isave(44)
      double precision f, factr, pgtol, norm2, &
                       x(nmax), l(nmax), u(nmax), g(nmax), dsave(29), &
                       wa(2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax)

!     Declare a few additional variables for the sample problem.
      integer          i,k,t
 
!     We suppress the default output.

      iprint = -1

!     We suppress both code-supplied stopping tests because the
!        user is providing his own stopping criteria.

      factr=0.0d0
      pgtol=0.0d0

!     We specify the dimension n of the sample problem and the number
!        m of limited memory corrections stored.  (n and m should not
!        exceed the limits nmax and mmax respectively.)
 
      if(mode.ge.1.and.mode.le.3)then
         if(mode.eq.1)then
            n=ndim*tdim
         else if(mode.eq.2)then
            n=ndim
         else if(mode.eq.3)then
            n=1
         end if
         if(iscf.eq.1)then
            n=n*2
         end if
      else if(mode.ge.10.and.mode.le.30)then
         n = ndim
      else
         stop'unknown mode for BFGS'
      end if

      m=nevamax
      if(id_proc.eq.0)then
         write(*,'(1x,a,99f8.3)')'>> BFGS Start from ', &
         ((theta(i,t),i=1,ndim),t=1,tdim)
      end if

!     We now provide nbd which defines the bounds on the variables:
!                    l   specifies the lower bounds,
!                    u   specifies the upper bounds.
!       nbd(i)=0 if x(i) is unbounded,
!              1 if x(i) has only a lower bound,
!              2 if x(i) has both lower and upper bounds, 
!              3 if x(i) has only an upper bound. 

      do 10 i=1,n
         nbd(i)=2
  10  continue
      do i=1,ndim
         l(i) = bound(ndim+1)
         u(i) = bound(i)
      end do
      if(iscf.eq.1.and.mode.le.3)then
       do i=1,ndim
         l(ndim+i) = 1.d0
         u(ndim+i) = 2.d0
       end do
      end if


!     We now define the starting point.

      if(mode.eq.1)then
        do t=1,tdim
          do i=1,ndim
            x(          (t-1)*ndim+i) = theta(i,t)
            if(iscf.eq.1) &
            x(ndim*tdim+(t-1)*ndim+i) = power(i,t)
          end do
        end do
      else if(mode.eq.2)then
        do i=1,ndim
          x(     i) = theta(i,1)
          if(iscf.eq.1) &
          x(ndim+i) = power(i,1)
        end do
      else if(mode.eq.3)then
        x(1) = theta(1,1)
        x(2) = power(1,1)
      else
        do i=1,ndim
          x(i) = theta(i,1)
        end do
      end if

!     We now write the heading of the output.

!      write (6,16)
!  16  format(/,5x, 'Minimize Rosenbrocks function:',
!     +       /,5x, 'f=0 at the optimal solution for x_i=1 for all i!',/)               
               

!     We start the iteration by initializing task.
! 
      task = 'START'

!        ------- the beginning of the loop ----------
 
 111  continue
      
!     This is the call to the L-BFGS-B code.
      if(id_proc.eq.0)then
        call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint, &
                    csave,lsave,isave,dsave)
!       write(*,*)'*task = ',task(1:20)
      end if
      call MPI_BCAST(task,50,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      if (task(1:2) .eq. 'FG') then
!        the minimization routine has returned to request the
!        function f and gradient g values at the current x.

!        Compute function value f for the sample problem.
!        Compute gradient g for the sample problem.

         call eva_bfgs(mode,n,x,f,g)

         if(id_proc.eq.0.and.ndebug.eq.1)then
         write(*,'(11x,a,99f8.3)') &
         '>> Func/Grad Call in BFGS',f,(g(i),i=1,n)
         end if

!          go back to the minimization routine.
         goto 111
      endif
!
      if (task(1:5) .eq. 'NEW_X') then   
!     
!        the minimization routine has returned with a new iterate.
!        At this point have the opportunity of stopping the iteration 
!        or observing the values of certain parameters
!
!        First are two examples of stopping tests.

!        Note: task(1:4) must be assigned the value 'STOP' to terminate  
!          the iteration and ensure that the final results are
!          printed in the default format. The rest of the character
!          string TASK may be used to store other information.

!        1) Terminate if the total number of f and g evaluations
!             exceeds 99.

         call MPI_BCAST(isave(34),1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         if (isave(34) .ge. 999) &
            task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'

!        2) Terminate if  |proj g|/(1+|f|) < 1.0d-10, where 
!           "proj g" denoted the projected gradient

         norm2=0.d0
         do i=1,n
            norm2=norm2+g(i)**2
         end do
         norm2=sqrt(norm2)
         call MPI_BCAST(norm2,1,MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)

         if (norm2 .le. 1.d-10) &
            task='STOP: THE GRADIENT IS SUFFICIENTLY SMALL'

!        We now wish to print the following information at each
!        iteration:
!        
!          1) the current iteration number, isave(30),
!          2) the total number of f and g evaluations, isave(34),
!          3) the value of the objective function f,
!          4) the norm of the projected gradient,  dsve(13)
!
!        See the comments at the end of driver1 for a description
!        of the variables isave and dsave.

         if(id_proc.eq.0) &
         write (*,'(6x,2(a,i5),2(a,e15.5))') &
         '>>',isave(30),' nfg =',isave(34),' f =',f,' |g| =',norm2

!          go back to the minimization routine.
         goto 111

      endif

!     If the run is to be terminated, we print also the information
!     contained in task as well as the final value of x.

      if(id_proc.eq.0)then
        write (*,'(6x,2a)')'>> ',task  
        write (*,'(6x,a,99f8.3)')'>> Final DV by BFGS =',(x(i),i=1,n)
      end if

      call MPI_BCAST(x(1),n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(mode.eq.1)then
        do t=1,tdim
        do i=1,ndim
           theta(i,t) =   x(          (t-1)*ndim+i)
           if(iscf.eq.1)then
             power(i,t) = x(tdim*ndim+(t-1)*ndim+i)
           else
             power(i,t) = 2.d0
           end if
        end do
        end do
      else if(mode.eq.2)then
        do t=1,tdim
        do i=1,ndim
           theta(i,t) =   x(     i)
           if(iscf.eq.1)then
             power(i,t) = x(ndim+i)
           else
             power(i,t) = 2.d0
           end if
        end do
        end do
      else if(mode.eq.3)then
        do t=1,tdim
        do i=1,ndim
           theta(i,t) =   x(1)
           if(iscf.eq.1)then
             power(i,t) = x(2)
           else
             power(i,t) = 2.d0
           end if
        end do
        end do
      else if(mode.ge.11.and.mode.le.30)then
        do i=1,ndim
          power(i,1) = x(i)
        end do
      end if

!           ---------- the end of the loop -------------
 
!     If task is neither FG nor NEW_X we terminate execution.
!     write(*,*)'return from BFGS',id_proc

      return
      end
