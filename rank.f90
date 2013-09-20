        subroutine Rank_New_One
        use dimKrig
        implicit none
        double precision, dimension(ndim)  :: xnew
        double precision, dimension(nfunc) :: ynew,ycon
        integer :: go_to_gradient, limit, irank
        integer :: i,l
        double precision :: pena,penew,fnew,fold,resd
        character(len=6) :: Cmod
        character(len=2) :: Cl
        character(len=29) :: Cdum
        character(len=10) :: Cd

        if(nOPT.ne.1)stop'Rank is only for Mono-objective'
        go_to_gradient =0 ! 1 means to proceed to adjoint analysis

        open(10,file='rank.inp',form='formatted',status='old')
        read(10,*)limit
        close(10)
        if(limit.le.0)limit = 100

        call find_Optimal

        open(10,file='newsample.dat',form='formatted',status='old')
        read(10,*)Cmod
        read(10,*)(xnew(i),i=1,ndim),(ynew(i),i=1,nfunc)
        close(10)
        do l=1,9
          write(Cl,101)l
          if(Cmod(1:2).eq.Cl(1:2))go to 102
        end do
        stop'unknown level in newsample.dat in rank'
101     format(i1,'_')
102     continue

        irank = nsample
        penew = 0.d0
        if(ynew(nfunc).ne.0.d0)go to 999 ! Dummy
        if(l.lt.lopt)then      ! higher fidelity sample
          irank = 1
          go_to_gradient = 1
          go to 999
        else if(l.gt.lopt)then ! lower fidelity sample
          go to 999
        else                   ! same fidelity sample

        end if
        
        call constraint(ynew,penew)

        irank = 1
        do 100 i=1,nsample
           if(info(i)(1:2).ne.Cmod(1:2))go to 100
           if(func(i,nfunc).ne.0.d0)go to 100
           ycon(:) = func(i,:)
           call constraint(ycon,pena)
           if(penew.eq.0.d0.and.pena.ne.0.d0)then
             go to 100
           else if(penew.ne.0.d0.and.pena.eq.0.d0)then
             irank = irank + 1
           else if(penew.ne.0.d0.and.pena.ne.0.d0)then
             fnew = ynew(  nfOPT(1)) + penew
             fold = func(i,nfOPT(1)) + pena
             if(fold.lt.fnew)irank = irank + 1
           else if(penew.eq.0.d0.and.pena.eq.0.d0)then
             fnew = ynew(  nfOPT(1))
             fold = func(i,nfOPT(1))
             if(fold.lt.fnew)irank = irank + 1
           end if
100     continue

        if(penew.eq.0.d0.and.irank.le.limit)then
          go_to_gradient = 1
        else if(penew.ne.0.d0.and.irank.le.limit)then
          go_to_gradient = 0
        end if

999     continue
!       if(ynew(nfOPT(1)).gt.0.02d0)go_to_gradient = 0

        open(10,file='conver.dat',form='formatted',status='unknown')
        read(10,'(a)')Cdum
        read(10,'(a29)')Cdum
        if(Cdum.eq." >> Convergence Problem for F")then
          read(10,*)Cd,Cd,Cd,i,Cd,resd
          write(*,'(6x,a,e15.5)')">> Residual = ",resd
          if(resd.gt.1.d-8)go_to_gradient = 0
        end if
        close(10)

        if(id_proc.eq.0)then
          write(*,'(1x,a,i5)') &
          '>> Rank Check for New Design, Rank = ',irank
          write(*,'(6x,a,2f15.8)') &
          '>> New Sample Info = ',ynew(nfOPT(1)),penew
          write(*,'(6x,a,5x,a)') &
          '>> Fidelity Level  = ',Cmod
          open(10,file='flag_grad.dat',form='formatted',status='unknown')
          write(10,'(i5)')go_to_gradient
          close(10)
        end if

        end subroutine Rank_New_One
