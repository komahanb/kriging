        subroutine Variable_Range
        use dimKrig
        implicit none
        double precision, dimension(ndim)  :: xnew
        double precision, dimension(nfunc) :: ynew,ycon
        integer :: limit, irank
        integer :: i,inew,nsnew,nsold
        double precision :: pena,pnew,fnew,fold
        double precision :: border,vmin,vmax,upfac,downfac
        double precision :: xin,xinmin,vrange,vold

        if(nOPT.ne.1)stop'Rank is only for Mono-objective'

        open(10,file='rank.inp',form='formatted',status='old')
        read(10,*)limit
        close(10)
        if(limit.le.0)limit = 100
        open(10,file='vrange.inp',form='formatted',status='old')
        read(10,*)vrange,nsold
        close(10)
        nsnew = nsample
        if(nsnew.le.nsold  )stop'nsnew<=nsold'
        if(nsnew.gt.nsold+5)stop'nsnew>>nsold+5'

        nsample = nsold
        call find_Optimal  ! Get Optimal without new deisgns
        nsample = nsnew

! Find Best within New Designs
        fnew = 1.d10
        pnew = 1.d0
        inew = 0
        do 100 i=nsold+1,nsnew
           if(func(i,nfunc).ne.0.d0)go to 100
           ycon(:) = func(i,:)
           call constraint(ycon,pena)
           if(pnew.ne.0.d0.and.pena.ne.0.d0)then
             if(func(i,nfOPT(1)).lt.fnew)then
               inew = i
               fnew = func(i,nfOPT(1))
               pnew = pena
             end if
           else if(pnew.eq.0.d0.and.pena.ne.0.d0)then
           else if(pnew.ne.0.d0.and.pena.eq.0.d0)then
               inew = i
               fnew = func(i,nfOPT(1))
               pnew = pena
           else if(pnew.eq.0.d0.and.pena.eq.0.d0)then
             if(func(i,nfOPT(1)).lt.fnew)then
               inew = i
               fnew = func(i,nfOPT(1))
               pnew = pena
             end if
           end if
100     continue
        if(inew.le.nsold.or.inew.gt.nsnew)then
          write(*,*)'*no design in new designs'
          go to 999
        end if
        if(pnew.ne.0.d0)then
          write(*,*)'*unfeasible new design'
          go to 999
        end if
        if(id_proc.eq.0) &
        write(*,'(6x,a,3(i4,a))') &
        '>> Best New Design = ',inew,' E [',nsold+1,' :',nsnew,']'

        xnew(:) = sample(inew,:)
        ynew(:) = func(  inew,:)
        if(ynew(nfunc).ne.0.d0)stop'inew is dummy?'

        irank = 1
        do 200 i=1,nsold
           if(func(i,nfunc).ne.0.d0)go to 200
           ycon(:) = func(i,:)
           call constraint(ycon,pena)
           if(pnew.eq.0.d0.and.pena.ne.0.d0)then
             go to 200
           else if(pnew.ne.0.d0.and.pena.eq.0.d0)then
             irank = irank + 1
           else if(pnew.ne.0.d0.and.pena.ne.0.d0)then
             fnew = ynew(  nfOPT(1)) + pnew
             fold = func(i,nfOPT(1)) + pena
             if(fold.lt.fnew)irank = irank + 1
           else if(pnew.eq.0.d0.and.pena.eq.0.d0)then
             fnew = ynew(  nfOPT(1))
             fold = func(i,nfOPT(1))
             if(fold.lt.fnew)irank = irank + 1
           end if
200     continue

        if(id_proc.eq.0)then
          write(*,'(1x,a,i5)') &
          '>> Rank Check for New Design, Rank = ',irank
          write(*,'(6x,a,2f15.8)') &
          '>> New Sample Info = ',ynew(nfOPT(1)),pnew
        end if

! Trust Region Method for Vrange
        if(IDopt(1).lt.1.or.IDopt(1).gt.nsold)then
           if(id_proc.eq.0)&
           write(*,*)'*No current optimal, cant define vrange'
           go to 999
        end if
        vold   = vrange
        xinmin = 1.d0
        do i=1,ndim
           xin = 0.5d0 + (xnew(i)-sample(IDopt(1),i))/vrange
           if(xin.lt.-0.01d0.or.xin.gt.1.01d0)then
              write(*,'(4f15.8)')xin,xnew(i),sample(IDopt(1),i),vrange
              stop'*out of range, xin'
           end if
           if(xin.lt.0.d0)xin=0.d0
           if(xin.gt.1.d0)xin=1.d0
           if(xin.ge.0.5d0)then
              xinmin = min(xinmin, 1.d0-xin)
           else
              xinmin = min(xinmin, xin)
           end if
        end do
        if(xinmin.lt.0.d0.or.xinmin.gt.0.5d0)then
          write(*,'(f15.8)')xinmin
          stop'range of xinmin'
        end if

        open(10,file='vrangeset.inp',form='formatted',status='old')
        read(10,*)border
        read(10,*)vmin,vmax
        read(10,*)upfac,downfac
        close(10)
        if(border.lt.0.d0.or.border.ge.0.5d0)stop'border for vrangeset'
        if(vmin.ge.vmax)stop'vmin>=vmax for vrangeset'
        if(vmin.le.0.d0)stop'vmin<=0 for vrangeset'
        if(vmax.ge.1.d0)stop'vmax>=1 for vrangeset'
        if(upfac.le.1.d0)stop'upfac<=1 for vrangeset'
        if(downfac.ge.1.d0)stop'downfac>=1 for vrangeset'

        if(xinmin.le.border.and.irank.le.limit)then
           vrange = vrange * upfac
        else if(irank.gt.limit)then
           vrange = vrange * downfac
        end if
        if(vrange.lt.vmin)vrange = vmin
        if(vrange.gt.vmax)vrange = vmax

        if(id_proc.eq.0)then
          if(vrange.ne.vold)then
          write(*,'(6x,2(a,f12.7))') &
          '>> Vrange is modified from',vold,' to',vrange
          else
          write(*,'(6x,1(a,f12.7))') &
          '>> Vrange is preserved to',vrange
          end if
          open(10,file='vrange.inp',form='formatted',status='unknown')
          write(10,*)vrange,nsnew
          close(10)
        end if
        return

!+++ for Error cases
999     continue
        if(id_proc.eq.0)then
          open(10,file='vrange.inp',form='formatted',status='unknown')
          write(10,*)vrange,nsnew
          close(10)
        end if
        end subroutine Variable_Range
