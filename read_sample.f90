        subroutine Read_Sample
        use dimKrig
        implicit none
        integer :: i,j,k,l,ll
        integer :: level,ict
        character(len=6) :: cdum
        character(len=2) :: cl,cll
        double precision :: dummy,prd
        double precision, dimension(ndim,ndim) :: H
        double precision, dimension(ndim,1) :: HV,V

        open(10,file='sample.dat',form='formatted',status='old')
        read(10,*)ndim,zsample,nfunc
        nhes = (ndim*ndim-ndim)/2 + ndim
        if(id_proc.eq.0)then
          write(filenum,*)'>> Reading sample.dat'
          write(filenum,'(6x,a,3i5)') &
          '>> # of DV, Func, Samples = ',ndim,nfunc,zsample
        end if
        if(nfunc-1.ne.nKRI+nDCK+nICK)stop'nfunc != nKRI+nCOK'
        if(nfunc-1.lt.nDMF          )stop'nfunc < nDMF'
        if(nfunc-1.lt.nCON          )stop'nfunc < nCON'
        if(nfunc-1.lt.nOPT          )stop'nfunc < nOPT'

        call Allocating

        ict_sample = 0
        ict_dummy  = 0
        do 100 i=1,zsample
           read(10,*)info(i)
           backspace(10)

           ! for previous format
           if(     info(i).eq.'0 ')then
             info(i) = '1_F   '
           else if(info(i).eq.'1 ')then
             info(i) = '1_FG  '
           else if(info(i).eq.'2 ')then
             info(i) = '1_FH  '
           else if(info(i).eq.'3 ')then
             info(i) = '1_FGH '
           else if(info(i).eq.'4 ')then
             info(i) = '1_FHv '
           else if(info(i).eq.'5 ')then
             info(i) = '1_FGHv'
           end if

           ! Fidelity Level
           do l=1,9
             write(cl,101)l
101          format(i1,'_')
             if(info(i)(1:2).eq.cl)then
               level = l
               go to 102
             end if
           end do
           write(*,*)'*Unknown level',i,info(i)
           call stop_all
102        continue

           ! Derivative Information
           if(info(i)(3:6).eq.'F   ')then              ! only func
              read(10,*,err=1000,end=1000)   &
              cdum,(sample(i,j           ),j=1,ndim),           &
              (       func(i,j           ),j=1,nfunc)
              ict_sample(level,1) = ict_sample(level,1) + 1
           else if(info(i)(3:6).eq.'FG  ')then         ! func and grad
              read(10,*,err=1000,end=1000)   &
              cdum,(sample(i,j           ),j=1,ndim),           &
              (       func(i,j           ),j=1,nfunc),          &
              ((     gfunc(i,nfCOK(j),k  ),k=1,ndim),j=1,nCOK)
              ict_sample(level,2) = ict_sample(level,2) + 1
           else if(info(i)(3:6).eq.'FH  ')then         ! func and hess
              read(10,*,err=1000,end=1000)   &
              cdum,(sample(i,j           ),j=1,ndim),           &
              (       func(i,j           ),j=1,nfunc),          &
              ((     dummy                ,k=1,ndim),j=1,nCOK), &
              (((    hfunc(i,nfCOK(j),k,l),l=1,ndim),k=1,ndim),j=1,nCOK)
              ict_sample(level,3) = ict_sample(level,3) + 1
           else if(info(i)(3:6).eq.'FGH ')then         ! func,grad and hess
              read(10,*,err=1000,end=1000)                      &
              cdum,(sample(i,j           ),j=1,ndim),           &
              (       func(i,j           ),j=1,nfunc),          &
              ((     gfunc(i,nfCOK(j),k  ),k=1,ndim),j=1,nCOK), &
              (((    hfunc(i,nfCOK(j),k,l),l=1,ndim),k=1,ndim),j=1,nCOK)
              ict_sample(level,4) = ict_sample(level,4) + 1
           else if(info(i)(3:6).eq.'FHv ')then         ! func and Hv
              read(10,*,err=1000,end=1000)                      &
              cdum,(sample(i,j           ),j=1,ndim),           &
              (       func(i,j           ),j=1,nfunc),          &
              ((     dummy                ,k=1,ndim),j=1,nCOK), &
              ((     hvect(i,nfCOK(j),k,1),k=1,ndim),j=1,nCOK), &
              ((     hvect(i,nfCOK(j),k,2),k=1,ndim),j=1,nCOK)
              ict_sample(level,5) = ict_sample(level,5) + 1
           else if(info(i)(3:6).eq.'FGHv')then         ! func,grad and Hv
              read(10,*,err=1000,end=1000)                      &
              cdum,(sample(i,j           ),j=1,ndim),           &
              (       func(i,j           ),j=1,nfunc),          &
              ((     gfunc(i,nfCOK(j),k  ),k=1,ndim),j=1,nCOK), &
              ((     hvect(i,nfCOK(j),k,1),k=1,ndim),j=1,nCOK), &
              ((     hvect(i,nfCOK(j),k,2),k=1,ndim),j=1,nCOK)
              ict_sample(level,6) = ict_sample(level,6) + 1
           else
              write(*,*)'*Unknown information',i,info(i)
              call stop_all
           end if

           if(func(i,nfunc).ne.0.d0)then
             ict_dummy(level) = ict_dummy(level) + 1
           end if
100     continue ! zsample loop
        close(10)

        ! Make Sum
        lmax = 0
        do l=1,9
          do i=1,6
             ict_sample(l,0) = ict_sample(l,0) + ict_sample(l,i)
          end do
          if(Cmode(:6).eq.'Update'.or.Cmode(:4).eq.'Rank')then
            if(ict_sample(l,0).ne.0)lmax = l
          else
            if(ict_sample(l,0).ne.0)lmax = lmax + 1
          end if
        end do

        ! H to Hv for special case
        if(nMXV.ne.0.and.nMXH.eq.0)then
          write(filenum,'(a)')'*all G/H samples are changed to G/Hv...'
          do i=1,zsample
            if(info(i)(3:6).eq.'FGH ')then
              do j=1,nCOK
                H(:,:) = hfunc(i,nfCOK(j),:,:)
                prd = 0.d0
                do k=1,ndim
                   prd = prd + gfunc(i,nfCOK(j),k)**2
                end do
                do k=1,ndim
                  V(k,1) = gfunc(i,nfCOK(j),k) / dsqrt(prd)
                  V(k,1) = 1.d0 / dsqrt(dble(ndim))
                end do
                HV = matmul(H,V)
                do k=1,ndim
                  hvect(i,nfCOK(j),k,1) =  V(k,1)
                  hvect(i,nfCOK(j),k,2) = HV(k,1)
                end do
              end do
              info(i)(3:6) = 'FGHv'
            end if
          end do
          do l=1,9
            ict_sample(l,6) = ict_sample(l,6) + ict_sample(l,4)
            ict_sample(l,4) = 0
          end do
        end if

        ! Cram up Fidelity Levels
        if(Cmode(:6).eq.'Update')go to 230
        if(Cmode(:4).eq.'Rank')go to 230
200     continue
        do l=1,9
          if(l.le.lmax)then
            if(ict_sample(l,0).eq.0)go to 210
          else
            if(ict_sample(l,0).ne.0)go to 210
          end if
        end do
        go to 230
210     continue
        if(id_proc.eq.0) &
        write(filenum,'(a)')'*Cram up Fidelity Levels'
        do l=1,9
          if(ict_sample(l,0).eq.0)then
            do ll=l+1,9
              ict_sample(ll-1,:) = ict_sample(ll,:)
              ict_dummy(ll-1)    = ict_dummy(ll)
            end do
            do 220 i=1,zsample
              do ll=l+1,9
                write(cl,101)ll
                if(info(i)(1:2).eq.cl)then
                  write(cll,101)ll-1
                  info(i)(1:2) = cll(1:2)
                end if
              end do
220         continue
            go to 200
          end if
        end do
230     continue

        ict_sample(0,:) = 0
        do i=1,6
         do l=1,lmax
          ict_sample(0,i) = ict_sample(0,i) + ict_sample(l,i)
         end do
        end do
        ict_dummy(0) = 0
        do l=1,lmax
          ict_dummy(0) = ict_dummy(0) + ict_dummy(l)
        end do

        if(id_proc.eq.0)then
          write(filenum,'(6x,a,i5)')'>> # of Fidelity Levels   = ',lmax
          do l=1,lmax
            write(filenum,'(6x,a,i2)')'>> Fidelity Level : ',l
            if(ict_sample(l,1).ne.0) &
            write(filenum,'(11x,a,i5)')'>> # of Samples(F     )    = ',ict_sample(l,1)
            if(ict_sample(l,2).ne.0) &
            write(filenum,'(11x,a,i5)')'>> # of Samples(F,G   )    = ',ict_sample(l,2)
            if(ict_sample(l,3).ne.0) &
            write(filenum,'(11x,a,i5)')'>> # of Samples(F,  H )    = ',ict_sample(l,3)
            if(ict_sample(l,4).ne.0) &
            write(filenum,'(11x,a,i5)')'>> # of Samples(F,G,H )    = ',ict_sample(l,4)
            if(ict_sample(l,5).ne.0) &
            write(filenum,'(11x,a,i5)')'>> # of Samples(F,  Hv)    = ',ict_sample(l,5)
            if(ict_sample(l,6).ne.0) &
            write(filenum,'(11x,a,i5)')'>> # of Samples(F,G,Hv)    = ',ict_sample(l,6)
            if(ict_dummy(l).ne.0) &
            write(filenum,'(11x,a,i5)')'>> # of Samples with Dummy = ',ict_dummy(l)
          end do
        end if

        ict = 0
        do l=1,lmax
          if(ict_sample(l,3).ne.0)then ! temporary
            write(*,*)'*Only Hessian is not effective yet'
            call stop_all
          end if
          if(ict_sample(l,5).ne.0)then ! temporary
            write(*,*)'*Only Hessian Vector is not effective yet'
            call stop_all
          end if
          ict = ict + ict_sample(l,0)
        end do
        if(ict.ne.zsample)stop'ict_sample != zsample'

        nsample = zsample 
        tdim    = lmax*(lmax+1)/2
        return

1000    continue
        write(*,*)'*Stop by Unexpected Data Arrangement in sample.dat',i
        call stop_all
        end subroutine Read_Sample
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine Allocating
        use dimKrig
        implicit none
        integer :: iflag,istat

        iflag = 0
        allocate( info(zsample),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate( inf( zsample),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate( icone(zsample),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1

        allocate( sample(zsample,ndim),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate( sampl( zsample,ndim),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1

        allocate(  func(zsample,nfunc),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate(  fun( zsample      ),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1

        allocate( gfunc(zsample,nfunc,ndim),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate( gfun( zsample      ,ndim),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1

        allocate( hfunc(zsample,nfunc,ndim,ndim),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate( hfun( zsample      ,ndim,ndim),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate( hvect(zsample,nfunc,ndim,2),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate( hvec( zsample      ,ndim,2),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1

        allocate( mreg(kreg_orig,ndim),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
       
        allocate( nparent(zsample),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate( dxadd(zsample),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1

        if(nICK.ne.0)then
        allocate( tdxx(zsample),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        allocate( tdx( zsample),stat = istat )
                  if(istat.ne.0)iflag = iflag + 1
        end if
        if(iflag.ne.0)stop'*Allocation Miss in read_samples'

        info    = '      '
        inf     = '      '
        icone   = 0
        sample  = 0.d0
        sampl   = 0.d0
        func    = 0.d0
        fun     = 0.d0
        gfunc   = 0.d0
        gfun    = 0.d0
        hfunc   = 0.d0
        hfun    = 0.d0
        hvect   = 0.d0
        hvec    = 0.d0
        mreg    = 0
        nparent = 0
        dxadd   = 0.d0
        if(nICK.ne.0)then
          tdxx  = 0.d0
          tdx   = 0.d0
        end if

        end subroutine Allocating
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine Deallocating
        use dimKrig
        implicit none
        integer :: ideb
        ideb = 0
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'1'
        deallocate(info,inf)
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'2'
        deallocate(icone)
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'2-2'
        deallocate(nparent,dxadd)
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'3'
        deallocate(sample,sampl)
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'4'
        deallocate( func, fun)
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'5'
        deallocate(gfunc,gfun)
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'6'
        deallocate(hfunc,hfun,hvect,hvec)
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'7'
        if(nICK.ne.0)then
        deallocate(tdxx,tdx)
        end if
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'8'
        deallocate(mreg)
        if(id_proc.eq.0.and.ideb.eq.1)write(*,*)'9'
        end subroutine Deallocating
