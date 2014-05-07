        subroutine reduce_data(ifunc,igrad,dlim,nlim)
        use dimKrig
        implicit none
! to reduce dummy, by using nMXS, dlim, nlim
! to limit Grad/Hess info by nMXG, nMXH, nMXV
! make sampl from sample
        integer, intent(in)  :: ifunc,nlim
        double precision, intent(in) :: dlim
        integer, intent(out) :: igrad
        integer :: i

        if(id_proc.eq.0.and.Cmode(:4).eq.'Make')then
        write(filenum,'(1x,a,i3)')'>> Reduce Sample Points for Function',ifunc
        end if
        igrad = 0
        do i=1,nDCK
          if(ifunc.eq.nfDCK(i))igrad = 1
        end do
        do i=1,nICK
          if(ifunc.eq.nfICK(i))igrad = 2
        end do

        if(igrad.eq.2)then
          call read_trustdx(ifunc)
        end if

        do i=1,nDMF
           if(ifunc.eq.nfDMF(i))then
              rsample = nsample - ict_dummy(0)
              go to 100
           end if
        end do
        rsample = nsample
100     continue

        call REDUCE_DUMMY(ifunc,igrad)
        call REDUCE_SAMPLE(igrad,dlim,nlim)

        if( (nMXG.ge.0.and.nMXG.lt.rsample).or. &
            (nMXH.ge.0.and.nMXH.lt.rsample).or. &
            (nMXV.ge.0.and.nMXV.lt.rsample) )then
          call REDUCE_GRADHESS(igrad)
        end if

        call FIDELITY_CHECK(ifunc,igrad)
        call report_data(ifunc,igrad)

        end subroutine reduce_data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine FIDELITY_CHECK(ifunc,igrad)
        use dimKrig
        implicit none
        integer, intent(in) :: ifunc,igrad
        integer :: i,j,k,ict
        double precision :: eps,dd
        integer,          dimension(rsample) :: iflg
        character(len=6), dimension(rsample)           :: inft
        integer, dimension(rsample)                    :: iconet
        double precision, dimension(rsample,ndim)      :: samplt
        double precision, dimension(rsample)           :: funt,tdxt
        double precision, dimension(rsample,ndim)      :: gfunt
        double precision, dimension(rsample,ndim,ndim) :: hfunt
        double precision, dimension(rsample,ndim,2)    :: hvet

        do i=1,nVFM
          if(nfVFM(i).eq.ifunc)return
        end do

        ! not variable fidelity function
        lmax = 1
        tdim = 1

        iflg = 0
        eps  = 1.d-5
        do 100 i=1,rsample
          inf(i)(1:2) = '1_'
          do 110 j=i+1,rsample
            if(i.eq.j)go to 110
            dd = 0.d0
            do k=1,ndim
               dd = dd + (sampl(i,k)-sampl(j,k))**2
            end do
            dd = sqrt(dd)
            if(dd.le.eps)then
              iflg(i) = 1
            end if
110       continue
100     continue

        ict = 0
        do 200 i=1,rsample
          if(iflg(i).eq.1)go to 200
          ict = ict + 1
          inft(ict)        = inf(i)
          iconet(ict)      = icone(i)
          samplt(ict,:)    = sampl(i,:)
          funt(ict)        = fun(i)
          gfunt(ict,:)     = gfun(i,:)
          hfunt(ict,:,:)   = hfun(i,:,:)
          hvet(ict,:,:)    = hvec(i,:,:)
          if(igrad.eq.2) &
          tdxt(ict)        = tdx(i)
200     continue
        inf   = '      '
        icone = 0
        sampl = 0.d0
        fun   = 0.d0
        gfun  = 0.d0
        hfun  = 0.d0
        hvec  = 0.d0
        do 210 i=1,ict
          inf(i)      = inft(i)
          icone(i)    = iconet(i)
          sampl(i,:)  = samplt(i,:)
          fun(i)      = funt(i)
          gfun(i,:)   = gfunt(i,:)
          hfun(i,:,:) = hfunt(i,:,:)
          hvec(i,:,:) = hvet(i,:,:)
210     continue
        if(ict.ne.rsample)then
          if(id_proc.eq.0) &
          write(filenum,'(6x,a,i4,a,i4,a)') &
          '>> Reduce Samples   from ',rsample,' to',ict, &
          ' for Non-VFM Case'
        end if

        end subroutine FIDELITY_CHECK
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine REDUCE_GRADHESS(igrad)
        use dimKrig
        implicit none
        integer, intent(in) :: igrad
        integer :: i,j,k,itarg,ict,trank,ngh,ng,nh,nv,iter,imin
        double precision :: dd,ddmin
        integer,          dimension(rsample) :: rank,isub
        double precision, dimension(rsample) :: dist

        ng = 0
        nh = 0
        nv = 0
        do i=1,rsample
           if(inf(i)(3:6).eq.'FG  ')then
             ng = ng + 1
           else if(inf(i)(3:6).eq.'FGH ')then
             nh = nh + 1
           else if(inf(i)(3:6).eq.'FGHv')then
             nv = nv + 1
           end if
        end do
        if(ng.le.nMXG.and.nh.le.nMXH.and.nv.le.nMXV)return
        if(nOPT.ne.1)stop'REDUCE_SAMPLE is for mono-objective'
        if(IDopt(1).le.0.or.IDopt(1).gt.nsample)return
        ngh = ng + nh + nv

! optimal design
        itarg = 0
        do i=1,rsample
           if(icone(i).eq.IDopt(1))then
             itarg = i
           end if
        end do
        if(itarg.le.0.or.itarg.gt.rsample)stop'unknown itarg'

! rank by distance
        ict = 0
        do 100 i=1,rsample
           if(inf(i)(3:6).eq.'F   ')go to 100
           dd = 0.d0
           do 110 k=1,ndim
              dd = dd + ( sampl(i,k)-sampl(itarg,k) )**2
110        continue
!          dd = (fun(nfOPT(1))-OFopt(1))**2
           ict = ict + 1
           dist(ict) = dd
           rank(ict) = i
100     continue
        do 120 i=1,ngh-1
        do 120 j=i+1,ngh
          if(dist(i).gt.dist(j))then
            trank   = rank(j)
            rank(j) = rank(i)
            rank(i) = trank
            dd      = dist(j)
            dist(j) = dist(i)
            dist(i) = dd
          end if
120     continue
        do 130 i=1,ngh-1
           if(dist(i).gt.dist(i+1))stop'sorting in reduce_gradhess'
130     continue

! eliminating
        if(id_proc.eq.0.and.Cmode(:4).eq.'Make')then
          if(nMXH.ge.0) &
          write(filenum,'(6x,a,i4,a,i4,a)') &
          '>> Reduce Grad/Hess from ',nh,' to',nMXH,' by distance'
          if(nMXV.ge.0) &
          write(filenum,'(6x,a,i4,a,i4,a)') &
          '>> Reduce Grad/Hvec from ',nv,' to',nMXV,' by distance'
          if(nMXG.ge.0) &
          write(filenum,'(6x,a,i4,a,i4,a)') &
          '>> Reduce Grad      from ',ng+(nh-nMXH)+(nv-NMXV), &
                                         ' to',nMXG,' by distance'
        end if 
        ! Limit H
        if(nMXH.ne.0)then
          ict = 0
          do 200 i=1,ngh
            if(inf(rank(i))(3:6).eq.'FGH ')then
              if(ict.lt.nMXH)then
               ict = ict + 1
              else
               inf(rank(i))(3:6) = 'FG  '
              end if
            end if
200       continue
        else
          do 210 i=1,rsample
            if(inf(i)(3:6).eq.'FGH ')inf(i)(3:6)='FG  '
210       continue
        end if
        ! Limit Hv
        if(nMXV.ne.0)then
          ict = 0
          do 220 i=1,ngh
            if(inf(rank(i))(3:6).eq.'FGHv')then
              if(ict.lt.nMXV)then
               ict = ict + 1
              else
               inf(rank(i))(3:6) = 'FG  '
              end if
            end if
220       continue
        else
          do 230 i=1,rsample
            if(inf(i)(3:6).eq.'FGHv')inf(i)(3:6)='FG  '
230       continue
        end if
        ! Limit G
        if(nMXG.ne.0)then
          ict = 0
          do 300 i=1,ngh
            if(inf(rank(i))(3:6).eq.'FG  ')then
              if(ict.lt.nMXG)then
               ict = ict + 1
              else
               inf(rank(i))(3:6) = 'F   '
              end if
            end if
300       continue
        else
          do 310 i=1,rsample
            if(inf(i)(3:6).eq.'FG ')inf(i)(3:6)='F   '
310       continue
        end if
!       gfun,hfun,hvec,tdx are preserved

        return
        end subroutine REDUCE_GRADHESS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine REDUCE_SAMPLE(igrad,dlim,nlim)
        use dimKrig
        implicit none
        integer, intent(in) :: igrad,nlim
        double precision, intent(in) :: dlim
        integer :: i,j,k,itarg,ict,trank
        integer :: nct,gct,hct,vct,ng,nh,nv
        character(len=4) :: ni
        double precision :: dd,ddlim
        integer,          dimension(rsample) :: rank
        double precision, dimension(rsample) :: dist
        character(len=6), dimension(rsample)           :: inft
        integer, dimension(rsample)                    :: iconet
        double precision, dimension(rsample,ndim)      :: samplt
        double precision, dimension(rsample)           :: funt,tdxt
        double precision, dimension(rsample,ndim)      :: gfunt
        double precision, dimension(rsample,ndim,ndim) :: hfunt
        double precision, dimension(rsample,ndim,2)    :: hvet

        if(nMXS.le.0.and.dlim.le.0.d0.and.nlim.le.0)return
! nMXS : use till nMXS-th nearest samples from current optimal
! dlim : use samples within the range of dlim
! nlim : specify the size of correlation matrix
        if(nOPT.ne.1)stop'REDUCE_SAMPLE is for mono-objective'
        if(IDopt(1).le.0.or.IDopt(1).gt.nsample)return


! optimal design
        itarg = 0
        do i=1,rsample
           if(icone(i).eq.IDopt(1))then
             itarg = i
           end if
        end do
        if(itarg.le.0.or.itarg.gt.rsample)stop'unknown itarg red_sam'

! rank by distance
        do 100 i=1,rsample
           dd = 0.d0
           do 110 k=1,ndim
              dd = dd + ( sampl(i,k)-sampl(itarg,k) )**2
110        continue
!          dd = (fun(nfOPT(1))-OFopt(1))**2
           dist(i) = dd
           rank(i) = i
100     continue
        do 120 i=1,rsample-1
        do 120 j=i+1,rsample
          if(dist(i).gt.dist(j))then
            trank   = rank(j)
            rank(j) = rank(i)
            rank(i) = trank
            dd      = dist(j)
            dist(j) = dist(i)
            dist(i) = dd
          end if
120     continue
        do 130 i=1,rsample-1
           if(dist(i).gt.dist(i+1))stop'sorting in reduce_sample'
!          write(*,'(2i5,2f15.8)')i,rank(i),dist(i),sampl(rank(i),1)
130     continue

        if(nMXS.ne.0.and.nMXS.le.rsample)then
          ddlim = dist(nMXS)
        else if(dlim.ne.0.d0)then
          ddlim = dlim
        else if(nlim.ne.0)then
          go to 250
        else
          stop'ddlim in reduce_sample'
        end if

! eliminating by distance
        ict   = 0
        do 200 i=1,rsample
          do 210 j=1,rsample
            if(rank(j).eq.i)go to 220
210       continue
          stop'no target sample in reduce_sample'
220       continue
          if(dist(j).gt.ddlim)go to 200

          ict = ict + 1
          inft(ict)        = inf(i)
          iconet(ict)      = icone(i)
          samplt(ict,:)    = sampl(i,:)
          funt(ict)        = fun(i)
          if(inf(i)(3:6).eq.'FG  '.or.inf(i)(3:6).eq.'FGH ' &
                                  .or.inf(i)(3:6).eq.'FGHv')then
            gfunt(ict,:)   = gfun(i,:)
          else
            gfunt(ict,:)   = 0.d0
          end if
          if(inf(i)(3:6).eq.'FH  '.or.inf(i)(3:6).eq.'FGH ')then
            hfunt(ict,:,:) = hfun(i,:,:)
          else
            hfunt(ict,:,:) = 0.d0
          end if
          if(inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv')then
            hvet(ict,:,:) = hvec(i,:,:)
          else
            hvet(ict,:,:) = 0.d0
          end if
          if(igrad.eq.2)then
            tdxt(ict)      = tdx(i)
          end if
200     continue
        go to 299

250     continue
! eliminating by nlim
        ddlim = 0.d0
        ict = 0
        nct = 0
        gct = 0
        hct = 0
        vct = 0
        if(igrad.eq.1)then
          ng = ndim
          nv = ndim
          if(mode_dck.eq.1)then
             nh = ndim + ndim
          else
             nh = nhes + ndim
          end if
        else if(igrad.eq.2)then
          if(ngput.eq.-2)then
            ng = ndim
          else if(ngput.eq.-1)then
            ng = 1
          end if
          if(nhput.eq.-4)then
            nh = nhes + ndim
          else if(nhput.eq.-3)then
            nh = ndim*2
          else if(nhput.eq.-2)then
            nh = ndim
          else if(nhput.eq.-1)then
            nh = 1
          else
            nh = nhput
          end if
          nv = 0  ! temporary
        end if
        do 260 j=1,rsample
           i   = rank(j)
           if(nct.ge.nlim)go to 270
           if(     inf(i)(3:6).eq.'FGHv'.and.nct.le.nlim-(1+ng+nv))then
              ni = 'FGHv'
           else if(inf(i)(3:6).eq.'FGHv'.and.nct.le.nlim-(1+ng))then
              ni = 'FG  '
           else if(inf(i)(3:6).eq.'FGHv'.and.nct.le.nlim-(1))then
              ni = 'F   '
           else if(inf(i)(3:6).eq.'FHv '.and.nct.le.nlim-(1+nv))then
              ni = 'FHv '
           else if(inf(i)(3:6).eq.'FHv '.and.nct.le.nlim-(1))then
              ni = 'F   '
           else if(inf(i)(3:6).eq.'FGH '.and.nct.le.nlim-(1+nh))then
              ni = 'FGH '
           else if(inf(i)(3:6).eq.'FGH '.and.nct.le.nlim-(1+ng))then
              ni = 'FG  '
           else if(inf(i)(3:6).eq.'FGH '.and.nct.le.nlim-(1))then
              ni = 'F   '
           else if(inf(i)(3:6).eq.'FG  '.and.nct.le.nlim-(1+ng))then
              ni = 'FG  '
           else if(inf(i)(3:6).eq.'FG  '.and.nct.le.nlim-(1))then
              ni = 'F   '
           else if(inf(i)(3:6).eq.'F   ')then
              ni = 'F   '
           else
              write(*,'(2a,i7)')'*unknown set for reduce ',ni,nct,nlim
              call stop_all
           end if
           
           if(ni.eq.'FGH ')then
              hct = hct + 1
              if(hct.gt.nMXH)ni = 'FG  '
           end if
           if(ni.eq.'FGHv')then
              vct = vct + 1
              if(vct.gt.nMXV)ni = 'FG  '
           end if
           if(ni.eq.'FG ')then
              gct = gct + 1
              if(gct.gt.nMXG)ni = 'F   '
           end if

           if(     ni.eq.'F   ')then
             nct = nct + 1
           else if(ni.eq.'FG  ')then
             nct = nct + 1 + ng
           else if(ni.eq.'FGH ')then
             nct = nct + 1 + nh
           else if(ni.eq.'FHv ')then
             nct = nct + 1 + nv
           else if(ni.eq.'FGHv')then
             nct = nct + 1 + nv + ng
           end if

           ict = ict + 1
           inft(ict)(1:2)   = inf(i)(1:2)
           inft(ict)(3:6)   = ni(1:4)
           iconet(ict)      = icone(i)
           samplt(ict,:)    = sampl(i,:)
           funt(ict)        = fun(i)
           if(ni.eq.'FG  '.or.ni.eq.'FGH '.or.ni.eq.'FGHv')then
             gfunt(ict,:)   = gfun(i,:)
           else
             gfunt(ict,:)   = 0.d0
           end if
           if(ni.eq.'FH  '.or.ni.eq.'FGH ')then
             hfunt(ict,:,:) = hfun(i,:,:)
           else
             hfunt(ict,:,:) = 0.d0
           end if
           if(ni.eq.'FHv '.or.ni.eq.'FGHv')then
             hvet(ict,:,:) = hvec(i,:,:)
           else
             hvet(ict,:,:) = 0.d0
           end if
           if(igrad.eq.2)then
             tdxt(ict)      = tdx(i)
           end if
           ddlim = max(ddlim,dist(j))
260     continue
270     continue
        if(ddlim.le.0.d0)stop'ddlim in reduce_nlim'
        go to 299

299     continue
        if(dlim.ne.0.d0.or.nlim.ne.0)then
          if(id_proc.eq.0)then
            open(10,file='vrange.inp',form='formatted',status='unknown')
            if(ddlim.ge.0.01d0)then
             write(10,'(f15.8)')ddlim*2.d0
            else
             write(10,'(f15.8)')0.02d0
            end if
            close(10)
          end if
        end if

        inf   = '      '
        icone = 0
        sampl = 0.d0
         fun  = 0.d0
        gfun  = 0.d0
        hfun  = 0.d0
        hvec  = 0.d0
        if(igrad.eq.2)then
        tdx   = 0.d0
        end if

        do 300 i=1,ict
           inf(i)      = inft(i)
           icone(i)    = iconet(i)
           sampl(i,:)  = samplt(i,:)
            fun(i)     =  funt(i)
           gfun(i,:)   = gfunt(i,:)
           hfun(i,:,:) = hfunt(i,:,:)
           hvec(i,:,:) = hvet(i,:,:)
           if(igrad.eq.2)then
           tdx(i)      = tdxt(i)
           end if
300     continue

        if(id_proc.eq.0.and.Cmode(:4).eq.'Make')then
          if(nMXS.ne.0.and.nMXS.le.rsample)then
            write(filenum,'(6x,a,i4,a,i4)') &
            '>> Reduce Samples   from ',rsample,' to',ict
            if(ict.ne.nMXS)stop'ict.ne.nMXS'
          else if(dlim.ne.0.d0)then
            write(filenum,'(6x,a,i4,a,i4,a,f8.4)') &
            '>> Reduce Samples   from ',rsample,' to',ict,' by d=',dlim
          else if(nlim.ne.0)then
            write(filenum,'(6x,a,i4,a,i4,a,i5)') &
            '>> Reduce Samples   from ',rsample,' to',ict,' by n=',nlim
          end if
        end if 
        rsample = ict
        end subroutine REDUCE_SAMPLE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine REDUCE_DUMMY(ifunc,igrad)
        use dimKrig
        implicit none
        integer, intent(in) :: ifunc,igrad
        integer :: i,ict,iCOK

        iCOK = 0
        do i=1,nCOK
          if(nfCOK(i).eq.ifunc)iCOK = i
        end do

        if(nsample.eq.rsample)then ! no dummy
           inf         = info
           sampl       = sample
           fun(:)      = func(:,ifunc)
           if(iCOK.ne.0)then
             gfun(:,:)   = gfunc(:,iCOK,:)
             hfun(:,:,:) = hfunc(:,iCOK,:,:)
             hvec(:,:,:) = hvect(:,iCOK,:,:)
           else
             do i=1,nsample
                inf(i)(3:6) = 'F   '
             end do
           end if
           if(igrad.eq.2)then
             tdx   = tdxx
           end if
           do i=1,nsample
             icone(i) = i
           end do
           return
        else if(nsample.lt.rsample)then
           stop'rsample>nsample'
        end if

        if(id_proc.eq.0.and.Cmode(:4).eq.'Make')then
          write(filenum,'(6x,a,i4,a)') &
          '>> Reduce',ict_dummy(0),' Dummy Points'
        end if 

        ict = 0
        do 100 i=1,nsample
           if(func(i,nfunc).ne.0.d0)go to 100
           ict = ict + 1

           inf(ict)         = info(i)
           icone(ict)       = i
           sampl(ict,:)     = sample(i,:)
           fun(ict)         = func(i,ifunc)

           if(igrad.eq.0)then
             inf(ict)(3:6) = 'F   '
             go to 100
           end if
           if((info(i)(3:6).eq.'FG  '.or.info(i)(3:6).eq.'FGH '.or. &
               info(i)(3:6).eq.'FGHv') .and. (iCOK.ne.0))then
             gfun(ict,:)    = gfunc(i,iCOK,:)
           else
             gfun(ict,:)    = 0.d0
           end if
           if((info(i)(3:6).eq.'FH  '.or.info(i)(3:6).eq.'FGH ') &
                                       .and.(iCOK.ne.0))then
             hfun(ict,:,:)  = hfunc(i,iCOK,:,:)
           else
             hfun(ict,:,:)  = 0.d0
           end if
           if((info(i)(3:6).eq.'FHv '.or.info(i)(3:6).eq.'FGHv') &
                                       .and.(iCOK.ne.0))then
             hvec(ict,:,:)  = hvect(i,iCOK,:,:)
           else
             hvec(ict,:,:)  = 0.d0
           end if
           if(igrad.eq.2)then
             tdx(ict)       = tdxx(i)
           end if
100     continue

        if(ict.ne.rsample)then
          write(*,*)'*ict = ',ict,' rsample = ',rsample
          call stop_all
        end if
        end subroutine REDUCE_DUMMY
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine report_data(ifunc,igrad)
        use dimKrig
        implicit none
        integer, intent(in) :: ifunc,igrad
        integer :: i,j,k,l,level,itarg
        integer, dimension(0:10,0:10) :: ict
        character(len=2) :: cl

        itarg = 0
        do i=1,nDMF
           if(ifunc.eq.nfDMF(i))itarg = i
        end do

        ict = 0
        do 100 i=1,rsample
          if(itarg.ne.0.and.func(icone(i),nfunc).ne.0.d0) &
          stop'still have dummy data'
          if(icone(i).le.0.or.icone(i).gt.nsample)stop'icone_report'

          do l=1,lmax
            write(cl,101)l
101         format(i1,'_')
            if(inf(i)(1:2).eq.cl)then
              level = l
              go to 102
            end if
          end do
          stop'unknown level in report_data'
102       continue

          if(     inf(i)(3:6).eq.'F   ')then
            ict(l,1) = ict(l,1) + 1
          else if(inf(i)(3:6).eq.'FG  ')then
            ict(l,2) = ict(l,2) + 1
          else if(inf(i)(3:6).eq.'FH  ')then
            ict(l,3) = ict(l,3) + 1
          else if(inf(i)(3:6).eq.'FGH ')then
            ict(l,4) = ict(l,4) + 1
          else if(inf(i)(3:6).eq.'FHv ')then
            ict(l,5) = ict(l,5) + 1
          else if(inf(i)(3:6).eq.'FGHv')then
            ict(l,6) = ict(l,6) + 1
          else
            stop'mode in report_data'
          end if

          if((inf(i)(3:6).eq.'FG  '.or.inf(i)(3:6).eq.'FGH ' &
                                   .or.inf(i)(3:6).eq.'FGHv') &
              .and.(igrad.ne.0))then
            do j=1,ndim
              if(gfun(i,j).ne.0.d0)go to 110
            end do
            write(filenum,150)'*All gradient comps are zero in sample # = ',i,icone(i)
          end if
110       continue
          if((inf(i).eq.'FH '.or.inf(i).eq.'FGH ').and.(igrad.ne.0))then
            do j=1,ndim
            do k=1,ndim
              if(hfun(i,j,k).ne.hfun(i,k,j))stop'hessian assymmetry'
            end do
            end do
            do j=1,ndim
            do k=1,ndim
              if(hfun(i,j,k).ne.0.d0)go to 120
            end do
            end do
            write(filenum,150)'*All Hessian  comps are zero in sample # = ',i,icone(i)
          end if
120       continue
          if((inf(i).eq.'FHv '.or.inf(i).eq.'FGHv').and.(igrad.ne.0))then
            do j=1,ndim
              if(hvec(i,j,2).ne.0.d0)go to 130
            end do
            write(filenum,150)'*All H Vector comps are zero in sample # = ',i,icone(i)
          end if
130       continue
100     continue
150     format(a,2i5)

        ! Summing up
        do l=1,lmax
         do i=1,6
          ict(0,i) = ict(0,i) + ict(l,i)
          ict(l,0) = ict(l,0) + ict(l,i)
         end do
        end do

        if(id_proc.eq.0)then
          if(ndebug.eq.1) &
          write(filenum,'(11x,2(a,i4))')'>> # of Samples,',nsample,' ->',rsample
          do l=1,lmax
          write(filenum,'(11x, a,i1,a,9i4)') &
          '>> # of Modes in Fid-',l,' =',(ict(l,i),i=1,6),ict(l,0)
          end do
          write(filenum,'(11x, a,9i4)') &
          '>> # of Modes (sum)    =',(ict(0,i),i=1,6)
        end if

        itarg = 0
        do i=1,rsample
           if(icone(i).eq.IDopt(1))then
             itarg = i
           end if
        end do
        if(itarg.lt.0.or.itarg.gt.rsample)stop'optimal is eliminated?'
        if(nMXS.gt.0.and.rsample.gt.nMXS)stop'nMXS'
        if(nMXG.gt.0)then
          if(ict(0,2).gt.nMXG)stop'nMXG'
        end if
        if(nMXH.gt.0)then
          if(ict(0,4).gt.nMXH)stop'nMXH'
        end if
        if(nMXV.gt.0)then
          if(ict(0,6).gt.nMXV)stop'nMXV'
        end if
        if(nsample.lt.rsample)stop'nsample.lt.rsample'

        ict_sample = ict

        return
        do i=1,rsample
          if(id_proc.eq.0) &
          write(filenum,'(2i5,1x,a,2f8.3)')i,icone(i),inf(i),sampl(i,1),fun(i)
        end do
        end subroutine report_data
