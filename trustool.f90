        subroutine read_trustdx(ifunc)
        use dimKrig
        implicit none
        integer, intent(in) :: ifunc
        character(len=13) :: Cfile
        integer :: i,icl,icu
        double precision :: dxmin,dxmax

        if(itrust.eq.0)go to 1001

        if(ifunc.le.9)then
          write(Cfile,21)ifunc
        else if(ifunc.ge.10.and.ifunc.le.99)then
          write(Cfile,22)ifunc
        end if
21      format('trustdx0',i1,'.dat')
22      format('trustdx',i2,'.dat')
        if(tdxmin.ge.tdxmax)stop'tdxmin/max'

        dxmin =  1.d10
        dxmax = -1.d10
        icl   =  0
        icu   =  0
        open(10,file=Cfile,form='formatted',status='old',err=1000)
        do 100 i=1,nsample*2
          read(10,*,end=101)tdxx(i)
          if(tdxx(i).le.tdxmin)then
             tdxx(i) = 0.d0
             icl    = icl + 1
          else if(tdxx(i).ge.tdxmax)then
             tdxx(i) = dxmax
             icu    = icu + 1
          end if
          dxmin = min(dxmin,tdxx(i))
          dxmax = max(dxmax,tdxx(i))
100     continue
        stop'wrong trustdx.dat'
101     continue
        if(i.ne.nsample+1)then
          stop'wrong number of lines of trust.dx'
        end if
        close(10)

        if(id_proc.eq.0.and.ndebug.eq.1)then
          write(*,'(11x,a,2f12.7)')'>> Min/Max of Trust Region = ',dxmin,dxmax
          if(icl.ne.0.or.icu.ne.0)then
          write(*,'(11x,a,2i5)')'>> # of LU Modified Trust Region',icl,icu
          end if
        end if

        return
!+++++++++++++
1000    continue
        close(10)
1001    continue
        if(id_proc.eq.0)then
          write(*,'(11x,a,2f12.7)')'>> Set All dx_trust = ',(tdxinit(i),i=1,2)
        end if
        do i=1,nsample
           if(     info(i)(3:6).eq.'F   ')then
              tdxx(i) = 0.d0
           else if(info(i)(3:6).eq.'FG  ')then
              tdxx(i) = tdxinit(1)
           else if(info(i)(3:6).eq.'FH  ')then
              tdxx(i) = tdxinit(2)
           else if(info(i)(3:6).eq.'FGH ')then
              tdxx(i) = tdxinit(2)
           else
              stop'Hvec for indirect?'
           end if
        end do

        if(itrust.eq.1.and.id_proc.eq.0)then
          call write_trustdx(ifunc,-1.d0)
        end if
        
        end subroutine read_trustdx
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine write_trustdx(ifunc,dxnew)
        use dimKrig
        implicit none
        integer, intent(in) :: ifunc
        double precision, intent(in) :: dxnew
        character(len=13) :: Cfile
        integer :: i

        if(ifunc.le.9)then
          write(Cfile,21)ifunc
        else if(ifunc.ge.10.and.ifunc.le.99)then
          write(Cfile,22)ifunc
        end if
21      format('trustdx0',i1,'.dat')
22      format('trustdx',i2,'.dat')

        open(10,file=Cfile,form='formatted',status='unknown')
        do 100 i=1,nsample
           write(10,*)tdxx(i)
100     continue
        if(dxnew.ge.0.d0.and.dxnew.le.1.d0)then
        write(10,*)dxnew
        end if
        close(10)

        end subroutine write_trustdx
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine check_overlap(nlap,olap)
        use dimKrig
        implicit none
        integer, dimension(rsample),           intent(out) :: nlap
        integer, dimension(rsample,rsample,2), intent(out) :: olap
! olap(mine,nlap,1) = target
! olap(mine,nlap,2) = 1  : mine swallows target
! olap(mine,nlap,2) = 2  : mine is swallowed by target
! olap(mine,nlap,2) = 3  : deep
! olap(mine,nlap,2) = 4  : shallow
        integer :: i,j,k
        integer :: n1,n2
        integer :: ici,icd,ics
        double precision :: d1,d2,dd

        nlap =  0
        olap = -1

        ici = 0
        icd = 0
        ics = 0
        do 100 n1=1,rsample-1
        do 110 n2=n1+1,rsample
          d1 = tdx(n1)
          d2 = tdx(n2)
          dd = 0.d0
          do 120 k=1,ndim
            dd = dd + ( sampl(n1,k)-sampl(n2,k) )**2
120       continue
          dd = dsqrt(dd)

          if(d1.ge.(dd+d2))then              ! n1 swallows n2
              ici = ici + 1
              olap(n1,nlap(n1)+1,2) = 1
              olap(n2,nlap(n2)+1,2) = 2
!             if(ndebug.eq.1)write(*,'(a,2i5,3f12.7)') &
!             '*swallow ',n1,n2,d1,d2,dd
          else if(d2.ge.(dd+d1))then         ! n2 swallows n1
              ici = ici + 1
              olap(n1,nlap(n1)+1,2) = 2
              olap(n2,nlap(n2)+1,2) = 1
          else if(d1.ge.dd.or.d2.ge.dd)then  ! deep overlap
              icd = icd + 1
              olap(n1,nlap(n1)+1,2) = 3
              olap(n2,nlap(n2)+1,2) = 3
          else if(d1+d2.ge.dd)then           ! shallow overlap
              ics = ics + 1
              olap(n1,nlap(n1)+1,2) = 4
              olap(n2,nlap(n2)+1,2) = 4
          else if(dd.gt.d1+d2)then           ! may be one more inch
              go to 110
          else
              write(*,*)d1,d2,dd
              stop'unknown cases for overlap'
          end if
          nlap(n1) = nlap(n1) + 1
          olap(n1,nlap(n1),1) = n2
          nlap(n2) = nlap(n2) + 1
          olap(n2,nlap(n2),1) = n1
110     continue
100     continue

        if(id_proc.eq.0.and.Cmode(:4).eq.'Make')then
        write(*,'(11x,a,3i5)')'>> # of Overlap, (swal,deep,shlw) = ',ici,icd,ics
        end if

        do i=1,rsample
          if(nlap(i).ge.rsample)                           stop'nlap'
          do j=1,nlap(i)
             if(olap(i,j,1).lt.0.or.olap(i,j,1).gt.rsample)stop'olap-1'
             if(olap(i,j,2).lt.0.or.olap(i,j,2).gt.4)      stop'olap-2'
          end do
        end do

        end subroutine check_overlap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine reduce_additional_pts(nall,nlap,olap,rnew2)
        use dimKrig
        implicit none
! nall    : # of real/add sample points 
! rsample : # of real     sample points
! rnew2   : # of real/add sample points after elimination
        integer,                               intent(in) :: nall ! rnews
        integer, dimension(rsample),           intent(in) :: nlap
        integer, dimension(rsample,rsample,2), intent(in) :: olap
        integer,                               intent(out) :: rnew2

        integer, dimension(nall)       :: ndeath
        integer, dimension(nall,0:100) :: ndebt
        integer, dimension(nall)       :: tparent
        double precision, dimension(nall,ndim) :: tsampl
        double precision, dimension(nall)      :: tfun,tdxadd,fprev

        double precision, dimension(   1,   1) :: d1y,d2y
        double precision, dimension(   1,ndim) :: vdx
        double precision, dimension(ndim,   1) :: vg
        double precision, dimension(ndim,ndim) :: vh

        integer :: i,j,k
        integer :: ict
        integer :: n1,n2,m1,m2,m3
        double precision :: dd,d1,d2,tdd1,tdd2

        ndeath = 1 ! 0 means death
        ndebt  = 1
        do i=1,rsample
           ndebt(i,1) = i
        end do

! death by |dx|, bound
        ict = 0
        do i=rsample+1,nall
           if(nparent(i).eq.0)then
              ndeath(i) = 0 ! for |dx|=0
              ict = ict + 1
           end if
           ndebt(i,1) = nparent(i)
        end do
        if(id_proc.eq.0.and.Cmode(:4).eq.'Make')then
        write(*,'(11x,a,i5)') &
        '>> # of Eliminated Additional Pts (dx, bound) = ',ict
        end if

! deterministic death for swallow and deep overlap
        do 100 n1 = 1,rsample
          do 110 i=1,nlap(n1)
             n2    = olap(n1,i,1)
             m1    = olap(n1,i,2)
             if(     m1.eq.1)then
                m2 = 2
             else if(m1.eq.2)then
                m2 = 1
             else if(m1.eq.3)then
                m2 = 3
             else if(m1.eq.4)then
                m2 = 4
             else
                stop'unknown olap2'
             end if
             dd = 0.d0
             do k=1,ndim
                dd = dd + (sampl(n1,k)-sampl(n2,k))**2
             end do
             dd = dsqrt(dd)
             do 200 j=rsample+1,nall     ! additional sample loop
               if(ndeath(j).eq.0)go to 200
               if(nparent(j).ne.n1.and.nparent(j).ne.n2)go to 200
               d1 = 0.d0
               d2 = 0.d0
               do k=1,ndim
                 d1 = d1 + (sampl(n1,k)-sampl(j,k))**2
                 d2 = d2 + (sampl(n2,k)-sampl(j,k))**2
               end do
               d1 = dsqrt(d1)
               d2 = dsqrt(d2)
               if(     m1.eq.1.and.m2.eq.2)then        ! n1 swallows n2
                 if(nparent(j).eq.n2)then
                    ndeath(j) = 0
                    if(ndebug.eq.1)write(*,'(a,3i5,100f12.7)') &
                    '*Death by swallow',n1,n2,j,(sampl(j,k),k=1,ndim)
                 else if(nparent(j).eq.n1.and.inf(n2)(3:6).ne.'F   ')then
                    if(d2.le.d1)call add_ndebt(n2,j,rsample,nall,ndebt)
                 end if
               else if(m1.eq.2.and.m2.eq.1)then        ! n2 swallows n1
                 if(nparent(j).eq.n1)then
                    ndeath(j) = 0
                    if(ndebug.eq.1)write(*,'(a,3i5,100f12.7)') &
                    '*Death by swallow',n2,n1,j,(sampl(j,k),k=1,ndim)
                 else if(nparent(j).eq.n2.and.inf(n1)(3:6).ne.'F   ')then
                    if(d1.le.d2)call add_ndebt(n1,j,rsample,nall,ndebt)
                 end if
               else if(m1.eq.3.and.m2.eq.3)then        ! deep overlap
                 if(nparent(j).eq.n1)then
                    if(tdx(n2).ge.d2)then
                       ndeath(j) = 0
                       if(ndebug.eq.1)write(*,'(a,3i5,100f12.7)') &
                       '*Death by deep   ',n1,n2,j,(sampl(j,k),k=1,ndim)
                    end if
                 else if(nparent(j).eq.n2)then
                    if(tdx(n1).ge.d1)then
                       ndeath(j) = 0
                       if(ndebug.eq.1)write(*,'(a,3i5,100f12.7)') &
                       '*Death by deep   ',n2,n1,j,(sampl(j,k),k=1,ndim)
                    end if
                 else
                    stop'nparent?'
                 end if
               else if(m1.eq.4.and.m2.eq.4)then        ! shallow overlap
210              continue  ! deal after
               else
                 write(*,*)m1,m2
                 stop'unknown set of m1,m2'
               end if
200          continue  ! additional sample loop (j)
110       continue  ! nlap loop (i)
100     continue  ! real pts loop (n1)

        ict = 0
        do i=1,nall
          if(ndeath(i).eq.0)then
            if(i.le.rsample)stop'death of real pts ?'
            ict = ict + 1
          end if
        end do
        if(id_proc.eq.0.and.Cmode(:4).eq.'Make')then
        write(*,'(11x,a,i5)') &
        '>> # of Eliminated Additional Pts (swal,deep) = ',ict
        end if

! distance thresholds based death for shallow overlap and one more inch
        do 300 n1=rsample+1,nall
        do 310 n2=1,nall
          if(ndeath(n1).eq.0)go to 310
          if(ndeath(n2).eq.0)go to 310
          if(n1.eq.n2)go to 310
          dd = 0.d0
          do 320 k=1,ndim
            dd = dd + (sampl(n1,k)-sampl(n2,k))**2
320       continue
          dd = dsqrt(dd)

          if(n2.le.rsample)then        ! n1=additional,n2=real
            m1   = nparent(n1)
            if(m1.eq.n2)go to 310 ! parent(n2) and child(n1)
            if(tdx(n2).eq.0.d0)then ! substitute tdx(m1)
               tdd1 = max((tdx(m1)*tdd_thres(2)),tdd_thres(1))
            else
               tdd1 = max((tdx(n2)*tdd_thres(2)),tdd_thres(1))
            end if
            if(dd.ge.tdd1)go to 310
            ndeath(n1) = 0
            if(ndebug.eq.1)write(*,'(a,3i5,100f12.7)') &
            '*Death by r-a     ',n2,n1,m1,dd,tdd1,tdx(n2),tdx(m1)
            go to 310

          else if(n2.gt.rsample)then   ! n1=n2=additional
            m1 = nparent(n1)
            m2 = nparent(n2)
            if(m1.eq.m2)go to 310 ! n1 and n2 are brothers
            tdd1 = 0.5d0*(tdx(m1)+tdx(m2))

            tdd2 = max( (tdd1*tdd_thres(3)),tdd_thres(1) )
            do 330 i=1,nlap(m1)
               m3 = olap(m1,i,1)  ! m3=real, overlap with m1(n1's parent)
               if(m3.ne.m2)go to 330
               if(     olap(m1,i,2).eq.1)then        ! m1 swallows m2
                 stop'why m2s child still survives?'
               else if(olap(m1,i,2).eq.2)then        ! m2 swallows m1
                 stop'why m1s child still survives?'
               else if(olap(m1,i,2).eq.3.or.olap(m1,i,2).eq.4)then ! overlap
                 if(dd.ge.tdd2)go to 310
                 if(tdx(m1).le.tdx(m2))then
                   ndeath(n2) = 0
                   if(ndebug.eq.1)write(*,'(a,4i5,100f12.7)') &
                   '*Death by a-a over',n1,n2,m1,m2,dd,tdd1,tdd2
                   call add_ndebt(n2,n1,rsample,nall,ndebt)
                 else
                   ndeath(n1) = 0
                   if(ndebug.eq.1)write(*,'(a,4i5,100f12.7)') &
                   '*Death by a-a over',n2,n1,m2,m1,dd,tdd1,tdd2
                   call add_ndebt(n1,n2,rsample,nall,ndebt)
                 end if
                 go to 310
               else
                 stop'unknown olap2'
               end if
330         continue ! overlap loop

            tdd2 = max( (tdd1*tdd_thres(4)),tdd_thres(1) )
            if(dd.ge.tdd2)go to 310           ! one more inch
            if(tdx(m1).le.tdx(m2))then
              ndeath(n2) = 0
              call add_ndebt(n2,n1,rsample,nall,ndebt)
              if(ndebug.eq.1)write(*,'(a,4i5,100f12.7)') &
              '*Death by a-a inch',n1,n2,m1,m2,dd,tdd1,tdd2
            else
              ndeath(n1) = 0
              call add_ndebt(n1,n2,rsample,nall,ndebt)
              if(ndebug.eq.1)write(*,'(a,4i5,100f12.7)') &
              '*Death by a-a inch',n2,n1,m2,m1,dd,tdd1,tdd2
            end if
            go to 310
          end if ! n1,n2 are additional
310     continue ! all sample loop (n2)
300     continue ! additional sample loop (n1)

        ict = 0
        do i=1,nall
          if(ndeath(i).eq.0)then
            if(i.le.rsample)stop'death of real pts ?'
            ict = ict + 1
          end if
        end do
        if(id_proc.eq.0.and.Cmode(:4).eq.'Make') &
        write(*,'(11x,a,i5)') &
        '>> # of Eliminated Additional Pts ( finally ) = ',ict
        rnew2 = nall - ict

! Parameter Based Update of ndebt by dpar
        do 400 n1 = 1,rsample
          do 410 i=1,nlap(n1)
             n2 = olap(n1,i,1)
             m1 = olap(n1,i,2)
             dd = 0.d0
             do k=1,ndim
                dd = dd + (sampl(n1,k)-sampl(n2,k))**2
             end do
             dd = dsqrt(dd)
             do 420 j=rsample+1,nall
                if(ndeath(j).eq.0)go to 420
                if(nparent(j).ne.n1.and.nparent(j).ne.n2)go to 420
                d1 = 0.d0
                d2 = 0.d0
                do k=1,ndim
                   d1 = d1 + (sampl(n1,k)-sampl(j,k))**2
                   d2 = d2 + (sampl(n2,k)-sampl(j,k))**2
                end do
                d1 = dsqrt(d1)
                d2 = dsqrt(d2)
                if(m1.eq.3.or.m1.eq.4)then ! n1/n2(reals) are overlap
                  if(d1.le.0.d0.or.d2.le.0.d0)stop'd1,d2'
                  if(nparent(j).eq.n1)then
                    if( d2/d1.le.dpar )then
                      call add_ndebt(n2,j,rsample,nall,ndebt)
                    end if
                  else if(nparent(j).eq.n2)then
                    if( d1/d2.le.dpar )then
                      call add_ndebt(n1,j,rsample,nall,ndebt)
                    end if
                  else
                      stop'nparent at dpar'
                  end if
                end if
420          continue ! additional pt loop (j)
410       continue ! nlap loop (i)
400     continue ! real pt loop (n1)

! check ndebt
        ict = 0
        do 450 i=1,nall
          if(i.le.rsample)then
            if(ndebt(i,0).ne.1)stop'chc ndebt0 real'
            if(ndebt(i,1).ne.i)stop'chc ndebt1 real'
          else
            if(ndebt(i,0).lt.1.or.ndebt(i,0).gt.100)stop'chc ndebt0 add'
            do j=1,ndebt(i,0)-1
              do k=j+1,ndebt(i,0)
                 n1 = ndebt(i,j)
                 n2 = ndebt(i,k)
                 if(n1.lt.1.or.n1.gt.rsample)stop'chc ndebt1 real'
                 if(n2.lt.1.or.n2.gt.rsample)stop'chc ndebt1 real'
                 if(n1.eq.n2)stop'chc ndebt1 double'
              end do
            end do
            if(ndebt(i,0).gt.1.and.ndeath(i).ne.0)then
              ict = ict + 1
!             write(*,'(100f12.7)')(sampl(i,k),k=1,ndim)
            end if
          end if
450     continue
        if(id_proc.eq.0.and.Cmode(:4).eq.'Make') &
        write(*,'(11x,a,i5)') &
        '>> # of Additional Pts by Weighted Approx.    = ',ict

! weighting approximation
        fprev = fun
        call weighting_taylor(nall,ndeath,ndebt)

! Arrangement and output
        if(id_proc.eq.0)then
        open(20,file='idcok_real.dat',form='formatted',status='unknown')
        open(21,file='idcok_addt.dat',form='formatted',status='unknown')
        open(22,file='idcok_addw.dat',form='formatted',status='unknown')
        open(23,file='idcok_dead.dat',form='formatted',status='unknown')
25      format(100f15.8)
        end if
        tsampl = 0.d0
        tfun   = 0.d0
        tparent= 0
        tdxadd = 0.d0
        ict    = 0
        do 490 i=1,nall
          if(i.gt.rsample.and.nparent(i).eq.0)go to 490 ! |dx|=0

          if(id_proc.eq.0)then
          if(i.le.rsample)then
            write(20,25)(sampl(i,k),k=1,ndim),fun(i),tdx(i)
          else
            if(ndeath(i).eq.1)then
              if(ndebt(i,0).eq.1)then
                write(21,25)(sampl(i,k),k=1,ndim),fun(i)
              else if(ndebt(i,0).gt.1)then
                write(22,25)(sampl(i,k),k=1,ndim),fun(i),fprev(i)
              else
                stop'unknown ndebt0'
              end if
            else if(ndeath(i).eq.0)then
              write(23,25)(sampl(i,k),k=1,ndim),fun(i)
            else
              stop'ndeath'
            end if
          end if
          end if

          if(ndeath(i).eq.0)go to 490
          ict = ict + 1
          tsampl(ict,:) = sampl(i,:)
          tfun(ict)     = fun(i)
          tparent(ict)  = nparent(i)
          tdxadd(ict)   = dxadd(i)
490     continue
        if(id_proc.eq.0)then
        close(20)
        close(21)
        close(22)
        close(23)
        end if
        if(ict.ne.rnew2)then
           write(*,*)ict,rnew2,nall
           stop'ict.ne.rnew2 in elimination'
        end if
        sampl   = tsampl
        fun     = tfun
        nparent = tparent
        dxadd   = tdxadd

! Special Treatment for Current Optimal
        if(nOPT.ne.1)go to 500
        if(IDopt(1).le.0.or.IDopt(1).gt.nsample)go to 500
        n1 = 0
        do i=1,rsample
           if(icone(i).eq.IDopt(1))then
             n1 = i
           end if
        end do
        if(n1.le.0.or.n1.gt.rsample)stop'optimal was eliminated?'
        if(inf(n1)(3:6).eq.'F   '.or.inf(n1)(3:6).eq.'FG  ')go to 500
        if(it4opt.ne.1)go to 500
        if(tdx(n1).ge.tdxinit(1)*10.d0)then
          if(id_proc.eq.0)then
         open(24,file='idcok_spcl.dat',form='formatted',status='unknown')
          end if
          ict = 0
          vg(:,1) = gfun(n1,:)
          vh(:,:) = hfun(n1,:,:)
          do 510 i=1,1 !ndim
            if(rnew2+ict.lt.nall+1)then
              vdx      = tdxinit(1)/dsqrt(dble(ndim))
!             vdx      = 0.d0
!             vdx(1,i) = tdxinit(1)
              d1y      = matmul(vdx,vg)
              d2y      = matmul( matmul(vdx,vh),transpose(vdx) )
              ict      = ict + 1
              sampl(  rnew2+ict,:) = sampl(n1,:) + vdx(1,:)
              fun(    rnew2+ict  ) = fun(n1) + d1y(1,1) + 0.5d0*d2y(1,1)
              nparent(rnew2+ict  ) = -1 * n1 ! special case
              dxadd(  rnew2+ict  ) = tdxinit(1)

              do 520 j=1,rnew2 ! all real and additional
                 if(j.eq.n1)go to 520
                 dd = 0.d0
                 do k=1,ndim
                    dd = dd + (sampl(rnew2+ict,k)-sampl(j,k))**2
                 end do
                 dd = dsqrt(dd)
                 if(dd.le.tdd_thres(1))then
                    ict = ict - 1
                    go to 510
                 end if
520           continue
              if(id_proc.eq.0)then
              write(24,25)(sampl(rnew2+ict,k),k=1,ndim),fun(rnew2+ict)
              end if
            end if
510       continue
          if(ict.ne.0.and.id_proc.eq.0.and.Cmode(:4).eq.'Make')then
            write(*,'(11x,a,i3,a)') &
            '>> Special Treatment for Current Optimal, Add ',ict,' Pts'
          end if
          rnew2 = rnew2 + ict
          if(id_proc.eq.0)then
            close(24)
          end if
        end if
500     continue

        if(id_proc.eq.0)&
        write(*,'(6x,a,i5)') &
        '>> # of Final Sample Points for Indirect = ',rnew2
        end subroutine reduce_additional_pts
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine weighting_taylor(nall,ndeath,ndebt)
        use dimKrig
        implicit none
        integer, intent(in) :: nall
        integer, dimension(nall), intent(in) :: ndeath
        integer, dimension(nall,0:100), intent(in) :: ndebt

        double precision, dimension(1,1)        :: d1y,d2y
        double precision, dimension(1,ndim)     :: vdx
        double precision, dimension(ndim,1)     :: vg
        double precision, dimension(ndim,ndim)  :: vh
        integer :: i,j,k,nreal
        double precision :: yold,ynewd,ynewr,ynew
        double precision :: yj,rj,dj,dd,syr,syd,sr,sd,ww,ratio

        do 100 i=rsample+1,nall
           if(ndeath(i).eq.0)go to 100
           if(ndebt(i,0).le.0)stop'ndebt0-w'
           if(ndebt(i,0).eq.1)go to 100
           yold = fun(i)
           syr  = 0.d0
           syd  = 0.d0
           sr   = 0.d0
           sd   = 0.d0
           do 200 j=1,ndebt(i,0)
              nreal    = ndebt(i,j)
              if(nreal.le.0.or.nreal.gt.rsample)stop'ndebt1-w'
              if(tdx(nreal).le.0.d0)then
                 if(inf(nreal)(3:6).eq.'F   ')go to 200
                 write(*,'(i5,99f12.7)')    i,(sampl(    i,k),k=1,ndim)
                 write(*,'(i5,99f12.7)')nreal,(sampl(nreal,k),k=1,ndim)
                 stop'dx-w'
              end if
              if(inf(nreal)(3:6).eq.'F   ')then
                go to 200
              else if(inf(nreal)(3:6).eq.'FG  ')then
                if(j.eq.1)then
                  vg(:,1) = gfun(nreal,:)
                  vh      = 0.d0
                else
                  go to 200
                end if
              else if(inf(nreal)(3:6).eq.'FGH ')then
                vg(:,1) = gfun(nreal,:)
                vh(:,:) = hfun(nreal,:,:)
              else
                stop'unknown rfGRD in trust'
              end if
              vdx(1,:) = sampl(i,:) - sampl(nreal,:)
              d1y      = matmul(vdx,vg)
              d2y      = matmul( matmul(vdx,vh),transpose(vdx) )
              yj       = fun(nreal) + d1y(1,1) + 0.5d0*d2y(1,1)

              dd = 0.d0
              do k=1,ndim
                 dd = dd + vdx(1,k)**2
              end do
              dd  = dsqrt(dd)
              rj  = 1.d0 / tdx(nreal)
              dj  = 1.d0 / dd
              sr  = sr  + rj
              sd  = sd  + dj
              syr = syr + yj*rj
              syd = syd + yj*dj
200        continue
           if(sd.eq.0.d0.or.sr.eq.0.d0)stop'sd/sr'
           ynewd = syd / sd
           ynewr = syr / sr
           ww    = 1.d0
           ynew  = ww*ynewd + (1.d0-ww)*ynewr
           fun(i)= ynew

           if(yold.eq.0.d0)stop'yold in weight_taylor'
           ratio = abs( (ynew-yold)/(yold) )
           if(ratio.ge.0.2d0.and.id_proc.eq.0) &
           write(*,'(a,f12.7,2i5)') &
           '*Value much changed by weighted approx.',ratio,i,ndebt(i,0)
           if(ndebug.eq.1)then
             write(*,'(a,i3,a,10i5)') &
             '*Coorperation of ',ndebt(i,0),' Pts ',(ndebt(i,j),j=1,ndebt(i,0))
             write(*,'(1x,3(a,f12.7),100f12.7)') &
             'from',yold,' to',ynew,' @',(sampl(i,k),k=1,ndim)
           end if
100     continue

        end subroutine weighting_taylor
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine add_ndebt(ndead,nadpt,rsample,nall,ndebt)
        implicit none
! ndead (> rsample) : dead additional pt
! ndead (<=rsample) : swallowed real pt
! nadpt             : survived additional pt which gets new ndebt
        integer, intent(in) :: ndead,nadpt,rsample,nall
        integer, dimension(nall,0:100), intent(inout) :: ndebt
        integer :: i,j,m,n,nadd,nalr

        if(ndead.le.rsample)then !real pt specification
          if(ndebt(ndead,0).ne.1)stop'ndebt0 of real'
          if(ndebt(ndead,1).ne.ndead)stop'ndebt of real'
        end if

        m = ndebt(ndead,0)
        n = ndebt(nadpt,0)
        if(m.le.0.or.n.le.0)stop'ndebt0<=0'
        if(m.ge.100.or.n.ge.100)stop'many debts'

        do 100 i=1,m ! i=1 is parent of ndead
          nadd = ndebt(ndead,i)
          n    = ndebt(nadpt,0)
          do 110 j=1,n
            nalr = ndebt(nadpt,j)
            if(i.eq.1.and.j.eq.1)then
               if(nalr.eq.nadd)stop'brother killer?'
            end if
            if(nalr.eq.nadd)go to 100
110       continue
          if(n+1.gt.100)stop'ndebt0>100'
          ndebt(nadpt,0)   = n+1
          ndebt(nadpt,n+1) = nadd
100     continue

        end subroutine add_ndebt
