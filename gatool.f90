        subroutine check_param
        use dimGA
        use dimKrig, only:filenum
        implicit none
  
        if(id_proc.eq.0)then
         write(filenum,*)'>> Check Again the Parameters Setting'
! sharing
         if(Ctype_sh(:12).eq.'conventional'.or. &
            Ctype_sh(:12).eq.'Conventional') then
            write(filenum,'(6x,a)')'>> Use conventional sharing function'
         else if(Ctype_sh(:3).eq.'new'.or.Ctype_sh(:3).eq.'New') then
            write(filenum,'(6x,a)')'>> Use new sharing function (index-based)'
         end if
         if(ntype_sh.gt.0)then
            write(filenum,'(11x,a)') &
            '>> Use size of pool as N for sharing radius'
         else if(ntype_sh.lt.0)then
            write(filenum,'(11x,a)') &
            '>> Specify N for sharing radius'
         end if
         if(ngen_sh.eq.0)then
            write(filenum,'(11x,a)')'>> Use all design for sharing'
         else if(ngen_sh.gt.2)then
            write(filenum,'(11x,a,i3,a)') &
            '>> Use',ngen_sh,' recent generations for sharing'
         end if
         if(alpha_sh.ne.1.d0)then
            write(filenum,'(11x,a,f8.3)')'>> Alpha for sharing = ',alpha_sh
         end if
         if(Cpareto(:7).ne.'include'.and.Cpareto(:7).ne.'Include')then
            write(filenum,'(11x,a)')'>> Not include all pareto in pool'
         end if

! mating and selection
         if(Cprob(:3).eq.'new'.or.Cprob(:3).eq.'New')then
            write(filenum,'(6x,a)')'>> Use new fitness function'
         else if(Cprob(:3).eq.'old'.or.Cprob(:3).eq.'Old')then
            write(filenum,'(6x,a)')'>> Use old fitness function (1/rank)'
         end if
         if(fac_prob.le.0.d0.or.fac_prob.ge.1.d0)then
            write(filenum,'(11x,a)')'>> Power = 1.d0'
         else if(fac_prob.ne.0.2d0)then
            write(filenum,'(11x,a,f8.3)') &
            '>> Variable power with factor = ',fac_prob
         end if
         if(ngen_mat.eq.0)then
            write(filenum,'(11x,a)')'>> Use all design for mating'
         else if(ngen_mat.gt.2)then
            write(filenum,'(11x,a,i3,a)') &
            '>> Use',ngen_sh,' recent generations for mating'
         end if
         if(id_sus.le.2.or.id_sus.ge.10)then
            write(filenum,'(11x,a,i2)')'>> Number of division for SUS = ',id_sus
         end if
! crossover
         if(Cross(:3).eq.'BLX'.or.Cross(:3).eq.'blx')then
            write(filenum,'(6x,a)')'>> Use blended crossover'
         else if(Cross(:3).eq.'SBX'.or.Cross(:3).eq.'sbx')then
            write(filenum,'(6x,a)')'>> Use simulated binary crossover'
         end if
         if(if_crs.eq.1)then
            write(filenum,'(11x,a)')'>> Use only 1 random for a crossover'
         end if
         if(Cross(:3).eq.'BLX'.or.Cross(:3).eq.'blx')then
           if(beta_cr.ne.0.5d0)then
            write(filenum,'(11x,a,f8.3)')'>> BLX parameter = ',beta_cr
           end if
         else if(Cross(:3).eq.'SBX'.or.Cross(:3).eq.'sbx')then
           if(beta_cr.lt.1.0d0.or.beta_cr.gt.5.d0)then
            write(filenum,'(11x,a,f8.3)')'>> SBX parameter = ',beta_cr
           end if
         end if
         if(p_drxs.ne.0.d0.or.p_drxe.ne.0.d0)then
           write(filenum,'(11x,a,2f8.3)') &
           '>> Probability of directional crossover = ',p_drxs,p_drxe
         end if
         if(p_necs.ne.0.d0.or.p_nece.ne.0.d0)then
           write(filenum,'(11x,a,2f8.3)') &
           '>> Probability of neighborhood cultivation = ',p_necs,p_nece
         end if

! gradient
         if(ngen_grad.ne.ngen*2)then
           write(filenum,'(6x,a,i3,a)') &
           '>> Gradient-based evolution per',ngen_grad,' generations'
         end if

! mutation
         if(Cmut(:3).eq.'new'.or.Cmut(:3).eq.'New')then
            write(filenum,'(6x,a)')'>> Use new polynomial mutation'
         else if(Cmut(:3).eq.'old'.or.Cmut(:3).eq.'Old')then
            write(filenum,'(6x,a)')'>> Use old uniform mutation'
         end if
         if(p_muts.ne.0.1d0.or.p_mute.ne.0.1d0)then
            write(filenum,'(11x,a,2f8.3)') &
            '>> Probability of mutation = ',p_muts,p_mute
         end if
         if(r_mut.ne.0.1d0)then
            write(filenum,'(11x,a,2(f8.3,a))') &
            '>> Range of mutation = [',-1.d0*r_mut,' :',r_mut,' ]'
         end if
         if(Cmut(:3).eq.'new'.or.Cmut(:3).eq.'New')then
           if(eta_mut.lt.10.d0.or.eta_mut.gt.30.d0)then
             write(filenum,'(11x,a,f8.3)') &
             '>> Eta for mutation = ',eta_mut
           end if
         end if

! convergence
         if(ngen_stop.ne.0)then
            write(*,'(6x,a,i3)') &
            '>> Stop generation criteria      by ',ngen_stop
         end if
!        if(ngen_mut.ne.0)then
!           write(*,'(6x,a,i3,f8.3)') &
!           '>> Variable mutation probability by ',ngen_mut,fac_mut
!        end if
         if(neva_max.ne.0)then
            write(filenum,'(6x,a,i8)') &
            '>> Max number of evaluations      = ',neva_max
         end if
        end if

        end subroutine check_param


        subroutine convergence_check(gen,istop)
        use dimGA
        use dimKrig, only:filenum
        implicit none
        include 'mpif.h'
        integer, intent(in) :: gen
        integer, intent(out):: istop
        integer :: ierr
        integer :: i,j,g1,n1
        integer :: irank,last,lgen1,lgen2,ict
        double precision :: p_old
        character(len=14) :: Cfile1,Cfile2,Cfile3,Cfile4
! we have data until gen-1 generation
        istop = 0

! midterm report
        if(mod(gen-1,ngen_out).eq.0)then
          if(gen-1.le.9)then
             write(Cfile1,111)gen-1
!            write(Cfile2,121)gen-1
!            write(Cfile3,131)gen-1
!            write(Cfile4,141)gen-1
          else if(gen-1.ge.10.and.gen-1.le.99)then
             write(Cfile1,112)gen-1
!            write(Cfile2,122)gen-1
!            write(Cfile3,132)gen-1
!            write(Cfile4,142)gen-1
          else if(gen-1.ge.100.and.gen-1.le.999)then
             write(Cfile1,113)gen-1
!            write(Cfile2,123)gen-1
!            write(Cfile3,133)gen-1
!            write(Cfile4,143)gen-1
          else
            write(*,*)'*Cfile in convergence'
            call stop_all
          end if
111       format('_pareto00',i1,'.dat')
112       format('_pareto0',i2,'.dat')
113       format('_pareto',i3,'.dat')
121       format('_poolal00',i1,'.dat')
122       format('_poolal0',i2,'.dat')
123       format('_poolal',i3,'.dat')
131       format('_poolpr00',i1,'.dat')
132       format('_poolpr0',i2,'.dat')
133       format('_poolpr',i3,'.dat')
141       format('_poolsh00',i1,'.dat')
142       format('_poolsh0',i2,'.dat')
143       format('_poolsh',i3,'.dat')

          if(id_proc.eq.0)then
            open(11,file=Cfile1,form='formatted',status='unknown')
!           open(12,file=Cfile2,form='formatted',status='unknown')
!           open(13,file=Cfile3,form='formatted',status='unknown')
!           open(14,file=Cfile4,form='formatted',status='unknown')
            do g1=0,gen-1
            do n1=1,npop
               irank = int(rank(g1,n1))
               if(irank.eq.1)then
                 write(11,'(i6,999e15.5)')g1*npop+n1, &
                 (d(g1,n1,i),i=1,ndat),Frank(g1,n1),Rrank(g1,n1)
               end if
            end do 
            end do 
!           do j=1,ipool
!              g1 = pool(j,1)
!              n1 = pool(j,2)
!              write(12,'(i6,999e15.5)')g1*npop+n1, &
!              (d(g1,n1,i),i=1,ndat),Frank(g1,n1),Rrank(g1,n1)
!              irank = int(rank(g1,n1))
!              if(irank.eq.1)then
!                write(13,'(i6,999e15.5)')g1*npop+n1, &
!                (d(g1,n1,i),i=1,ndat),Frank(g1,n1),Rrank(g1,n1)
!              end if
!              if(Rrank(g1,n1).eq.1.d0)then
!                write(14,'(i6,999e15.5)')g1*npop+n1, &
!                (d(g1,n1,i),i=1,ndat),Frank(g1,n1),Rrank(g1,n1)
!              end if
!           end do
            close(11)
!           close(12)
!           close(13)
!           close(14)
          end if
        end if

! check the recent performance
        last = 0
        do i=1,nobj
           last = max(elite(i,1),last)
        end do
        lgen1 = gen - last
        if(nobj.eq.1) go to 100

        do g1=gen-1,0,-1
        do n1=1,npop
           irank = int(rank(g1,n1))
           if(irank.eq.1)then
              last = g1
              go to 100
           end if
        end do
        end do
        stop'no pareto optimal?'
100     continue
        lgen2 = gen - last

!       write(*,*)lgen1,lgen2
! lgen1 : number of generation recently elite designs have been uodated
! lgen2 : number of generation recently pareto designs have been uodated

! mutation prob
!       p_old = p_mut
!       if(lgen1.ge.ngen_mut.and.ngen_mut.gt.0)then
!          p_mut = p_mut * fac_mut
!          if(p_mut.ge.0.5d0)p_mut = 0.5d0
!       end if
!       if(lgen1.eq.1)then
!          p_mut = p_mut_kep  ! modify to original value
!       end if
!       if(p_old.ne.p_mut.and.id_proc.eq.0)then
!         write(*,'(11x,2(a,f7.3))') &
!         '>> Variable Prob. of Mutation',p_old,'->',p_mut
!       end if

! variable
        p_mut = p_muts + dble(gen-1)/dble(ngen-1)*(p_mute-p_muts)
        p_drx = p_drxs + dble(gen-1)/dble(ngen-1)*(p_drxe-p_drxs)
        p_nec = p_necs + dble(gen-1)/dble(ngen-1)*(p_nece-p_necs)
        if(id_proc.eq.0.and.ndebug.eq.1)then
          if(p_muts.ne.p_mute) &
          write(filenum,'(11x,a,f9.3)')'>> Prob. of Mutation      = ',p_mut
          if(p_drxs.ne.p_drxe) &
          write(filenum,'(11x,a,f9.3)')'>> Prob. of Direc. Cross. = ',p_drx
          if(p_necs.ne.p_nece) &
          write(filenum,'(11x,a,f9.3)')'>> Prob. of Neigh. Culti. = ',p_nec
        end if

! stop
        call MPI_Allreduce(ict_eva(1),ictg_eva(1),10,MPI_INTEGER, &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        ict = ictg_eva(1)+ictg_eva(2)+ictg_eva(3)
        if(lgen1.ge.ngen_stop.and.ngen_stop.gt.0)then
           istop = 1
           ngen_fin = gen-1
           if(id_proc.eq.0) &
           write(*,'(11x,a,i3,a)')'>> No Best Design within', &
           ngen_stop,' Recent Generations, STOP!!'
        end if
        if(ict.ge.neva_max.and.neva_max.gt.0)then
           istop = 1
           ngen_fin = gen-1
           if(id_proc.eq.0) &
           write(*,'(11x,a)') &
           '>> Number of Function Evaluations over the Threshold, &
               STOP!!'
        end if

        call MPI_BCAST(istop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        end subroutine convergence_check


        subroutine gradient_evolution(gen)
        use dimGA
        use dimKrig, only:filenum
        implicit none
        include 'mpif.h'
        integer, intent(in) :: gen
        integer, dimension(0:num_proc-1,1000)      :: ispl
        integer, dimension(0:num_proc-1)           :: ispm
        double precision, dimension(ndat)          :: xini
        double precision, dimension(nobj,ndv)      :: gini
        double precision, dimension(ndv)           :: xdir,xnew
        double precision, dimension(nobj,1)        :: wlam
        double precision, dimension(ndat)          :: ynew
        double precision, dimension(npop,ndat)     :: res

        integer :: i,j,k,l
        integer :: ista,iend,ierr,id,ict_l,ict_g
        double precision :: alpha,ratio,rat_l,rat_g
        if(id_proc.eq.0)then
          write(filenum,'(11x,a,i3,a)') &
          '>> Gradient-Based Evolution for ',gen,' Generation'
        end if

        ista = 1
        iend = npop

! decomposition of npop evaluation on nodes
        if(id_proc.eq.0) &
        call spliting(npop,num_proc,ista,ispl,ispm)
        call MPI_BCAST(ispl(0,1),num_proc*1000,MPI_INTEGER,0, &
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ispm(0),num_proc,MPI_INTEGER,0, &
                       MPI_COMM_WORLD,ierr)

! accord of parents
        call MPI_BCAST(parent(1,1,1),npop/2*3*ndat, &
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! gradient-based evolution
        res   = 0.d0
        rat_l = 0.d0
        ict_l = 0
        do i=1,ispm(id_proc)
           id = ispl(id_proc,i)
           if(id.lt.ista.or.id.gt.iend)stop'wrong id in ispl'
           j = id/2 + 1
           k = mod(id,2)
           if(k.eq.0)then
              j = j-1
              k = 2
           end if

           xini(:) = parent(j,k,:)
           call evagrad(xini,gini)
           call decide_direction(gini,wlam,xdir)
           call decide_stepsize(xini,gini,wlam,xdir,alpha,ratio,ynew)

           res(id,:) = ynew(:)
           if(alpha.gt.0.d0)then
           else if(alpha.eq.0.d0)then
               ratio = 0.d0
               res(id,:) = parent(j,k,:)
           else
               stop'negative alpha'
           end if

           do l=1,ndv
             if(res(id,l).lt.0.d0.or.res(id,l).gt.1.d0)then
               write(*,'(3e15.5)')xini(l),alpha,res(id,l)
               stop'res in grad'
             end if
           end do

           if(alpha.eq.0.d0)ict_l = ict_l + 1
           rat_l = rat_l + ratio
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(ict_l,ict_g,1,MPI_INTEGER, &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(rat_l,rat_g,1,MPI_DOUBLE_PRECISION, &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        rat_g = rat_g / dble(npop)
        if(id_proc.eq.0)then
          if(ict_g.gt.npop/10)  &
          write(filenum,'(11x,a,i6)') &
          '>> # of Designs without Improvement = ',ict_g
          write(filenum,'(11x,a,e10.3)') &
          '>> Average Improvement Ratio        = ',rat_g
        end if

! Sharing of Results
        do i=0,num_proc-1
          if(ispm(i).ne.0)then
            id = ispl(i,1)
            do j=1,ndat
              call MPI_BCAST( res(id,j), ispm(i), &
              MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
            end do
!           if(id_proc.eq.0) &
!           write(*,'(i3,100f6.2)')i,(res(j,1),j=ista,iend)
          end if
        end do

! Update dimension d[], evaluation was also finished !!
        do i=1,npop
          do j=1,ndat
             d(gen,i,j) = res(i,j)
          end do
        end do

        return
        if(id_proc.eq.0)then
          do i=1,npop/2
            do j=1,2
               write(106,'(999f15.8)')(parent(i,j,k),k=1,ndat)
            end do
          end do
          do i=1,npop
            write(107,'(999f15.8)')(d(gen,i,k),k=1,ndat)
          end do
          close(106)
          close(107)
        end if
        call stop_all

        end subroutine gradient_evolution
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine decide_stepsize(xini,gini,wlam,xdir,alpha,fwr,ynew)
        use dimGA
        use dimKrig, only:filenum
        implicit none
        double precision, dimension(ndat), intent(in)     :: xini
        double precision, dimension(nobj,ndv), intent(in) :: gini
        double precision, dimension(nobj,1), intent(in)   :: wlam
        double precision, dimension(ndv), intent(inout)   :: xdir
        double precision, intent(out)                     :: alpha,fwr
        double precision, dimension(ndat), intent(out)    :: ynew
        integer :: i,k,iter
        double precision :: det,alp,alpnew,alpini,alpmin,alpfac
        double precision :: fwnew,fwini,fwb,gwini
        double precision :: aa,bb,cc,dmax,dmin
        double precision, dimension(nobj) :: gw
        double precision, dimension(ndat) :: ymid
        double precision, dimension(0:100,3)  :: hist
        character(len=10) :: Ctag

        alpha = 0.d0
        ynew  = 0.d0
        fwr   = 0.d0

        alpini = 1.d-2  ! initial alpha
        alpmin = 1.d-10 ! lower threshold of alpha
        alpfac = 0.5d0  ! reduction factor of alpha
        Ctag = 'unknown'! debug tag

        det = 0.d0
        do i=1,ndv
           det = det + xdir(i)**2
        end do
        det = sqrt(det)
        if(det.le.0.d0)then
          write(filenum,'(a)')'*|xdir|=0 in decide_step'
          return
        end if
        xdir = xdir / det

        fwini = 0.d0
        gwini = 0.d0
        do i=1,nobj
           gw(i) = 0.d0
           do k=1,ndv
              gw(i) = gw(i) + gini(i,k)*xdir(k)
           end do
           ! fw is minimization function
           if(Cobj(i).eq.'Min'.or.Cobj(i).eq.'min')then
             fwini = fwini + wlam(i,1)*xini(ndv+i)
             gwini = gwini + wlam(i,1)*gw(i)
           else if(Cobj(i).eq.'Max'.or.Cobj(i).eq.'max')then
             fwini = fwini - wlam(i,1)*xini(ndv+i)
             gwini = gwini - wlam(i,1)*gw(i)
           else
             stop'unknown Cobj'
           end if
        end do

        iter   = 0
        hist(iter,1) = 0.d0
        hist(iter,2) = 0.d0
        hist(iter,3) = fwini/fwini
        alp    = alpini
100     continue
        iter = iter + 1
        call eva_1dsearch(xini,wlam,xdir,alp,alpnew,ymid,fwnew)
        hist(iter,1) = alp
        hist(iter,2) = alpnew
        hist(iter,3) = fwnew/fwini
        if(fwnew.lt.fwini.and.ymid(ndv+nobj+1).eq.0.d0)then
           fwb   = fwnew
           alpha = alp
           ynew  = ymid
        else
           if(alp.le.alpmin)then
             alpha = 0.d0
             fwb   = fwini
             if(iter.le.9)then
               write(Ctag,991)iter
             else
               write(Ctag,992)iter
             end if 
991          format('alpmin0',i1)
992          format('alpmin',i2)
             dmax = 0.d0
             dmin = 1.d0
             do i=1,ndv
                dmin = min(dmin,xini(i))
                dmax = max(dmax,xini(i))
             end do
!            write(*,'(i5,4f15.8)')iter,alp,alpnew,dmin,dmax
             go to 999 ! no improvement, return
           end if
           alp = alp * alpfac
           go to 100
        end if

        cc = fwini
        bb = gwini
        aa = (fwb-bb*alpha-cc)/(alpha**2)
!       if(aa.le.0.d0)write(*,*)'*aa<=0 in gradient'
        if(aa.gt.0.d0)then
          alp = (-1.d0*bb) / (2.d0*aa)
          if(alp.le.0.d0)then
             Ctag = 'aa>0,alp<0'
             go to 999
          end if

          iter = iter + 1
          call eva_1dsearch(xini,wlam,xdir,alp,alpnew,ymid,fwnew)
          hist(iter,1) = alp
          hist(iter,2) = alpnew
          hist(iter,3) = fwnew/fwini
          if(fwnew.lt.fwb.and.ymid(ndv+nobj+1).eq.0.d0)then
             fwb   = fwnew
             alpha = alp
             ynew  = ymid
             Ctag = 'aa>0,good1'
          else
             Ctag = 'aa>0,bad'
             alp  = 0.5d0*(alp+alpha)
             iter = iter + 1
             call eva_1dsearch(xini,wlam,xdir,alp,alpnew,ymid,fwnew)
             hist(iter,1) = alp
             hist(iter,2) = alpnew
             hist(iter,3) = fwnew/fwini
             if(fwnew.lt.fwb.and.ymid(ndv+nobj+1).eq.0.d0)then
               fwb   = fwnew
               alpha = alp
               ynew  = ymid
               Ctag = 'aa>0,good2'
             end if
          end if
        else
          Ctag = 'aa<0,bad'
          alp  = alpha * 1.5d0
          iter = iter + 1
          call eva_1dsearch(xini,wlam,xdir,alp,alpnew,ymid,fwnew)
          hist(iter,1) = alp
          hist(iter,2) = alpnew
          hist(iter,3) = fwnew/fwini
          if(fwnew.lt.fwb.and.ymid(ndv+nobj+1).eq.0.d0)then
             fwb   = fwnew
             alpha = alp
             ynew  = ymid
             Ctag = 'aa<0,good'
          end if
        end if

999     continue
        if(fwini.ne.0.d0)then
          fwr = (fwini-fwb)/dabs(fwini)
        else
          fwr = 0.d0
        end if
        if(ndebug.eq.1)then
        write(filenum,'(11x,a,i3,2e10.3,3x,a)') &
        '>> Improvement by Grad',iter,fwr,alpha,Ctag
!       do i=0,iter
!       write(108,'(i5,3e15.5)')i,hist(i,1),hist(i,2),hist(i,3)
!       end do
!       write(108,'(a)')' '
        end if

        end subroutine decide_stepsize
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine decide_direction(gini,wlam,xdir)
        use dimGA
        implicit none
        double precision, dimension(nobj,ndv), intent(in)  :: gini
        double precision, dimension(ndv), intent(out)      :: xdir
        double precision, dimension(nobj,1), intent(out)   :: wlam
        double precision, dimension(nobj,nobj)             :: A,Ainv,Achc
        double precision, dimension(nobj,1)                :: B
        integer :: i,j,k
        double precision :: det,chc

        if(nobj.eq.1)then ! mono-objective case
          det = 0.d0
          do i=1,ndv
             det = det + gini(1,i)**2
          end do
          det = sqrt(det)
          if(det.eq.0.d0)then
            do i=1,ndv
               xdir(i) =  1.d0
            end do
            wlam(1,1) = 1.d0
            return
!           stop'det=0 @ 1d case in gradient'
          end if
          if(Cobj(1).eq.'Min'.or.Cobj(1).eq.'min')then
            do i=1,ndv
               xdir(i) = -1.d0 * gini(1,i)/det
            end do
            wlam(1,1) = 1.d0
          else if(Cobj(1).eq.'Max'.or.Cobj(1).eq.'max')then
            do i=1,ndv
               xdir(i) = +1.d0 * gini(1,i)/det
            end do
            wlam(1,1) = 1.d0
          end if
          return
        end if

        A   = 0.d0
        do i=1,nobj
          do j=1,nobj
            if(Cobj(j).eq.'Min'.or.Cobj(j).eq.'min')then
              do k=1,ndv
               A(i,j) = A(i,j) -1.d0*(gini(i,k)*gini(j,k))
              end do
            else if(Cobj(j).eq.'Max'.or.Cobj(j).eq.'max')then
              do k=1,ndv
               A(i,j) = A(i,j) +1.d0*(gini(i,k)*gini(j,k))
              end do
            end if
          end do 
        end do

        B   = 0.d0
        det = 0.d0
        do i=1,nobj
          if(Cobj(i).eq.'Min'.or.Cobj(i).eq.'min')then
            if(Fmax(i)-Fmin(i).gt.0.d0)then
               B(i,1) = -1.d0 / (Fmax(i)-Fmin(i))
            else
               B(i,1) = -1.d0 / 0.0001d0
            end if
          else if(Cobj(i).eq.'Max'.or.Cobj(i).eq.'max')then
            if(Fmax(i)-Fmin(i).gt.0.d0)then
               B(i,1) = +1.d0 / (Fmax(i)-Fmin(i))
            else
               B(i,1) = +1.d0 / 0.0001d0
            end if
          end if
          det = det + B(i,1)**2
        end do
        det = sqrt(det)
        if(det.eq.0.d0)stop'det=0 @ multi-d case in gradient'
!       B = B / det

        if(nobj.eq.2)then
           call INV_MAT_2D(A,Ainv,det)
        else if(nobj.eq.3)then
           call INV_MAT_3D(A,Ainv,det)
        else
           stop'not considered yet, LU decomposition?'
        end if
        if(det.eq.0.d0)then
!         write(*,'(a)')'*Determinant = 0 for decision_direction'
          wlam = 1.d0
          go to 100
        end if

        Achc = matmul(Ainv,A)
        chc  = 0.d0
        do i=1,nobj
          do j=1,nobj
             if(i.eq.j)then
                chc = chc + (Achc(i,j)-1.d0)**2
             else
                chc = chc + (Achc(i,j)-0.d0)**2
             end if
          end do
        end do
        chc = sqrt(chc)
!       if(id_proc.eq.0.and.ndebug.eq.1)then
        if(chc.gt.1.d-5)then
          write(*,'(a,e15.5)')'*Error in A*A-1 = ',chc
        end if

        wlam = matmul(Ainv,B) ! nobjxnobj * nobjx1
 
        do i=1,nobj
          if(wlam(i,1).lt.0.d0)then
!           write(*,'(a,99e12.3)')'*wlam = ',(wlam(j,1),j=1,nobj)
            wlam(i,1) = 0.d0
!           wlam = 1.d0
!           go to 100
          end if
        end do
100     continue
        xdir = 0.d0
        do i=1,nobj
          if(Cobj(i).eq.'Min'.or.Cobj(i).eq.'min')then
            xdir(:) = xdir(:) - wlam(i,1) * gini(i,:)
          else if(Cobj(i).eq.'Max'.or.Cobj(i).eq.'max')then
            xdir(:) = xdir(:) + wlam(i,1) * gini(i,:)
          end if
        end do

        end subroutine decide_direction
        subroutine help
        use dimGA
        implicit none
        if(id_proc.eq.0)then
        write(*,'(1x,a)')'>> Welcome to Help of GA, &
                             Introduction of Input File Parameters'
        write(*,'(3x,a)')' 1: Dummy-Basic Parameters'
        write(*,'(3x,a)')' 2: number of population (even number)'
        write(*,'(3x,a)')' 3: number of max generation'
        write(*,'(3x,a)')' 4: number of design variables'
        write(*,'(3x,a)')' 5: number of objective functions, min/max'
        write(*,'(3x,a)')'    ex, "2 Min Max" means two objectives, &
                                   1st is minimize, 2nd is maximize'
        write(*,'(3x,a)')' 6: number of info in a design, >= ndv+nobj+1'
        write(*,'(3x,a)')' 7: kind of evaluation,  &
                              negative values are for test functions'
        write(*,'(3x,a)')' 8: random seed, (positive is active)'
        write(*,'(3x,a)')' 9: definition of initial population'
        write(*,'(3x,a)')'    "[Use[ [File Name]" : by a file'
        write(*,'(3x,a)')'     else               : by LHS sampling'
        write(*,'(3x,a)')' '
        write(*,'(3x,a)')'10: Dummy-Sharing'
        write(*,'(3x,a)')'11: type of sharing function, &
                              [Conventional] or [New]'
        write(*,'(3x,a)')'12: N for sharing radius, &
                              0:npareto, else+: pool size, &
                              else-: specify (rec=0)'
        write(*,'(3x,a)')'13: number of recent generations &
                              that are cosidered by sharing (rec=2)'
        write(*,'(3x,a)')'14: power for sharing function (rec=1)'
        write(*,'(3x,a)')'15: [include] all helpful pareto in the pool &
                              or not (rec=include)'
        write(*,'(3x,a)')' '
        write(*,'(3x,a)')'16: Dummy-Selection'
        write(*,'(3x,a)')'17: type of probability function, &
                              "Old" or "New"'
        write(*,'(3x,a)')'18: factor for the best design prob. function, &
                              (rec=0.0-0.2)'
        write(*,'(3x,a)')'19: number of recent generations that &
                              are used to make a mating pool (rec=2)'
        write(*,'(3x,a)')'20: number of division for SUS selection &
                              (rec=3-10)'
        write(*,'(3x,a)')' '
        write(*,'(3x,a)')'21: Dummy-Crossover'
        write(*,'(3x,a)')'22: type of crossover, "BLX" or "SBX"'
        write(*,'(3x,a)')'23: set within crossover (rec!=1)'
        write(*,'(3x,a)')'    1   : for a crossover,  a random  is used'
        write(*,'(3x,a)')'    else:                 ndv random are used'
        write(*,'(3x,a)')'24: parameter for crossover'
        write(*,'(3x,a)')'    BLX: 0.5 is recommended, &
                        <0.5: possibility of interpolation is bigger'
        write(*,'(32x,a)')'=0.5: possibilities are same'
        write(*,'(32x,a)')'>0.5: possibility of extrapolation is bigger'
        write(*,'( 3x,a)')'    SBX: [0;5] is recommended,     &
                                    possibilities are same at any value'
        write(*,'(12x,a)')'@ bigger value, child is closer to parent'
 
        write(*,'(12x,a)')'eta:         -1   0 .075 .3 .8    1    2    &
                              3    4    5   10  inf'
        write(*,'(12x,a)')'bet_in_ave:   0 0.5            0.66 0.74 &
                           0.79 0.82 0.85 0.90  1.0'
        write(*,'(12x,a)')'bet_ex_ave: inf   5    4  3  2 1.80 1.48 &
                           1.33 1.26 1.21 1.12  1.0'
        write(*,'(12x,a)')'when u(ran)=0.99, &
                           eta : 5.0 2.5 1.8 1.4 1.0 0.7'
        write(*,'(30x,a)')'bet : 2   3   4   5   7   10'
        write(*,'(3x,a)')'25: variable probability of &
                          directional crossover (rec=0.0-0.1)'
        write(*,'(3x,a)')'26: variable probability of &
                          neighborhood cultivation (rec~0.5)'
        write(*,'(3x,a)')' '
        write(*,'(3x,a)')'27: Dummy-Gradient-Based Evolution'
        write(*,'(3x,a)')'28: gradient evolution per ? generations'
        write(*,'(3x,a)')' '
        write(*,'(3x,a)')'29: Dummy-Mutation'
        write(*,'(3x,a)')'30: type of mutation, "New" or "Old"'
        write(*,'(3x,a)')'31: variable prob. of mutation, (rec~0.1)'
        write(*,'(3x,a)')'32: range of mutation [-r:r] (rec~0.1-0.2), &
                              eta only for "New" (rec=10-30)'
        write(*,'(7x,a)')'-1<eta<inf correspond to 1<dxave<0.02'
        write(*,'(7x,a)')'when eta=0, distribution is a straight line'
        write(*,'(7x,a)')'eta    :-1,  0,0.5,1.4,3.3,5.2,9.7,&
                                  11,14,17,21,29,44,89,inf'
        write(*,'(7x,a)')'dx_ave : 1,0.5,0.4,0.3,0.2,.15,0.1,&
                                  09,08,07,06,05,04,03,0.02'
        write(*,'(3x,a)')' '
        write(*,'(3x,a)')'33: Dummy-Covergence'
        write(*,'(3x,a)')'34: stopper when best design is &
                          not obtained recent ? generations'
        write(*,'(7x,a)')'0 is non-active'
!       write(*,'(3x,a)')'35: variable set of mutation prob., igen,fac'
!       write(*,'(7x,a)')'when best design is not obtained recent &
!                         ? generations,'
!       write(*,'(7x,a)')'mutation prob. is increased by the factor'
        write(*,'(3x,a)')'35: number of maximum function evaluations'
        write(*,'(3x,a)')'36: number of generation for output'
        write(*,'(3x,a)')'37: debug mode when 1'

        end if
        call stop_all
        end subroutine help
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine make_initial_pop
        use dimGA
        use dimKrig, only:filenum
        implicit none
        integer :: leng,i,j,nsample
        double precision, allocatable, dimension(:,:) :: bound,sample

        nsample = npop
        if(Cini.eq.'Use'.or.Cini.eq.'use')then
          leng = index(Cinifile,'  ') - 1
          if(id_proc.eq.0) &
          write(*,'(2a)')' >> Read Initial Population from ',Cinifile(:leng)
          go to 1001
        else
          if(id_proc.eq.0) &
          write(filenum,'(1a)')' >> Make Initial Population by LHS'
          go to 1002
        end if

1001    continue
        open(10,file=Cinifile,form='formatted',status='old',err=1002)
        do i=1,npop
          read(10,*,end=100,err=100)(d(0,i,j),j=1,ndv)
        end do
        close(10)
        return    ! all designs were given by Cinifile
100     continue
        close(10)
        nsample = npop - (i-1)
        if(id_proc.eq.0) &
        write(*,'(6x,a,i3,a)') &
        '>> Remaining',nsample,' designs are given by LHS'

1002    continue
        allocate(sample(ndv,nsample))
        allocate(bound(2,ndv))
        bound(1,:) = 0.d0
        bound(2,:) = 1.d0
        sample     = 0.d0
        call Latin_HypercubeGA(ndv,nsample,bound,sample)
        do i=1,nsample
           do j=1,ndv
              d(0,npop+1-i,j) = sample(j,i)
           end do
        end do
        deallocate(bound,sample)
        return

        do i=1,npop
           if(id_proc.eq.0) &
           write(101,'(99f15.8)')(d(0,i,j),j=1,ndv)
        end do
        end subroutine make_initial_pop
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine Latin_HypercubeGA(ndim,nsam,bound,sample)
        implicit none
        integer, intent(in) :: ndim,nsam
        double precision, dimension(2,ndim), intent(in) :: bound
        double precision, dimension(ndim,nsam), intent(out) :: sample

        integer :: i,j,k,l,key
        double precision :: ran,dv01
        double precision, dimension(0:nsam) :: split
        integer, dimension(ndim,nsam) :: iflag_split

        if(nsam.eq.0)return
        if(nsam.lt.0.or.ndim.le.0)stop'LHS dimension/nsample'
        sample = 0.d0

! Set Cartesian
        do 100 i=0,nsam
           split(i) = dble(i)/dble(nsam)
100     continue

! Arrangement
        do 200 i=1,nsam
           do 210 j=1,ndim
250           continue
              call random_number(ran)
              do 220 k=1,nsam
                 if(ran.ge.split(k-1).and.ran.le.split(k))then
                   do 230 l=1,i-1
                      if(iflag_split(j,l).eq.k)go to 250 ! already used, try again
230                continue    ! previous sample loop
                   iflag_split(j,i) = k
                   go to 210
                 end if
220           continue         ! split loop
              stop'out of order at split loop'
210        continue            ! dimension loop
200     continue               ! sample loop
        
! Check
        do 300 j=1,ndim
          do 310 i=1,nsam
             if(iflag_split(j,i).le.0.or.iflag_split(j,i).gt.nsam) &
             stop'iflag_split'
             do 320 k=i+1,nsam
                if(iflag_split(j,i).eq.iflag_split(j,k)) &
                stop'overlap of iflag_split'
320          continue
310       continue
300     continue

! Generation
        do 400 i=1,nsam
          do 410 j=1,ndim
             key  = iflag_split(j,i)
             call random_number(ran)      ! random set
             ran = 0.01d0 + 0.99d0*ran
!            ran  = 0.5d0                 ! center of cartesian
             dv01 =  split(key-1) + ran*(split(key)-split(key-1))
             sample(j,i) = bound(1,j) + (bound(2,j)-bound(1,j))*dv01
410       continue
400     continue

      end subroutine Latin_HypercubeGA



      subroutine INV_MAT_2D(A,AINV,det)
      implicit none
      real*8, dimension(2,2), intent(in)  :: A
      real*8, dimension(2,2), intent(out) :: AINV
      real*8,                 intent(out) :: det
      real*8 :: a11,a21,a12,a22

      a11 = A(1,1)
      a21 = A(2,1)
      a12 = A(1,2)
      a22 = A(2,2)
      det = a11*a22 - a12*a21
      if(det.eq.0.d0)return  !stop'det = 0, INV_MAT_2D'
      det = 1.d0/det

      AINV(1,1) = det * a22
      AINV(1,2) = det *(-1.d0*a12)
      AINV(2,1) = det *(-1.d0*a21)
      AINV(2,2) = det * a11
      det = 1.d0/det
      end subroutine INV_MAT_2D
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine INV_MAT_3D(A,AINV,det)
      implicit none
      real*8, dimension(3,3), intent(in)  :: A
      real*8, dimension(3,3), intent(out) :: AINV
      real*8,                 intent(out) :: det
      real*8 :: a11,a21,a31,a12,a22,a32,a13,a23,a33

      a11 = A(1,1)
      a21 = A(2,1)
      a31 = A(3,1)
      a12 = A(1,2)
      a22 = A(2,2)
      a32 = A(3,2)
      a13 = A(1,3)
      a23 = A(2,3)
      a33 = A(3,3)

      det = a11*a22*a33 + a21*a32*a13 + a31*a12*a23 &
          - a11*a32*a23 - a31*a22*a13 - a21*a12*a33
      if(det.eq.0.d0)return !stop'det = 0, INV_MAT_3D'
      det = 1.d0 / det
      AINV(1,1) = det * (a22*a33 - a23*a32)
      AINV(1,2) = det * (a13*a32 - a12*a33)
      AINV(1,3) = det * (a12*a23 - a13*a22)
      AINV(2,1) = det * (a23*a31 - a21*a33)
      AINV(2,2) = det * (a11*a33 - a13*a31)
      AINV(2,3) = det * (a13*a21 - a11*a23)
      AINV(3,1) = det * (a21*a32 - a22*a31)
      AINV(3,2) = det * (a12*a31 - a11*a32)
      AINV(3,3) = det * (a11*a22 - a12*a21)
      det = 1.d0 / det
      end subroutine INV_MAT_3D
        subroutine output
        use dimGA
        use dimKrig, only:filenum
        implicit none
        include 'mpif.h'
        integer :: ierr,i,j,n1,g1,irank

        call MPI_Allreduce(ict_eva(1),ictg_eva(1),10,MPI_INTEGER, &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        if(id_proc.eq.0)then
          write(filenum,'(a)')' >> Output Results....'
          open(11,file='_alldes.dat',form='formatted',status='unknown')
          open(12,file='_pareto.dat',form='formatted',status='unknown')
          do g1=0,ngen_fin
          do n1=1,npop
             irank = int(rank(g1,n1))
             write(11,10)g1*npop+n1, &
             (d(g1,n1,i),i=1,ndat),rank(g1,n1),Frank(g1,n1),Rrank(g1,n1)
             if(irank.eq.1)then
               write(12,10)g1*npop+n1, &
               (d(g1,n1,i),i=1,ndat),rank(g1,n1),Frank(g1,n1),Rrank(g1,n1)
             end if
          end do
          end do
          close(11)
          close(12)
10        format(i6,9999e20.10)

          open(14,file='_elite.dat',form='formatted',status='unknown')
          do i=1,nobj
             g1=elite(i,1)
             n1=elite(i,2)
             write(14,10)g1*npop+n1, &
             (d(g1,n1,j),j=1,ndat),rank(g1,n1),Frank(g1,n1),Rrank(g1,n1)
             write(filenum,'(6x,a,i2,a,2i4,a,10e12.4)') &
             '>> Best Design for ',i,' : (',g1,n1,') :', &
             (d(g1,n1,j),j=ndv+1,ndv+nobj)
          end do
          close(14)
          do i=1,10
            if(i.eq.1)then
              write(filenum,'(6x,a,i8)') &
              '>> # of Function Evaluations for New Ones  = ',ictg_eva(i)
            else if(i.eq.2)then
              write(filenum,'(6x,a,i8)') &
              '>> # of Function Evaluations for Gradient  = ',ictg_eva(i)
            else if(i.eq.3)then
              write(filenum,'(6x,a,i8)') &
              '>> # of Function Evaluations for 1D Search = ',ictg_eva(i)
            else
              if(ictg_eva(i).ne.0)write(filenum,*)'*ict>0 for i=',i
            end if
          end do
          write(filenum,'(6x,a,i8)') &
          '>> # of Function Evaluations in Total      = ', &
          ictg_eva(1)+ictg_eva(2)+ictg_eva(3)
        end if

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        end subroutine output



        subroutine mating(gen)
        use dimGA
        use dimKrig, only:filenum
        implicit none
        include 'mpif.h'
        integer, intent(in) :: gen
        integer, dimension((ngen+1)*npop)      :: sort
        integer, dimension((ngen+1)*npop,nobj) :: inec
        double precision, dimension(0:(ngen+1)*npop) :: prob
        integer :: g1,n1,g2,n2,g3,n3,id,id1,id2,id3,nec,nec1,nec2,nec3
        integer :: i,j,k,ierr,itmp,iter
        integer :: ict_mat,igen,ict_hlf,ict,irank
        integer :: itarget,ipmax,mp
        double precision :: FM,RM,psum,ran,power,f2,f3
! we have data until gen-1 generation

! make pool by ngen_mat
        if(ngen_mat.eq.0)then
           igen = 0
        else if(ngen_mat.gt.0)then
           igen = (gen-1) - (ngen_mat-1)
           if(igen.lt.0)igen = 0
        else
           stop'negative ngen_mat'
        end if
        ict_mat = 0
        pool    = 0
        FM      = 0.d0
        RM      = 0.d0
        do 100 g1=igen,gen-1   ! order is not important here
        do 100 n1=1,npop
           ict_mat = ict_mat + 1
           pool(ict_mat,1) = g1
           pool(ict_mat,2) = n1
           FM = max(Frank(g1,n1),FM)
           RM = max(Rrank(g1,n1),RM)
           if(Frank(g1,n1).lt.0.d0)stop'Frank.lt.0'
           if(Rrank(g1,n1).lt.0.d0)stop'Rrank.lt.0'
100     continue
! add elite
        ict_hlf = ict_mat
        do 200 i=1,nobj
           g1 = elite(i,1)
           n1 = elite(i,2)
           if(g1.lt.igen)then
             do 210 j=1,i-1
                g2 = elite(j,1)
                n2 = elite(j,2)
                if(g1.eq.g2.and.n1.eq.n2)go to 200
210          continue
             ict_mat = ict_mat + 1
             pool(ict_mat,1) = g1
             pool(ict_mat,2) = n1
             Frank(g1,n1) = FM
             Rrank(g1,n1) = RM
             if(Frank(g1,n1).lt.0.d0)stop'Frank.lt.0e'
             if(Rrank(g1,n1).lt.0.d0)stop'Rrank.lt.0e'
           end if
200     continue
! add pareto
        if(Cpareto(:7).eq.'include'.or.Cpareto(:7).eq.'Include')then
          do 300 g1=0,igen-1
          do 300 n1=1,npop
             irank = int(rank(g1,n1))
             if(irank.ne.1)go to 300
             if(Rrank(g1,n1).ne.1.d0)go to 300
             do 310 i=1,nobj
                g2 = elite(i,1)
                n2 = elite(i,2)
                if(g1.eq.g2.and.n1.eq.n2)go to 300
310          continue
             ict_mat = ict_mat + 1
             pool(ict_mat,1) = g1
             pool(ict_mat,2) = n1
             Frank(g1,n1) = FM
             Rrank(g1,n1) = RM
300       continue
        end if

        ipool = ict_mat
        if(ict_hlf.ne.ict_mat.and.id_proc.eq.0.and.ndebug.eq.1)then
           write(filenum,'(11x,a,i3,a)') &
           '>> Add ',(ict_mat-ict_hlf),' Elite Points in Mating Pool'
        end if
!       if(id_proc.eq.0)then
!       do i=1,ipool
!        write(*,'(a,3i5)')'pool',i,pool(i,1),pool(i,2)
!       end do
!       end if

! selection for gradient-based evolution (only by sort)
        if(mod(gen,ngen_grad).eq.0)then
          do i=1,ict_mat
            sort(i) = i
          end do
          do i=1,ict_mat-1
            do j=i+1,ict_mat
              g1 = pool(sort(i),1)
              n1 = pool(sort(i),2)
              g2 = pool(sort(j),1)
              n2 = pool(sort(j),2)
              if(Frank(g1,n1).lt.Frank(g2,n2))then
                 itmp    = sort(i)
                 sort(i) = sort(j)
                 sort(j) = itmp
              end if
            end do
          end do
          do i=1,ict_mat-1
            g1 = pool(sort(i),1)
            n1 = pool(sort(i),2)
            g2 = pool(sort(i+1),1)
            n2 = pool(sort(i+1),2)
            if(Frank(g1,n1).lt.Frank(g2,n2))then
              stop'sorting in pool'
            end if
          end do
          ict = 0
          parent  = -1.d0
          do i=1,npop/2
            do j=1,2
               ict = ict + 1
               if(ict.gt.ict_mat)stop'ict in sort_pool'
               g1  = pool(sort(ict),1)
               n1  = pool(sort(ict),2)
               parent(i,j,:)          = d(g1,n1,:)
            end do
          end do
!         write(*,'(11x,2(a,i5))') &
!         '>> Choose',npop,' pts from pool of',ict_mat
          return
        end if

! sort for neighborhood
        if(p_nec.ne.0.d0)then
          do i=1,ict_mat
            inec(i,:) = i
          end do
          do k=1,nobj
           do i=1,ict_mat-1
            do j=i+1,ict_mat
              g1 = pool(inec(i,k),1)
              n1 = pool(inec(i,k),2)
              g2 = pool(inec(j,k),1)
              n2 = pool(inec(j,k),2)
              if(d(g1,n1,ndv+k).gt.d(g2,n2,ndv+k))then
                 itmp      = inec(i,k)
                 inec(i,k) = inec(j,k)
                 inec(j,k) = itmp
              end if
            end do
           end do
           do i=1,ict_mat-1
            g1 = pool(inec(i,  k),1)
            n1 = pool(inec(i,  k),2)
            g2 = pool(inec(i+1,k),1)
            n2 = pool(inec(i+1,k),2)
            if(d(g1,n1,ndv+k).gt.d(g2,n2,ndv+k))then
              stop'sorting in pool for nec'
            end if
           end do
          end do
        end if

! power for probability
        if(fac_prob.le.0.d0.or.fac_prob.ge.1.d0)then
          power = 1.d0
        else
          if(id_proc.eq.0) &
          call decide_power(FM,RM,power)
          call MPI_BCAST(power,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        end if

! probability
        psum = 0.d0
        do i=1,ict_mat
           g1 = pool(i,1)
           n1 = pool(i,2)
           if(Cprob(:3).eq.'old'.or.Cprob(:3).eq.'Old')then
             psum = psum + Rrank(g1,n1)**power
           else if(Cprob(:3).eq.'new'.or.Cprob(:3).eq.'New')then
             psum = psum + Frank(g1,n1)**power
           end if
        end do
        if(psum.le.0.d0)stop'psum<=0'
        prob = 0.d0
        do i=1,ict_mat
           g1 = pool(i,1)
           n1 = pool(i,2)
           if(Cprob(:3).eq.'old'.or.Cprob(:3).eq.'Old')then
             prob(i) = prob(i-1) + (Rrank(g1,n1)**power)/psum
           else if(Cprob(:3).eq.'new'.or.Cprob(:3).eq.'New')then
             prob(i) = prob(i-1) + (Frank(g1,n1)**power)/psum
           end if
        end do
        if(prob(ict_mat).lt.0.99d0.or.prob(ict_mat).gt.1.01d0)then
          write(*,'(a,2f10.5)')'*psum,prob(imax) = ',psum,prob(ict_mat)
          stop'sum of probability?'
        end if
        do i=ict_mat,(ngen+1)*npop
          prob(i) = 1.d0
        end do
        if(id_proc.eq.0.and.ndebug.eq.1)then
           if(Cprob(:3).eq.'old'.or.Cprob(:3).eq.'Old')then
           write(filenum,'(11x,a,a3,a,f8.3,a,f8.3)') &
           '>> Use ',Cprob,' Prob. Func., Pmax = ',(RM**power)/psum, &
           ' with power = ',power
           else if(Cprob(:3).eq.'new'.or.Cprob(:3).eq.'New')then
           write(filenum,'(11x,a,a3,a,f8.3,a,f8.3)') &
           '>> Use ',Cprob,' Prob. Func., Pmax = ',(FM**power)/psum, &
           ' with power = ',power
           end if
        end if

! selection of npop parents by SUS
        itarget =  1
        parent  = -1.d0
        do 500 i=1,npop/2
           id3 = -1
           g3  = -1
           n3  = -1

           ! farther
           call SUS(itarget,id_sus,ict_mat,(ngen+1)*npop,prob,id1)
           if(id1.le.0.or.id1.gt.ict_mat)stop'id1 by SUS'
           g1 = pool(id1,1)
           n1 = pool(id1,2)
           pool(id1,3) = pool(id1,3) + 1
           parent(i,1,:)  = d(g1,n1,:)  ! farther
           parent(i,1,ndv+nobj+1) = Frank(g1,n1)

           ! neighborhood cultivation
           if(p_nec.ne.0.d0)then
             call random_number(ran)
             if(ran.lt.p_nec*0.5d0.or.ran.gt.1.d0-p_nec*0.5d0)then
               iaxs_nec = iaxs_nec + 1
               if(iaxs_nec.lt.   1)iaxs_nec = 1
               if(iaxs_nec.gt.nobj)iaxs_nec = 1  ! [1:nobj]
               do j=1,ict_mat
                 if(inec(j,iaxs_nec).eq.id1)then
                   nec1 = j
                   go to 510
                 end if
               end do
               stop'no target for neighborhood cultivation'
510            continue
               if(ict_mat.le.20)stop'modify the selection for nec'
               nec2 = 0
               nec3 = 0
               f2   = 0.d0
               f3   = 0.d0
               do 520 j=int(-0.10d0*dble(ict_mat)),int(0.10d0*dble(ict_mat))
                  nec = nec1 + j
                  if(nec.eq.nec1)go to 520
                  if(nec.le.0.or.nec.gt.ict_mat)go to 520
                  id  = inec(nec,iaxs_nec)
                  g2  = pool(id,1)
                  n2  = pool(id,2)
                  if(Frank(g2,n2).ge.f2)then
                    f3   = f2
                    nec3 = nec2
                    f2   = Frank(g2,n2)
                    nec2 = nec
                  else if(Frank(g2,n2).lt.f2.and.Frank(g2,n2).ge.f3)then
                    f3   = Frank(g2,n2)
                    nec3 = nec
                  end if
!                 write(*,*)g2,n2,nec1,nec2,nec3
520            continue
               if(nec2.eq.0.or.nec3.eq.0)stop'nec2/nec3'
               if(nec2.eq.nec3)stop'nec2=nec3'

               id2 = inec(nec2,iaxs_nec)
               g2  = pool(id2,1)
               n2  = pool(id2,2)
               pool(id2,3) = pool(id2,3) + 1
               parent(i,2,:)  = d(g2,n2,:)  ! mother
               parent(i,2,ndv+nobj+1) = Frank(g2,n2)
               ! neighbor && directional
               call random_number(ran)
               if(ran.lt.p_drx*0.5d0.or.ran.gt.1.d0-p_drx*0.5d0)then
                 id3 = inec(nec3,iaxs_nec)
                 g3  = pool(id3,1)
                 n3  = pool(id3,2)
                 pool(id3,3) = pool(id3,3) + 1
                 parent(i,3,:)  = d(g3,n3,:)  ! flirt
                 parent(i,3,ndv+nobj+1) = Frank(g3,n3)
               end if
                  
!               call random_number(ran)  ! ict_mat x [-0.05:0.05]
!               nec2 = nec1 + int( dble(ict_mat)*(ran-0.5d0)*0.1d0 )
!               if(nec2.le.0      )nec2 = 1
!               if(nec2.gt.ict_mat)nec2 = ict_mat
!               if(nec2.eq.nec1)go to 533
!               id2 = inec(nec2,iaxs_nec)
!               g2  = pool(id2,1)
!               n2  = pool(id2,2)
!               if(g1.eq.g2.and.n1.eq.n2)stop'gn_mother for nec'
!               pool(id2,3) = pool(id2,3) + 1
!               parent(i,2,:)  = d(g2,n2,:)  ! mother
!               parent(i,2,ndv+nobj+1) = Frank(g2,n2)
!               ! neighbor && directional
!               call random_number(ran)
!               if(ran.lt.p_drx*0.5d0.or.ran.gt.1.d0-p_drx*0.5d0)then
!544              continue
!                 call random_number(ran)
!                 nec3 = nec1 + int( dble(ict_mat)*(ran-0.5d0)*0.1d0 )
!                 if(nec3.le.0      )nec3 = 1
!                 if(nec3.gt.ict_mat)nec3 = ict_mat
!                 if(nec3.eq.nec1.or.nec3.eq.nec2)go to 544
!                 id3 = inec(nec3,iaxs_nec)
!                 g3  = pool(id3,1)
!                 n3  = pool(id3,2)
!                 if(g1.eq.g3.and.n1.eq.n3)stop'gn_flirt for nec'
!                 if(g2.eq.g3.and.n2.eq.n3)stop'gn_flirt for nec'
!                 pool(id3,3) = pool(id3,3) + 1
!                 parent(i,3,:)  = d(g3,n3,:)  ! flirt
!                 parent(i,3,ndv+nobj+1) = Frank(g3,n3)
!               end if
               go to 599
             end if
           end if

           ! mother in usual
           iter = 0
555        continue
           iter = iter + 1
           call SUS(itarget,id_sus,ict_mat,(ngen+1)*npop,prob,id2)
           if(id2.le.0.or.id2.gt.ict_mat)stop'id2 by SUS'
           g2 = pool(id2,1)
           n2 = pool(id2,2)
           pool(id2,3) = pool(id2,3) + 1
           if(g1.eq.g2.and.n1.eq.n2)then
             if(iter.le.100)go to 555 ! farther is mother
             g2 = g1
             n2 = n1+1
             if(n2.gt.npop)n2=1
           end if
           parent(i,2,:)  = d(g2,n2,:)  ! mother
           parent(i,2,ndv+nobj+1) = Frank(g2,n2)

           ! directional in usual
           call random_number(ran)
           if(ran.lt.p_drx*0.5d0.or.ran.gt.1.d0-p_drx*0.5d0)then
             iter = 0
566          continue
             iter = iter + 1
             call SUS(itarget,id_sus,ict_mat,(ngen+1)*npop,prob,id3)
             if(id3.le.0.or.id3.gt.ict_mat)stop'id3 by SUS'
             g3 = pool(id3,1)
             n3 = pool(id3,2)
             pool(id3,3) = pool(id3,3) + 1
             if( (g1.eq.g3.and.n1.eq.n3).or.(g2.eq.g3.and.n2.eq.n3) )then
               if(iter.le.100)go to 566 
               g3 = g2
               n3 = n2+1
               if(n3.gt.npop)n3=1
             end if
             parent(i,3,:)  = d(g3,n3,:) ! flirt
             parent(i,3,ndv+nobj+1) = Frank(g3,n3)
           end if

599        continue
           if((g1.eq.g2.and.n1.eq.n2).or.(g1.eq.g3.and.n1.eq.n3).or. &
              (g2.eq.g3.and.n2.eq.n3))then
              write(*,'(a,3i5)')'G123 = ',g1,g2,g3
              write(*,'(a,3i5)')'N123 = ',n1,n2,n3
              write(*,'(a,3i5)')'NEC  = ',nec1,nec2,nec3
              stop'stop, id123 for parent'
           end if
500     continue

        ipmax = 0
        do i=1,ict_mat
           g1 = pool(i,1)
           n1 = pool(i,2)
           if(pool(i,3).gt.ipmax)then
             ipmax = pool(i,3)
             mp    = i
           end if
        end do
        if(id_proc.eq.0.and.ndebug.eq.1)then
          g1 = pool(mp,1)
          n1 = pool(mp,2)
          write(filenum,'(11x,a,2i3,a,i3,a,i5)') &
          '>> Design ',g1,n1,' Became Parent',ipmax,' times',ict_mat
        end if

        return
        if(id_proc.eq.0)then
        write(104,'(9999f15.8)')(d(g1,n1,ndv+k),k=1,nobj),FM*1.5d0
        do i=1,npop/2
        do j=1,2
           write(104,'(9999f15.8)') &
             (parent(i,j,ndv+k),k=1,nobj),parent(i,j,ndv+nobj+1)
        end do
        end do
        close(104)
        end if
        end subroutine mating
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine SUS(it,idiv,ict,ndim,prob,id)
        implicit none
        integer, intent(inout) :: it
        integer, intent(in)    :: idiv,ndim,ict
        double precision, dimension(0:ndim) :: prob
        integer, intent(out) :: id
        integer :: i
        double precision :: ran,step,star,targ

        id = 0
        if(idiv.le.0)stop'idiv in SUS'
        if(it.le.0.or.it.gt.idiv)stop'itarget in SUS'
        step = 1.d0 / dble(idiv)
        star = dble(it-1) / dble(idiv)

        call random_number(ran)
        targ = star + ran*step

        do i=1,ict
          if(prob(i-1).le.targ.and.prob(i).ge.targ)then
            id = i
            go to 100
          end if
        end do
        stop'no selection in SUS'

100     continue
        it = it + 1
        if(it.gt.idiv)it = 1
        
        end subroutine SUS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine decide_power(FM,RM,power)
        use dimGA
        implicit none
        double precision, intent(in) :: FM,RM
        double precision, intent(out) :: power
        integer :: iter,g1,n1,i
        double precision :: poi,eps,f,df

        iter = 0
        eps  = 1.d-5
        if(Cprob(:3).eq.'old'.or.Cprob(:3).eq.'Old')then
           poi  = 1.d0
        else if(Cprob(:3).eq.'new'.or.Cprob(:3).eq.'New')then
           poi  = 10.d0
        end if

100     continue
        iter = iter + 1

        if(Cprob(:3).eq.'old'.or.Cprob(:3).eq.'Old')then
           f =  RM**poi
          df = (RM**poi)*log(RM)
        else if(Cprob(:3).eq.'new'.or.Cprob(:3).eq.'New')then
           f =  FM**poi
          df = (FM**poi)*log(FM)
        end if
        do i=1,ipool
           g1 = pool(i,1)
           n1 = pool(i,2)
           if(Cprob(:3).eq.'old'.or.Cprob(:3).eq.'Old')then
             f =  f - fac_prob*  (Rrank(g1,n1)**poi)
            df = df - fac_prob*( (Rrank(g1,n1)**poi)*log(Rrank(g1,n1)) )
           else if(Cprob(:3).eq.'new'.or.Cprob(:3).eq.'New')then
             f =  f - fac_prob*  (Frank(g1,n1)**poi)
            df = df - fac_prob*( (Frank(g1,n1)**poi)*log(Frank(g1,n1)) )
           end if
        end do
!       write(*,'(i5,3e15.5)')iter,f,df,poi
        if(dabs(f).lt.eps)go to 120
        if(iter.ge.100)go to 110
        if(df.eq.0.d0)go to 110   !stop'df=0 in pool'
        poi = poi - f/df
        go to 100
110     continue
        poi = 1.d0
120     continue
        power = poi

        if(power.ge.0.01d0.and.power.le.10.d0)then
        else
           power = 1.d0
        end if
        end subroutine decide_power
        subroutine pregnancy(gen)
        use dimGA
        implicit none
        integer, intent(in) :: gen
        integer :: i,k,ict
        double precision :: Fp,Fm,Ff
        double precision, dimension(ndv)   :: papa,mama,flirt
        double precision, dimension(ndv,2) :: baby
! we have data until gen-1 generation

        if(mod(gen,ngen_grad).eq.0)then
          call gradient_evolution(gen) ! include evaluation of children
          return
        end if

        ict = 0
        do 100 i=1,npop/2
           do 110 k=1,ndv
             papa( k) = parent(i,1,k)
             mama( k) = parent(i,2,k)
             flirt(k) = parent(i,3,k)
110        continue
           if(flirt(1).lt.0.d0)then
             if(Cross(:3).eq.'BLX'.or.Cross(:3).eq.'blx')then
               call crossover_blx(ndv,if_crs,beta_cr,papa,mama,baby)
             else if(Cross(:3).eq.'SBX'.or.Cross(:3).eq.'sbx')then
               call crossover_sbx(ndv,if_crs,beta_cr,papa,mama,baby)
             else
               stop'unknown crossover'
             end if
           else ! directional crossover
             Fp = parent(i,1,ndv+nobj+1)
             Fm = parent(i,2,ndv+nobj+1)
             Ff = parent(i,3,ndv+nobj+1)
             call crossover_drx(ndv,papa,mama,flirt,Fp,Fm,Ff,baby)
           end if
           call boundbaby(1,ndv,baby)
           call mutation(Cmut,p_mut,r_mut,eta_mut,ndv,baby)
           call boundbaby(2,ndv,baby)

           ict = ict + 1
           do k=1,ndv
              d(gen,ict,k) = baby(k,1)
           end do
           ict = ict + 1
           do k=1,ndv
              d(gen,ict,k) = baby(k,2)
           end do
100     continue
        if(ict.ne.npop)stop'ict.ne.npop in pregnancy'

        end subroutine pregnancy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine boundbaby(mode,ndv,baby)
        implicit none
        integer, intent(in) :: mode,ndv
        double precision, dimension(ndv,2), intent(inout) :: baby
        integer :: i,j
        if(mode.ne.1.and.mode.ne.2)stop'mode of boundbaby'

        do 100 i=1,2
        do 110 j=1,ndv
           if(baby(j,i).lt.0.d0)then
             if(mode.eq.1)then
                baby(j,i) =  0.d0
             else if(mode.eq.2)then
                baby(j,i) = -1.d0 * baby(j,i)
                if(baby(j,i).gt.1.d0)baby(j,i) = 0.d0
             end if
           else if(baby(j,i).gt.1.d0)then
             if(mode.eq.1)then
                baby(j,i) =  1.d0
             else if(mode.eq.2)then
                baby(j,i) =  2.d0 - baby(j,i)
                if(baby(j,i).lt.0.d0)baby(j,i) = 1.d0
             end if
           end if
110     continue
100     continue

        end subroutine boundbaby
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine mutation(Ctype,p,r,eta,ndv,baby)
        implicit none
        character(len=30), intent(in) :: Ctype
        double precision, intent(in)  :: p,r,eta
        integer, intent(in) :: ndv
        double precision, dimension(ndv,2), intent(inout) :: baby
        integer :: i,j
        double precision :: ran,dxm
! eta is active only for new mutation
! -1<eta<infinity, correspond to 1<delta_ave<0.02
! when eta=0, distribution of delta is just a line
! when eta increase, delta_ave will reduce
! eta       : -1,  0,0.5,1.4,3.3,5.2,9.7,11,14,17,21,29,44,89,inf
! delta_ave :  1,0.5,0.4,0.3,0.2,.15,0.1,09,08,07,06,05,04,03,0.02
        if(p.gt.1.d0)stop'unknown prob of mutation'
        if(Ctype(:3).eq.'New'.or.Ctype(:3).eq.'new')then
          if(eta.lt.0.d0)stop'unreasonable eta of mutation'
        end if

        do 100 i=1,2
        do 110 j=1,ndv
           call random_number(ran)
           if(ran.ge.0.5d0*p.and.ran.le.(1.d0-0.5d0*p))go to 110

           if(Ctype(:3).eq.'Old'.or.Ctype(:3).eq.'old')then
             call random_number(ran)
             dxm = r * 2.d0 * (ran-0.5d0)  ! [-r:r]
           else if(Ctype(:3).eq.'New'.or.Ctype(:3).eq.'new')then
             call random_number(ran)
             if(ran.le.0.5d0)then
               dxm = (2.d0*ran)**(1.d0/(eta+1))-1.d0
             else
               dxm = 1.d0 - (2.d0*(1.d0-ran))**(1.d0/(eta+1)) ![-1:1]
             end if
             if(dxm.gt.r)then
                dxm = r
             else if(dxm.lt.-1.d0*r)then
                dxm = -1.d0 * r            ! [-r:r]
             end if
           else
             stop'unknown type of mutation'
           end if
           baby(j,i) = baby(j,i) + dxm
110     continue
100     continue
     
        end subroutine mutation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine crossover_sbx(ndv,iflg,eta,papa,mama,baby)
        implicit none
        integer, intent(in) :: ndv,iflg
        double precision, intent(in) :: eta
        double precision, dimension(ndv), intent(in)    :: papa,mama
        double precision, dimension(ndv,2), intent(out) :: baby
        integer :: i
        double precision :: ran,bet
! when eta increases, baby will be closer to parents
! ran<0.5 : interpolation, ran>0.5 : extrapolation
! when eta=5, the range is almost same with blx-0.5 (not probability)
!      eta=[1;5] is a compromised set?
! eta:         -1   0 .075 .3 .8    1    2    3    4    5   10  inf
! bet_in_ave:   0 0.5            0.66 0.74 0.79 0.82 0.85 0.90  1.0
! bet_ex_ave: inf   5    4  3  2 1.80 1.48 1.33 1.26 1.21 1.12  1.0
! when u(ran)=0.99
!       eta : 5.0 2.5 1.8 1.4 1.0 0.7
!       bet : 2   3   4   5   7   10
        if(eta.le.-0.5)stop'eta in SBX'

        if(iflg.eq.1)call random_number(ran)
        do i=1,ndv
           if(iflg.ne.1)call random_number(ran)
           if(ran.le.0.5d0)then
             bet = (2.d0*ran)**(1.d0/(eta+1.d0))
           else
             if(ran.gt.0.99d0)ran=0.99d0
             bet = (1.d0/(2.d0*(1.d0-ran)))**(1.d0/(eta+1.d0))
           end if
           baby(i,1) = 0.5d0*( (1.d0+bet)*papa(i)+(1.d0-bet)*mama(i) )
           baby(i,2) = 0.5d0*( (1.d0-bet)*papa(i)+(1.d0+bet)*mama(i) )
        end do

        end subroutine crossover_sbx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine crossover_blx(ndv,iflg,blx,papa,mama,baby)
        implicit none
        integer, intent(in) :: ndv,iflg
        double precision, intent(in) :: blx
        double precision, dimension(ndv), intent(in)    :: papa,mama
        double precision, dimension(ndv,2), intent(out) :: baby
        integer :: i
        double precision :: ran,gam
! blx = 0.5 is usually recommended (same probability of inter/extrapolation)
! blx = 0.0 is perfectly limited within parents

        if(iflg.eq.1) call random_number(ran)
        do i=1,ndv
           if(iflg.ne.1) call random_number(ran)
           gam = (1.d0 + 2.d0 * blx) * ran - blx
           baby(i,1) = gam*papa(i) + (1.d0-gam)*mama(i)
           baby(i,2) = gam*mama(i) + (1.d0-gam)*papa(i)
        end do
 
        end subroutine crossover_blx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine crossover_drx(ndv,papa,mama,flirt,Fp,Fm,Ff,baby)
        implicit none
        integer, intent(in) :: ndv
        double precision, dimension(ndv), intent(in)    :: papa,mama,flirt
        double precision, intent(in)                    :: Fp,Fm,Ff
        double precision, dimension(ndv,2), intent(out) :: baby

        integer :: i
        double precision :: ran,ran1,ran2,d1,d2
        double precision, dimension(ndv)  :: x0,x1,x2
        double precision                  :: f0,f1,f2

        call random_number(ran)
        if(ran.le.1./3.)then
           x0 = papa
           x1 = mama
           x2 = flirt
           f0 = Fp
           f1 = Fm
           f2 = Ff
        else if(ran.gt.1./3..and.ran.le.2./3.)then
           x0 = mama
           x1 = papa
           x2 = flirt
           f0 = Fm
           f1 = Fp
           f2 = Ff
        else
           x0 = flirt
           x1 = mama
           x2 = papa
           f0 = Ff
           f1 = Fp
           f2 = Fm
        end if
! F should be maximized

        do 100 i=1,2
          call random_number(ran1)
          call random_number(ran2)
          d1 = sign(1.d0,(f0-f1))
          d2 = sign(1.d0,(f0-f2))
          baby(:,i) = x0(:) &
                    + ran1 * d1 * (x0(:)-x1(:)) &
                    + ran2 * d2 * (x0(:)-x2(:))
100     continue

        end subroutine crossover_drx
        subroutine ranking(gen)
        use dimGA
        implicit none
        integer, intent(in) :: gen
        integer :: n1,n2,g1,g2,irank,np
! we have data until gen-1 generation

! within the generation
        rank(gen-1,:) = 1.d0
        do 100 n1 = 1,npop-1
          do 110 n2 = n1+1,npop
             call Pareto_Ranking(gen-1,n1,gen-1,n2)
110       continue
100     continue

! with previous generations
        g1 = gen-1
        do 200 n1=1,npop
          do 210 g2 = 0,gen-2
            do 220 n2=1,npop
               call Pareto_Ranking(g1,n1,g2,n2)
220         continue
210       continue
200     continue

! check
        np = 0
        do 250 g1=0,gen-1
          do 260 n1=1,npop
             irank = int( rank(g1,n1) )
             if(irank.eq.1)then
               np = np + 1
             else if(irank.lt.1)then
               stop'rank < 1'
             end if
             Rrank(g1,n1) = 1.d0 / rank(g1,n1)
260       continue
250     continue
!       write(*,*)'np =',np

        end subroutine ranking
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine Pareto_Ranking(g1,n1,g2,n2)
        use dimGA
        implicit none
        integer, intent(in) :: g1,g2,n1,n2
        integer :: isum1,isum2,itie,k
        double precision :: con1,con2,f1,f2

        con1 = d(g1,n1,ndv+nobj+1)
        con2 = d(g2,n2,ndv+nobj+1)
        if(     con1.eq.0.d0.and.con2.eq.0.d0)then  ! n1/n2 are feasible
          isum1 = 0
          isum2 = 0
          itie  = 0
          do k=1,nobj
             f1 = d(g1,n1,ndv+k)
             f2 = d(g2,n2,ndv+k)
             if(Cobj(k).eq.'Max'.or.Cobj(k).eq.'max')then
               if(f1.lt.f2)then
                 isum1 = isum1 + 1
               else if(f1.gt.f2)then
                 isum2 = isum2 + 1
               else if(f1.eq.f2)then
                 itie  = itie  + 1
               end if
             else if(Cobj(k).eq.'Min'.or.Cobj(k).eq.'min')then
               if(f1.lt.f2)then
                 isum2 = isum2 + 1
               else if(f1.gt.f2)then
                 isum1 = isum1 + 1
               else if(f1.eq.f2)then
                 itie  = itie  + 1
               end if
             else
               stop'unknown set of Cobj in ParetoRank'
             end if
          end do
          if(isum1.eq.nobj-itie.and.itie.ne.nobj)then
            if(isum2.ne.0)stop'isum in ParetoRank'
            rank(g1,n1) = rank(g1,n1) + 1.d0
          else if(isum2.eq.nobj-itie.and.itie.ne.nobj)then
            if(isum1.ne.0)stop'isum in ParetoRank'
            rank(g2,n2) = rank(g2,n2) + 1.d0
          end if

        else if(con1.gt.0.d0.and.con2.eq.0.d0)then  ! n1 is unfeasible
          rank(g1,n1) = rank(g1,n1) + 1.d0

        else if(con1.eq.0.d0.and.con2.gt.0.d0)then  ! n2 is unfeasible
          rank(g2,n2) = rank(g2,n2) + 1.d0

        else if(con1.gt.0.d0.and.con2.gt.0.d0)then  ! n1/n2 are unfeasible
          if(con1.gt.con2)then
            rank(g1,n1) = rank(g1,n1) + 1.d0
          else if(con1.lt.con2)then
            rank(g2,n2) = rank(g2,n2) + 1.d0
          end if
        else
          stop'unknown constraint set in Pareto Ranking'
        end if
        
        end subroutine Pareto_Ranking
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine read_param(Cfile)
        use dimGA
        use dimKrig, only:filenum
        implicit none
        character(len=30), intent(in) :: Cfile
        character(len=30) :: Cdum
        integer :: i,leng,iflag,istat,nseed
        double precision :: ran

        integer :: mdv,mobj,mdat,meva,mdebug
        character(len=3) :: Cmobj(10)

        open(10,file=Cfile, form='formatted',status='old')
! basic set
        read(10,*)Cdum
        read(10,*)npop
        read(10,*)ngen
        read(10,*)mdv  !ndv
        read(10,*)mobj !nobj
        backspace 10
        read(10,*)mobj,(Cmobj(i),i=1,mobj)
!       read(10,*)nobj,(Cobj(i),i=1,nobj)
        read(10,*)mdat !ndat
        read(10,*)meva !neva
        read(10,*)nseed
        do i=1,nobj
           if(Cobj(i).ne.'Min'.and.Cobj(i).ne.'min'.and. &
              Cobj(i).ne.'Max'.and.Cobj(i).ne.'max')     &
           stop'Unknown Objective Set Min/Max'
        end do
        if(mod(npop,2).ne.0)stop'npop should be even'
        if(nobj.gt.10)stop'too many objectives'
        ndat = max( ndat, ndv+nobj+1 )
        ngen_fin = ngen
        do i=1,nseed
           call random_number(ran)
        end do

        if(id_proc.eq.0)then
          if(mdv.ne.ndv)then
            write(filenum,*)'*mdv   != ndv in ',Cfile
          end if
          if(mobj.ne.nobj)then
            write(filenum,*)'*mobj  != nobj in ',Cfile
          end if
          do i=1,nobj
          if(Cmobj(i).ne.Cobj(i))then
            write(filenum,*)'*Cmobj != Cobj in ',Cfile
          end if
          end do
          if(mdat.ne.ndat)then
            write(filenum,*)'*mdat  != ndat in ',Cfile
          end if
          if(meva.ne.neva)then
            write(filenum,*)'*meva  != neva in ',Cfile
          end if
        end if
          
! initial pop
        read(10,'(a3)')Cini
        if(Cini.eq.'Use'.or.Cini.eq.'use')then
          backspace 10
          read(10,'(a,a)')Cini,Cinifile
        end if

! sharing
        read(10,*)Cdum
        read(10,*)Ctype_sh
        read(10,*)ntype_sh
        read(10,*)ngen_sh
        read(10,*)alpha_sh
        read(10,*)Cpareto
        if(Ctype_sh(:12).ne.'Conventional'.and. &
           Ctype_sh(:12).ne.'conventional'.and. &
           Ctype_sh( :3).ne.'New'.and. &
           Ctype_sh( :3).ne.'new') &
           stop'Type of Sharing, Conventional or New'

! mating and selection
        read(10,*)Cdum
        read(10,*)Cprob
        read(10,*)fac_prob
        read(10,*)ngen_mat
        read(10,*)id_sus
        if(ngen_mat.gt.ngen_sh.and.ngen_sh.ne.0) &
        stop'ngen_mat.gt.ngen_sh'
        if(id_sus.le.0)stop'id_sus<=0'
        if(Cprob(:3).ne.'New'.and.Cprob(:3).ne.'new'.and. &
           Cprob(:3).ne.'Old'.and.Cprob(:3).ne.'old')     &
           stop'Type of Fitness Function, New or Old'

! crossover
        read(10,*)Cdum
        read(10,*)Cross
        read(10,*)if_crs
        read(10,*)beta_cr
        read(10,*)p_drxs,p_drxe
        read(10,*)p_necs,p_nece
        if(Cross(:3).ne.'BLX'.and.Cross(:3).ne.'blx'.and. &
           Cross(:3).ne.'SBX'.and.Cross(:3).ne.'sbx')     &
           stop'Type of Crossover, BLX or SBX'
        if(p_drxs.lt.0.d0.or.p_drxs.gt.1.d0)stop'unknown p_drxs'
        if(p_drxe.lt.0.d0.or.p_drxe.gt.1.d0)stop'unknown p_drxe'
        if(p_necs.lt.0.d0.or.p_necs.gt.1.d0)stop'unknown p_necs'
        if(p_nece.lt.0.d0.or.p_nece.gt.1.d0)stop'unknown p_nece'

! gradient
        read(10,*)Cdum
        read(10,*)ngen_grad
        if(ngen_grad.le.0)ngen_grad = ngen*2
! mutation
        read(10,*)Cdum
        read(10,*)Cmut
        read(10,*)p_muts,p_mute
        read(10,*)r_mut,eta_mut
        if(Cmut(:3).ne.'New'.and.Cmut(:3).ne.'new'.and. &
           Cmut(:3).ne.'Old'.and.Cmut(:3).ne.'old')     &
           stop'Type of Mutation, New or Old'
        if(p_muts.lt.0.d0.or.p_muts.gt.1.d0)stop'unknown prob_mutation'
        if(p_mute.lt.0.d0.or.p_mute.gt.1.d0)stop'unknown prob_mutation'
        if(r_mut.lt.0.d0.or.r_mut.gt.1.d0)stop'unknown range_mutation'
        if(eta_mut.lt.0.d0)stop'unreasonable eta of mutation'

! convergence criteria / output
        read(10,*)Cdum
        read(10,*)ngen_stop
        read(10,*)neva_max
        read(10,*)ngen_out
        read(10,*)mdebug !ndebug
!       if(fac_mut.lt.1.d0)stop'unreasonable factor of mutation'
!       if(ngen_mut.gt.ngen_stop)stop'unreasonable set of ngen_mut/stop'
        close(10)

        leng = index(Cfile,'  ') - 1
        if(id_proc.eq.0)then
        write(filenum,'(a,a)')' >> Reading Parameter from ',Cfile(:leng)
        write(filenum,'(6x,a,i5)')'>> # of Population       = ',npop
        write(filenum,'(6x,a,i5)')'>> # of Generation       = ',ngen
        write(filenum,'(6x,a,i5)')'>> # of Design Variables = ',ndv
        write(filenum,'(6x,a,i5)')'>> # of Objectives       = ',nobj
        write(filenum,'(6x,a,10(a,x))')'>> Set of Objectives     :   ', &
        (Cobj(i),i=1,nobj)
        end if
        call check_neva(neva)
        call check_param

        iflag = 0
        allocate( d(0:ngen,npop,ndat), stat=istat )
                                         if(istat.ne.0)iflag = iflag + 1
        allocate( rank(0:ngen,npop), stat=istat )
                                         if(istat.ne.0)iflag = iflag + 1
        allocate( Frank(0:ngen,npop), stat=istat )
                                         if(istat.ne.0)iflag = iflag + 1
        allocate( Rrank(0:ngen,npop), stat=istat )
                                         if(istat.ne.0)iflag = iflag + 1
        allocate( parent(npop/2,3,ndat), stat=istat )
                                         if(istat.ne.0)iflag = iflag + 1
        allocate( pool((ngen+1)*npop,3), stat=istat )
                                         if(istat.ne.0)iflag = iflag + 1
        if(iflag.ne.0)stop'allocation in read_param'
        d        = 0.d0
        rank     = 0.d0
        Frank    = 0.d0
        Rrank    = 0.d0
        parent   = 0.d0
        pool     = 0
        ipool    = 0
        ict_eva  = 0

        iaxs_nec = 0
        p_mut    = p_muts
        p_drx    = p_drxs
        p_nec    = p_necs
 
        end subroutine read_param


        subroutine sharing(gen)
        use dimGA
        use dimKrig, only:filenum
        implicit none
        include 'mpif.h'
        integer, intent(in) :: gen
        double precision, dimension(0:ngen,npop) :: shnc
        integer, dimension((ngen+1)*npop,2) :: id_sh
        integer, dimension((ngen+1)*npop  ) :: rcount
        integer :: i,j,k,ir,igen,ict_sh
        integer :: i1st,npareto,iter
        integer :: n1,g1,n2,g2,irank,rmax
        integer :: n,isum,ierr
        double precision :: radius,radini,ri,f,df,eps
        double precision :: sh,shmin,shmax,f1,f2,dd
! we have data until gen-1 generation
! even for monoobjective, execute to get elite, fitness function

! get Fmin,Fmax,elite
        rmax    = 0
        i1st    = 0
        npareto = 0
        do 100 g1=0,gen-1
        do 100 n1=1,npop
          irank = int( rank(g1,n1) )
          rmax = max(rmax,irank)
          if(irank.ne.1)go to 100

          npareto = npareto + 1
          if(i1st.eq.0)then
            do i=1,nobj
              Fmin(i)    = d(g1,n1,ndv+i)
              Fmax(i)    = d(g1,n1,ndv+i)
              elite(i,1) = g1
              elite(i,2) = n1
            end do
            i1st = 1
          else
            do i=1,nobj
               if(d(g1,n1,ndv+i).lt.Fmin(i))then
                  Fmin(i)    = d(g1,n1,ndv+i)
                  if(Cobj(i).eq.'Min'.or.Cobj(i).eq.'min')then
                    elite(i,1) = g1
                    elite(i,2) = n1
                  end if
               end if
               if(d(g1,n1,ndv+i).gt.Fmax(i))then
                  Fmax(i)    = d(g1,n1,ndv+i)
                  if(Cobj(i).eq.'Max'.or.Cobj(i).eq.'max')then
                    elite(i,1) = g1
                    elite(i,2) = n1
                  end if
               end if
            end do
          end if
100     continue
        do i=1,nobj
           if(elite(i,1).lt.0)stop'elite-1'
           if(elite(i,2).le.0)stop'elite-2'
        end do
        if(id_proc.eq.0.and.mod(gen-1,ngen_out).eq.0)then
          write(filenum,'(6x,a,i4,a,i3,a,99e11.3)') &
          '>> ',npareto,' Pareto Designs until ', gen-1,' Gen.,', &
          (d(elite(i,1),elite(i,2),ndv+i),i=1,nobj)
        end if

! Sharing Population
        if(ngen_sh.eq.0)then
           igen = 0
        else if(ngen_sh.gt.0)then
           igen = (gen-1) - (ngen_sh-1)
           if(igen.lt.0)igen = 0
        else
           stop'negative ngen_sh'
        end if
        id_sh  = 0
        ict_sh = 0
        do 200 i=1,nobj          ! elite is the first
           g1 = elite(i,1)
           n1 = elite(i,2)
           if(g1.lt.igen)then
             ict_sh = ict_sh + 1
             id_sh(ict_sh,1) = g1
             id_sh(ict_sh,2) = n1
           end if
200     continue
        if(mod(gen-1,ngen_grad).eq.0.and.gen.ne.1)then
          ! previous generation is given by gradient-based evolution
          do 210 g1=gen-1,igen,-1  ! children, parent, and more
          do 210 n1=1,npop
             ict_sh = ict_sh + 1
             id_sh(ict_sh,1) = g1
             id_sh(ict_sh,2) = n1
210       continue
        else
          if(ngen_sh.ne.1.and.gen.ge.2)then
          do 220 g1=gen-2,gen-2    ! parent
          do 220 n1=1,npop
             ict_sh = ict_sh + 1
             id_sh(ict_sh,1) = g1
             id_sh(ict_sh,2) = n1
220       continue
          end if
          do 230 g1=gen-1,gen-1    ! children
          do 230 n1=1,npop
             ict_sh = ict_sh + 1
             id_sh(ict_sh,1) = g1
             id_sh(ict_sh,2) = n1
230       continue
          do 240 g1=gen-3,igen,-1  ! grandparent and more
          do 240 n1=1,npop
             ict_sh = ict_sh + 1
             id_sh(ict_sh,1) = g1
             id_sh(ict_sh,2) = n1
240       continue
        end if

        if(Cpareto(:7).eq.'Include'.or.Cpareto(:7).eq.'include')then
          do 300 g1=igen-1,0,-1
          do 300 n1=1,npop
             irank = int(rank(g1,n1))
             if(irank.ne.1)go to 300
             do 310 i=1,nobj
               g2 = elite(i,1)
               n2 = elite(i,2)
               if(g1.eq.g2.and.n1.eq.n2)go to 300
310          continue
             ict_sh = ict_sh + 1
             id_sh(ict_sh,1) = g1
             id_sh(ict_sh,2) = n1
300       continue
        end if

! Mono-objective Fitness
        rcount = 0
        rmax   = 0
        do 400 i=1,ict_sh
           g1    = id_sh(i,1)
           n1    = id_sh(i,2)
           irank = int(rank(g1,n1))
           rmax = max(rmax,irank)
           if(irank.lt.1)then
             if(id_proc.eq.0)&
             write(*,'(a,6i6)')'*rank<1 in share',id_proc,i,g1,n1,irank,gen
             call stop_all
           end if
           if(rmax.gt.(ngen+1)*npop)stop'rmax>(ngen+1)*npop'
           rcount(irank) = rcount(irank) + 1
400     continue
        Frank = -1.d0
        do 450 i=1,ict_sh
           g1    = id_sh(i,1)
           n1    = id_sh(i,2)
           irank = int(rank(g1,n1))
           isum  = 0
           do 460 j=1,irank-1
              isum = isum + rcount(j)
460        continue
           Frank(g1,n1) = dble(ict_sh)        &
                        - 0.5d0*(dble(rcount(irank))-1.d0) &
                        - dble(isum)
450     continue

        if(nobj.eq.1)return ! I only wanted to get elite,fitness
!       if(gen.eq.ngen+1)return ! Post-Process to get elite,fitness


! initial guess of r
!       radius = 0.d0
!       do i=1,nobj
!          radius = radius + (Fmax(i)-Fmin(i))**2
!       end do
        if(ntype_sh.eq.0)then
          N    = npareto
        else if(ntype_sh.lt.0)then
          if(ntype_sh.eq.-1)stop'ntype_sh'
          N    = min(npareto,ict_sh,-1*ntype_sh)
        else
          N    = ict_sh
        end if
        radius = dble(nobj)
        radius = sqrt(radius) / dble(N)
        radius = (radius) ** (1.d0/dble(nobj-1))
        radini = radius
        if(nobj.eq.2)then
          radini = 2.d0 / dble(N-1)
        end if

! Newton Method for Radius
        iter = 0
        eps  = 1.d-8
499     continue
        ri   = radini
500     continue
        iter = iter + 1
        f    = ((1.d0+ri)**nobj) - 1.d0 - dble(N)*((ri)**nobj)
        call MPI_BCAST(ri,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(f,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(iter,1,MPI_INTEGER,      0,MPI_COMM_WORLD,ierr)
        if(dabs(f).lt.eps)then
          if(ri.gt.eps)then
             go to 520
          else
             if(radini.ge.100.d0)then
               write(*,'(a,2i5,f10.4)')'*no sharing radius?',N,nobj,radius
               go to 510
             end if
             radini = radini * 2.d0
             go to 499
          end if
        end if
        if(iter.ge.1000)go to 510
        df  = dble(nobj)  *((1.d0+ri)**(nobj-1)) &
            - dble(N*nobj)*((     ri)**(nobj-1))
        if(df.eq.0.d0)stop'df=0.d0 in sharing'
!       if(id_proc.eq.0)write(103,'(i5,3e15.5)')iter,ri,f,df
        ri  = ri - f / df
        go to 500
510     continue
        ri = radius
520     continue
!       close(103)
        if(id_proc.eq.0.and.ndebug.eq.1) &
        write(filenum,'(11x,a,2f10.4)')'>> Sharing Radius      = ',radius,ri
        radius = ri
        if(radius.le.0.d0)stop'sharing radius <= 0.d0'

! Initial Set of Niche Count
        if(Ctype_sh(:12).eq.'Conventional'.or. &
           Ctype_sh(:12).eq.'conventional')then
           shnc = 0.d0
           go to 599
        else if(Ctype_sh( :3).eq.'New'.or. &
                Ctype_sh( :3).eq.'new')then
           shnc = 1.d0
           go to 699
        else
           stop'unknown sharing function'
        end if

! Conventional Sharing Function, Consider me and me
599     continue
        do 600 ir=1,rmax
          do 610 i=1,ict_sh
            g1 = id_sh(i,1)
            n1 = id_sh(i,2)
            irank = int( rank(g1,n1) )
            if(irank.ne.ir)go to 610
            do 620 j=1,ict_sh
              g2 = id_sh(j,1)
              n2 = id_sh(j,2)
              irank = int( rank(g2,n2) )
              if(irank.ne.ir)go to 620
              dd = 0.d0
              do 630 k=1,nobj
                 f1 = d(g1,n1,ndv+k)
                 f2 = d(g2,n2,ndv+k)
                 if(Fmax(k)-Fmin(k).gt.0.d0)then
                 dd = dd + ( (f1-f2)/(Fmax(k)-Fmin(k)) )**2
                 else
                 dd = dd + ( (f1-f2)/(    0.0001d0   ) )**2
                 end if
630           continue
              dd = sqrt(dd)
              if(dd.ge.radius)go to 620
              sh = 1.d0 - (dd/radius)**alpha_sh
              shnc(g1,n1) = shnc(g1,n1) + sh
620         continue
610       continue
600     continue
        go to 899

! New Sharing Function (AIAA-2009-1168)
699     continue
        do 700 ir=1,rmax
          do 710 i=1,ict_sh-1                ! note the difference
            g1 = id_sh(i,1)
            n1 = id_sh(i,2)
            irank = int( rank(g1,n1) )
            if(irank.ne.ir)go to 710
            do 720 j=i+1,ict_sh              ! note the difference
              g2 = id_sh(j,1)
              n2 = id_sh(j,2)
              irank = int( rank(g2,n2) )
              if(irank.ne.ir)go to 720
              dd = 0.d0
              do 730 k=1,nobj
                 f1 = d(g1,n1,ndv+k)
                 f2 = d(g2,n2,ndv+k)
                 if(Fmax(k)-Fmin(k).gt.0.d0)then
                 dd = dd + ( (f1-f2)/(Fmax(k)-Fmin(k)) )**2
                 else
                 dd = dd + ( (f1-f2)/(    0.0001d0   ) )**2
                 end if
730           continue
              dd = sqrt(dd)
              if(dd.ge.radius)go to 720
              sh = min( (dd/radius)**alpha_sh, 1.d0 )
              if(shnc(g1,n1).lt.shnc(g2,n2))then
                 shnc(g1,n1) = shnc(g1,n1) * sh
              else
                 shnc(g2,n2) = shnc(g2,n2) * sh
              end if
720         continue
710       continue
700     continue
        go to 899

! Final Check of Sharing Function
899     continue
        shmin = +10000.d0
        shmax = -10000.d0
        do 900 i = 1,ict_sh
           g1 = id_sh(i,1)
           n1 = id_sh(i,2)
           shmin = min(shnc(g1,n1),shmin)
           shmax = max(shnc(g1,n1),shmax)
           if(Ctype_sh(:12).eq.'Conventional'.or. &
              Ctype_sh(:12).eq.'conventional')then
              if(shnc(g1,n1).lt.1.d0)stop'niche_old is less than 1'
              Frank(g1,n1) = Frank(g1,n1) / shnc(g1,n1)
              Rrank(g1,n1) = Rrank(g1,n1) / shnc(g1,n1) 
           else if(Ctype_sh( :3).eq.'New'.or. &
                   Ctype_sh( :3).eq.'new')then
              if(shnc(g1,n1).lt.0.d0)stop'niche_new is less than 0'
              if(shnc(g1,n1).gt.1.d0)stop'niche_new is greater than 1'
              Frank(g1,n1) = Frank(g1,n1) * shnc(g1,n1)
              Rrank(g1,n1) = Rrank(g1,n1) * shnc(g1,n1)
           else
              stop'ctype_sh'
           end if
900     continue
        if(id_proc.eq.0.and.ndebug.eq.1)then
           write(filenum,'(11x,a,2f10.4)')'>> Max/Min Niche Count = ',shmax,shmin
        end if

        return
        if(id_proc.eq.0)then
          do g1=0,gen-1
          do n1=1,npop
             write(103,'(9999f15.8)') &
             (d(g1,n1,ndv+i),i=1,nobj),rank(g1,n1),Frank(g1,n1),Rrank(g1,n1)
          end do
          end do
        end if
        end subroutine sharing
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine check_neva(i)
        use dimGA
        use dimKrig, only:filenum
        implicit none
        integer, intent(in) :: i
        if(id_proc.eq.0)then
        if(     i.eq.-1)then
          write(*,'(6x,a)')'>> 2MO Test Function'
        else if(i.eq.-2)then
          write(*,'(6x,a)')'>> 2MO 2DV Test Function'
        else if(i.eq.-3)then
          write(*,'(6x,a)')'>> 2MO 3DV Test Function'
        else if(i.eq.-4)then
          write(*,'(6x,a)')'>> 2MO 2DV Test Function with Constraint'
        else if(i.eq.-5)then
          write(*,'(6x,a)')'>> Rastrigin Function'
        else if(i.eq.-6)then
          write(*,'(6x,a)')'>> Rosenbrock Function'
        else if(i.eq.-7)then
          write(*,'(6x,a)')'>> 2MO Test to Max/Min'
        else if(i.eq.-10)then
          write(*,'(6x,a)')'>> 2MO KUR'
        else if(i.eq.-11)then
          write(*,'(6x,a)')'>> 2MO ZDT1 (Convex)'
        else if(i.eq.-12)then
          write(*,'(6x,a)')'>> 2MO ZDT2 (Conacave)'
        else if(i.eq.-13)then
          write(*,'(6x,a)')'>> 2MO ZDT3 (Discontinuous)'
        else if(i.eq.-14)then
          write(*,'(6x,a)')'>> 2MO ZDT4 (Local Optimals)'
        else if(i.eq.-16)then
          write(*,'(6x,a)')'>> 2MO ZDT6 (Non-Uniform)'
        else if(i.eq.1.or.i.eq.3)then
          write(filenum,'(6x,a)')'>> Likelihood Maximization for Theta'
        else if(i.eq.2)then
          write(filenum,'(6x,a)')'>> Likelihood Maximization for Theta/Power'
        else if(i.ge.4.and.i.le.6)then
          write(filenum,'(6x,a)')'>> Likelihood Maximization for VFM Theta'
        else if(i.ge.11.and.i.le.14)then
          write(filenum,'(6x,a)')'>> New Sample Points Search on Metamodel'
        else if(i.eq.15)then
          write(filenum,'(6x,a)')'>> New Sample Points Search by MaxVar'
        end if
        end if

        if( (i.lt.0).or.(i.ge.1.and.i.le.6).or.(i.ge.11.and.i.le.15) )then
          return
        else
          write(*,*)i
          stop'unknown neva in check'
        end if
        end subroutine check_neva
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine testfunction(xnew,ynew)
        use dimGA
        double precision, dimension(ndv),      intent(in)  :: xnew
        double precision, dimension(ndat-ndv), intent(out) :: ynew
        double precision, dimension(ndv) :: xscale
        integer :: i
        double precision :: sum1,sum2,ddv,cons1,cons2,pi,ff,gg,tmp

        pi   = 4.d0*datan(1.d0)

        if(neva.eq.-1)then
! Minimize, -4 < x < 4
          if(nobj.ne.2)stop'testfunc1, nobj'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'testfunc1, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'testfunc1, Cobj'
          do i=1,ndv
            xscale(i) =  8.d0 * xnew(i) - 4.d0
          end do
          sum1 = 0.d0
          sum2 = 0.d0
          do i=1,ndv
            sum1 = sum1 + (xscale(i)       )**2
            sum2 = sum2 + (xscale(i) - 2.d0)**2
          end do
          ynew(1) = sum1 / dble(ndv)
          ynew(2) = sum2 / dble(ndv)
          ynew(3) = 0.d0

        else if(neva.eq.-2)then
! Minimize, -4 < x < 4
          if(ndv.ne.2) stop'testfunc2, ndv'
          if(nobj.ne.2)stop'testfunc2, nobj'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'testfunc2, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'testfunc2, Cobj'
          do i=1,ndv
            xscale(i) =  8.d0 * xnew(i) - 4.d0
          end do
          ddv  = dble(ndv)
          sum1 = 0.d0
          sum2 = 0.d0
          do i=1,ndv
            sum1 = sum1 + ( xscale(i)-(1.d0/sqrt(ddv)) )**2
            sum2 = sum2 + ( xscale(i)+(1.d0/sqrt(ddv)) )**2
          end do
          ynew(1) = 1.d0 - exp(-1.d0 * sum1)
          ynew(2) = 1.d0 - exp(-1.d0 * sum2)
          ynew(3) = 0.d0

        else if(neva.eq.-3)then
! Minimize, ndv =  3, -5 < x < 5
          if(ndv.ne.3) stop'testfunc3, ndv'
          if(nobj.ne.2)stop'testfunc3, nobj'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'testfunc3, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'testfunc3, Cobj'
          do i=1,ndv
            xscale(i) = 10.d0 * xnew(i) - 5.d0
          end do
          sum1 = 0.d0
          sum2 = 0.d0
          do i=1,ndv-1
            sum1 = sum1 &
            - 10.d0 * exp( -0.2d0*sqrt(xscale(i)**2+xscale(i+1)**2) )
          end do
          do i=1,ndv
            sum2 = sum2 &
            + (abs( xscale(i) ))**0.8d0 + 5.d0*sin( xscale(i)**3 )
          end do
          ynew(1) = sum1
          ynew(2) = sum2
          ynew(3) = 0.d0
 
        else if(neva.eq.-4)then
! Minimize, ndv = 2, 0.1 < x1 < 1, 0 < x2 < 5, 2 constraints
          if(ndv.ne.2) stop'testfunc4, ndv'
          if(nobj.ne.2)stop'testfunc4, nobj'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'testfunc4, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'testfunc4, Cobj'
          xscale(1) =  0.9d0 * xnew(1) + 0.1d0
          xscale(2) =  5.0d0 * xnew(2) + 0.0d0
          ynew(1) = xscale(1)
          ynew(2) = (1.d0 + xscale(2)) / xscale(1)
          cons1   = 6.d0 - xscale(2) - 9.d0 * xscale(1)
          cons2   = 1.d0 + xscale(2) - 9.d0 * xscale(1)
          ynew(3) = 0.d0
          if(cons1.gt.0.d0)then
            ynew(3) = ynew(3) + 1.d0
!           ynew(3) = ynew(3) + cons1
          end if
          if(cons2.gt.0.d0)then
            ynew(3) = ynew(3) + 1.d0
!           ynew(3) = ynew(3) + cons2
          end if

        else if(neva.eq.-5)then
! Minimize Rastrigin, -5.12 < x < 5.12
          if(nobj.ne.1)stop'testfunc5, nobj'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'testfunc5, Cobj'
          do i=1,ndv
            xscale(i) = 10.24d0 * xnew(i) - 5.12d0
          end do
          sum1 = 0.d0
          do i=1,ndv
             sum1 = sum1 &
             + xscale(i)**2 - 10.d0 * cos( 2.d0 * pi * xscale(i) )
          end do
          ynew(1) = 10.d0*dble(ndv) + sum1
          ynew(2) =  0.d0

        else if(neva.eq.-6)then
! Minimize Rosenbrock
          if(nobj.ne.1)stop'testfunc6, nobj'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'testfunc6, Cobj'
          do i=1,ndv
            xscale(i) = 5.d0 * xnew(i) - 2.d0
          end do
          sum1 = 0.d0
          do i=1,ndv-1
            sum1 = sum1 &
            + (1.d0-xscale(i))**2 + 100.d0*(xscale(i+1)-xscale(i)**2)**2
          end do
          ynew(1) =  sum1
          ynew(2) =  0.d0

        else if(neva.eq.-7)then
! Test to Min/Max
          if(nobj.ne.2)stop'testfunc7, nobj'
          if(ndv.ne.2) stop'testfunc7, ndvj'

          ynew(1) =  xnew(1)
          ynew(2) = (xnew(1)**2+xnew(2)**2)*0.5d0

        else if(neva.eq.-10)then
! KUR
          if(nobj.ne.2)stop'KUR, nobj'
          if(ndv.lt.3)stop'KUR, ndv'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'KUR, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'KUR, Cobj'
          do i=1,ndv
            xscale(i) = -5.d0 + 10.d0*xnew(i)
          end do
          ff = 0.d0
          gg = 0.d0
          do i=1,ndv
             if(i.ne.ndv)then
               tmp = sqrt(xscale(i)**2+xscale(i+1)**2)
               ff = ff -10.d0*(exp(-0.2d0*tmp))
             end if
             gg = gg + (dabs(xscale(i)))**(0.8d0) &
                     + 5.d0*sin(xscale(i)**3)
          end do
          ynew(1) = ff
          ynew(2) = gg
          ynew(3) = 0.d0

        else if(neva.eq.-11)then
! ZDT1
          if(nobj.ne.2)stop'ZDT1, nobj'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'ZDT1, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'ZDT1, Cobj'
          xscale(:) = xnew(:)
          ff = xscale(1)
          gg = 1.d0
          do i=2,ndv
             gg = gg + (9.d0/dble(ndv-1))*xscale(i)
          end do
          if(ff.lt.0.d0)ff = 0.d0
          if(gg.le.0.d0)stop'ZDT1'
          ynew(1) = ff
          ynew(2) = gg*(1.d0-sqrt(ff/gg))
          ynew(3) = 0.d0

        else if(neva.eq.-12)then
! ZDT2
          if(nobj.ne.2)stop'ZDT2, nobj'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'ZDT2, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'ZDT2, Cobj'
          xscale(:) = xnew(:)
          ff = xscale(1)
          gg = 1.d0
          do i=2,ndv
             gg = gg + (9.d0/dble(ndv-1))*xscale(i)
          end do
          if(ff.lt.0.d0)ff = 0.d0
          if(gg.le.0.d0)stop'ZDT2'
          ynew(1) = ff
          ynew(2) = gg*(1.d0-(ff/gg)**2)
          ynew(3) = 0.d0

        else if(neva.eq.-13)then
! ZDT3
          if(nobj.ne.2)stop'ZDT3, nobj'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'ZDT3, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'ZDT3, Cobj'
          xscale(:) = xnew(:)
          ff = xscale(1)
          gg = 1.d0
          do i=2,ndv
             gg = gg + (9.d0/dble(ndv-1))*xscale(i)
          end do
          if(ff.lt.0.d0)ff = 0.d0
          if(gg.le.0.d0)stop'ZDT3'
          ynew(1) = ff
          ynew(2) = gg*(1.d0-sqrt(ff/gg)) &
                  - (ff/gg)*sin(10.d0*pi*(ff/gg))
          ynew(3) = 0.d0

        else if(neva.eq.-14)then
! ZDT4
          if(nobj.ne.2)stop'ZDT4, nobj'
          if(ndv.ne.10)stop'ZDT4, ndv'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'ZDT4, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'ZDT4, Cobj'
          xscale(1)   = xnew(1) !-1.d0 +  2.d0*xnew(1)
          do i=2,ndv
            xscale(i) = -5.d0 + 10.d0*xnew(i)
          end do
          ff = xscale(1)
          gg = 1.d0 + 10.d0*dble(ndv-1)
          do i=2,ndv
             gg = gg + xscale(i)**2 - 10.d0*cos(4.d0*pi*xscale(i))
          end do
          if(ff.lt.0.d0)ff = 0.d0
          if(gg.le.0.d0)stop'ZDT4'
          tmp = ff/gg
          if(tmp.lt.0.d0)stop'ZDT4 ff/gg'
          ynew(1) = ff
          ynew(2) = gg*(1.d0-sqrt(tmp))
          ynew(3) = 0.d0

        else if(neva.eq.-16)then
! ZDT6
          if(nobj.ne.2)stop'ZDT6, nobj'
          if(ndv.ne.10)stop'ZDT6, ndv'
          if(Cobj(1).ne.'Min'.and.Cobj(1).ne.'min')stop'ZDT6, Cobj'
          if(Cobj(2).ne.'Min'.and.Cobj(2).ne.'min')stop'ZDT6, Cobj'
          xscale(:) = xnew(:)
          ff = 1.d0 -(exp(-4.d0*xscale(1)))*((sin(6.d0*pi*xscale(1)))**6)
          gg = 0.d0
          do i=2,ndv
             gg = gg + xscale(i)
          end do
          gg = 1.d0 + 9.d0*( (gg/dble(ndv-1))**(0.25d0) )
          if(gg.le.0.d0)stop'ZDT6'
          ynew(1) = ff
          ynew(2) = gg*(1.d0-((ff/gg)**2))
          ynew(3) = 0.d0

        else
          stop'unknown set of neva'
        end if
        end subroutine testfunction
