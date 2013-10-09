        subroutine find_Optimal
        use dimKrig
        implicit none
! IDopt is for sample(nsample,:)
! Assuming Minimization 
        integer :: i,j,istat,l
        double precision :: pena
        double precision, dimension(nfunc) :: ycon
        character(len=2) :: Cl

        IDopt = 0
        OFopt = 1.d10
        lopt  = 100
        do 100 i=1,nsample
           do l=1,9
             write(Cl,101)l
             if(info(i)(1:2).eq.Cl)then
               lopt = min(l,lopt)
             end if
           end do
100     continue
        write(Cl,101)lopt
101     format(i1,"_")


        istat = 0
        do 200 i=1,nsample
           if(info(i)(1:2).ne.Cl)go to 200
           if(func(i,nfunc).ne.0.d0)go to 200

           if(istat.eq.0)istat = 1
           ycon(:) = func(i,:)
           call constraint(ycon,pena)
           if(pena.ne.0.d0)go to 200

           istat = 2
           do j=1,nOPT
             if( OFopt(j).gt.func(i,nfOPT(j)) )then
               IDopt(j) = i
               OFopt(j) = func(i,nfOPT(j))
             end if
           end do
200     continue

!       if(istat.eq.0.and.Cmode(:6).ne.'Update')then
        if(istat.eq.0)then
          write(*,*)'*No Practical Samples?'
          call stop_all
        else if(istat.eq.1)then
          if(id_proc.eq.0)then
            write(filenum,'(6x,a)')'>> No design satisfies all constraints'
          end if
        else if(istat.eq.2)then
           if(id_proc.eq.0)then
              write(filenum,*)
              write(filenum,'(1x,a,i4,a,f12.7)') '>> Best design: Sample point',IDopt(1),' Fct value',OFopt(1)
              write(filenum,'(1x,a,50f12.7)')    '   with X value',sample(IDopt(1),:)
              write(filenum,*)
           end if
        end if
    
        end subroutine find_Optimal
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine constraint(ycon,pena)
        use dimKrig
        implicit none
        double precision, dimension(nfunc), intent(in)  :: ycon
        double precision,                   intent(out) :: pena
        integer :: i
        double precision :: val

        pena = 0.d0
        do i=1,nCON
           val = ycon(nfCON(i))
           if(     cCON1(i).eq.'>'.and.cCON2(i).eq.'add_one')then
              if(val.lt.vCON(i))pena = pena + 1.d0

           else if(cCON1(i).eq.'<'.and.cCON2(i).eq.'add_one')then
              if(val.gt.vCON(i))pena = pena + 1.d0

           else if(cCON1(i).eq.'>'.and.cCON2(i).eq.'add_dif')then
              if(val.lt.vCON(i))pena = pena + (vCON(i)-val)

           else if(cCON1(i).eq.'<'.and.cCON2(i).eq.'add_dif')then
              if(val.gt.vCON(i))pena = pena + (val-vCON(i))

           else
             write(*,*)'*constraint',cCON1(i),cCON2(i)
           end if
        end do
        end subroutine constraint
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine dist_other_sample(xnew,new,ndes,design,ddmin)
        use dimKrig
        implicit none
        double precision, dimension(ndim), intent(in) :: xnew
        integer, intent(in) :: new,ndes
        double precision, dimension(ndes,ndim), intent(in) :: design
        double precision, intent(out) :: ddmin
        integer :: i,k
        double precision :: dd

        ddmin = 1.d5
! just before
        do i=1,new-1
           dd = 0.d0
           do k=1,ndim
             dd = dd + ( xnew(k)-design(i,k) )**2
           end do
           dd    = dsqrt(dd)
           ddmin = min(dd,ddmin)
        end do
        do i=1,nsample
           dd = 0.d0
           do k=1,ndim
             dd = dd + ( xnew(k)-sample(i,k) )**2
           end do
           dd    = dsqrt(dd)
           ddmin = min(dd,ddmin)
        end do


        end subroutine dist_other_sample
