        subroutine Output_Design(mode,ndes,mdim,design)
        use dimKrig
        implicit none
        integer, intent(in) :: ndes,mdim,mode
        double precision, dimension(ndes,mdim) :: design
        integer :: i,new,iover
        double precision :: zero,pena
        double precision :: EI_thres,dd_thres,ddmin,difEI
        double precision, dimension(ndim) :: xnew
        double precision, dimension(nOPT*2+1+3*(nfunc-1)) :: ynew,yprev
        character(len=14) :: Cfile

        if(mEI.eq.0)then
          dd_thres = sqrt( dble(ndim)*((1.d-5)**2) )
        else
          dd_thres = sqrt( dble(ndim)*((1.d-5)**2) )
        end if
        EI_thres = 1.d-5

        if(ndim.ne.mdim)stop'mdim in output_design'
        do 100 new=1,ndes
          zero = 0.d0
          do i=1,ndim
            if(design(new,i).lt.0.d0.or.design(new,i).gt.1.d0)then
              if(id_proc.eq.0) &
              write(*,'(a,2i5,f15.8)') &
              '*New Design is outside of design space',new,i,design(new,i)
              go to 100
              call stop_all
            end if
            xnew(i) = design(new,i)
            zero    = zero + xnew(i)**2
          end do
          if(zero.eq.0.d0.and.id_proc.eq.0) &
          write(*,'(a)')'*New Design is at 0.0'

          call dist_other_sample(xnew,new,ndes,design,ddmin)
          if(mode.eq.15)then
            call meta_all(mode,nOPT,xnew,ynew)
            pena = ynew(nOPT+1)
          else
            call meta_all(mode,nOPT*2,xnew,ynew)
            pena = ynew(nOPT*2+1)
          end if

          if(id_proc.eq.0)then
           write(*,'(1x,a,i3)')'>> New Sample Point ------ ',new
           write(*,'(6x,a,99f8.3)')  '>> x       = ',(xnew(i),i=1,ndim)
          end if

          iover = new
! No.1 is with weak constraint to update the training data anyway
          if(new.eq.1)then
            if(ddmin.le.0.d0)then
              if(id_proc.eq.0)&
              write(*,'(6x,a)')'>> Passing by Distance = 0.d0'
              go to 100
            end if
            go to 150

! EI-based new design
          else if( (mod(new,2).eq.1.and.Cmode(11:13).eq.'GA2').or. &
                                       (Cmode(11:13).eq.'GA1').or. &
                                       (Cmode(11:13).eq.'GM2') )then
            if(ynew(new).eq.0.d0)then
              if(id_proc.eq.0)&
              write(*,'(6x,a)')'>> Passing by EI = 0.d0'
              go to 100
            end if
            if(ddmin.le.dd_thres)then
              if(id_proc.eq.0)&
              write(*,'(6x,a,e12.3)')'>> Passing by Distance = ',ddmin
              go to 100
            end if
            go to 150

! Yhat-based new design in GA
          else if(mod(new,2).eq.0.and.Cmode(11:13).eq.'GA2')then
            difEI = dabs(ynew(new-1)-yprev(new-1))
            if(ddmin.le.0.d0)then
              if(id_proc.eq.0)&
              write(*,'(6x,a)')'>> Passing by Distance = 0.d0'
              go to 100
            end if
            if(IDopt(new/2).eq.0)then
              iover = new-1
              if(id_proc.eq.0)&
              write(*,'(6x,a,i3,a)') &
              '>> Overwrite',iover,' by meaningless EI'
              go to 150
            end if
            if(difEI.le.EI_thres)then
              iover = new-1
              if(id_proc.eq.0)&
              write(*,'(6x,a,i3,a)') &
              '>> Overwrite',iover,' by comparable EI'
              go to 150
            end if
            if(ddmin.gt.0.d0.and.ddmin.le.dd_thres)then
              if(id_proc.eq.0)&
              write(*,'(6x,a,e12.3)')'>> Passing by Distance = ',ddmin
              go to 100
            end if
            go to 150
! Others/MaxVar
          else
            if(ddmin.ge.0.d0.and.ddmin.le.dd_thres)then
              if(id_proc.eq.0)&
              write(*,'(6x,a,e12.3)')'>> Passing by Distance = ',ddmin
              go to 100
            end if
            go to 150
          end if
          stop 'Something is wrong in output_des'

150       continue
          if(id_proc.eq.0)then
           if(mode.eq.15)then
             write(*,'(6x,a,99f8.3)')'>> Uncertainty = ',(ynew(i),i=1,nOPT)
           else if(nOPT.eq.1)then
             write(*,'(6x,a,2f15.8)')'>> EI/Yhat = ',(ynew(i),i=1,2)
           else
             write(*,'(6x,a,99f8.3)')'>> EI      = ',(ynew(i),i=1,nOPT*2-1,2)
             write(*,'(6x,a,99f8.3)')'>> Yhat    = ',(ynew(i),i=2,nOPT*2,2)
           end if
           if(pena.ne.0.d0)then
             write(*,'(11x,a,f8.3)')'>> Penalty Unfeasible',pena
           end if
           write(Cfile,'(a9,i1,a4)')'newsample',iover,'.dat'
           if(iover.eq.new)then
             write(*,'(11x,a,a)')'>> Output to ',Cfile
           else
             write(*,'(11x,a,a)')'>> Overwrite to ',Cfile
           end if
           open(10,file=Cfile,form='formatted',status='unknown')
           write(10,'(100f15.8)')(xnew(i),i=1,ndim),(ynew(i),i=1,2)
           close(10)
          end if
          yprev = ynew

100     continue

!!$        xnew = 0.6d0
!!$        call meta_all(mode,nOPT*2,xnew,ynew)
!!$        if(id_proc.eq.0)then
!!$          write(*,'(6x,a,2f15.8)')'>> EI/Yhat @ 0.6 = ',(ynew(i),i=1,2)
!!$        end if

        end subroutine Output_Design
