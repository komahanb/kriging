        subroutine Check_Sample
        use dimKrig
        implicit none
        integer :: i,j,k,l
        integer :: icheck
        integer :: il,jl
        double precision :: dl,dd,h1,h2
        character(len=2) :: cl

        do l=1,lmax
          write(cl,101)l
101       format(i1,'_')
          il  = 0
          jl  = 0
          dl  = 1.d10
          do i=1,nsample-1
            do j=i+1,nsample
              if(info(i)(1:2).eq.cl.and.info(j)(1:2).eq.cl)then
                dd = 0.d0
                do k=1,ndim
                   dd = dd + (sample(i,k)-sample(j,k))**2
                end do
                dd = dsqrt(dd)
                if(dd.lt.dl)then
                 il = i
                 jl = j
                 dl = dd
                end if
              end if
            end do
          end do
          if(id_proc.eq.0.and.ict_sample(l,0).gt.1)then
            write(filenum,'(6x,a,i1,2(a,i3),a,e10.2)') &
            '>> Min. Dist. between Fidelity-',l, &
            ' Samples (',il,',',jl,' ) =',dl
          end if
        end do

        icheck = 0
        do 100 i=1,nsample
        do 100 j=1,nCOK
           if(info(i)(3:6).eq.'FH  '.or.info(i)(3:6).eq.'FGH ')then
             do 110 k=1,ndim
             do 120 l=1,ndim
                if(k.eq.l)go to 120
                h1 = hfunc(i,nfCOK(j),k,l)
                h2 = hfunc(i,nfCOK(j),l,k)
                if(h1.ne.h2)then
                  write(*,'(6x,a,4i5,2e15.5)') &
                  '>> Hessian Asymmetry',i,j,k,l,h1,h2
                  icheck = 1
                end if
120          continue
110          continue
           end if
100     continue

        if(icheck.eq.1)then
          write(*,*)'*Stop by Hessian Asymmetry'
          call stop_all
        end if

        end subroutine Check_Sample
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine check
        use dimKrig
        implicit none
        double precision :: scf,dsdf,dsdb,d2sdfdb,d2sdf,d2sdb
        double precision :: xf,xb,pp,tt

        xf = 1.d0
        xb = 0.d0
        tt = 1.d0
        pp = 2.d0
        write(*,*)"*theta=1.d0, power=2.d0, input xf and xb"
        read(*,*)xf,xb
        write(*,'(a,2f15.8)')"*xf, xb = ",xf,xb
        call SCF_DF_DF(iscf,hstat,xf,1.d0,1.d0,xb,tt,pp,scf,dsdf,d2sdf)
        call SCF_DB_DB(iscf,hstat,xf,xb,1.d0,1.d0,tt,pp,scf,dsdf,d2sdb)
        call SCF_DF_DB(iscf,hstat,xf,1.d0,xb,1.d0,tt,pp,scf,dsdf,d2sdfdb)
        call SCF_DB(   iscf,hstat,xf,xb,1.d0,tt,pp,scf,dsdb)
        call SCF_OUT(  iscf,xf,xb,tt,pp,scf)
        write(*,'(i5,6f15.8)')iscf,scf,dsdf,dsdb,d2sdfdb,d2sdf,d2sdb

        call stop_all
        end subroutine check
 
