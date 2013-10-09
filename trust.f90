        subroutine Check_Trust_Region
        use dimKrig
        implicit none
        double precision, dimension(ndim)      :: xnew
        double precision, dimension(nfunc)     :: ynew
        integer, dimension(nsample)            :: nlist,near
        double precision, dimension(nsample)   :: dnear
        double precision, dimension(1,1)       :: d1y,d2y
        double precision, dimension(1,ndim)    :: vdx
        double precision, dimension(ndim,1)    :: vg
        double precision, dimension(ndim,ndim) :: vh

        character(len=1000000) :: Cline
        integer :: ndat
        character(len=9) :: Class
        integer :: i,j,k
        integer :: ilist,igrad,noone,ifunc
        double precision :: pk1,pk2,ycent,yaprx,yreal,dx,dd,dxnew
        character(len=6) :: mode

        if(itrust.eq.0)stop'itrust.eq.0'
! read new sample
        open(10,file='newsample.dat',form='formatted',status='old')
        read(10,*)mode
        read(10,'(a1000000)')Cline
        call Count_Data(Cline,ndat)
        backspace 10
        if(ndat.eq.(ndim+nfunc+nCOK*ndim+nCOK*ndim*ndim).and.nCOK.ne.0)then
          if(mode(3:6).ne.'FGH ')stop'check_trust'
        else if(ndat.eq.(ndim+nfunc+nCOK*ndim+nCOK*nhes).and.nCOK.ne.0)then
          if(mode(3:6).ne.'FGH ')stop'check_trust'
        else if(ndat.eq.(ndim+nfunc+nCOK*ndim).and.nCOK.ne.0)then
          if(mode(3:6).ne.'FG  ')stop'check_trust'
        else if(ndat.eq.(ndim+nfunc))then
          if(mode(3:6).ne.'F   ')stop'check_trust'
        else
          stop'unknown set of newsample.dat'
        end if
        read(10,*)(xnew(i),i=1,ndim),(ynew(i),i=1,nfunc)
        close(10)
        if(ynew(nfunc).ne.0.d0)stop'Dummy Sample'

! start trust region modification
        call find_Optimal
        do 100 ifunc=1,nfunc-1
           call reduce_data(ifunc,igrad,0.d0,0.d0)
           if(igrad.eq.2)then
!            call indirect(ifunc)
           else
             go to 100
           end if

           ilist = 0
           nlist = 0
           near  = 0
           dnear = 0.d0
           do 200 i=1,rsample
              dd = 0.d0
              do k=1,ndim
                dd = dd + (xnew(k)-sampl(i,k))**2
              end do
              dd = dsqrt(dd)
! within region
              if(dd.le.tdx(i)*mgfac)then
                ilist = ilist + 1
                nlist(ilist)  = i
              end if
! nearest
              if(i.eq.1)then
                 near(1) = i
                dnear(1) = dd
                go to 200
              else
                do 210 j=1,i-1
                  if(dd.gt.dnear(j))go to 210
                  do 220 k=i-1,j,-1
                     near(k+1) =  near(k)
                    dnear(k+1) = dnear(k)
220               continue
                   near(j) = i
                  dnear(j) = dd
                  go to 200
210             continue
                 near(i) = i
                dnear(i) = dd
              end if
200        continue
           do i=1,rsample-1
             if(dnear(i).gt.dnear(i+1))stop'dnear miss'
           end do

! Target Specification
           noone = 0
           if(ilist.eq.0)then
             noone = 1
             do 250 i=1,rsample
                ilist = ilist + 1
                nlist(ilist) = near(i)
                if(ilist.eq.3)go to 260
250          continue
           end if
260        continue
           if(ilist.eq.0)then
             write(*,'(1x,a)') &
             '>> No Target Sample Points for Trust Region'
           else if(noone.eq.1)then
             write(*,'(1x,a,i3,a)') &
             '>> New Pt is Outside of Any Trust Region, Try to',ilist, &
             ' Nearest Points'
           else
             write(*,'(1x,a,i3,a)') &
             '>> New Pt is Inside of',ilist,' Trust Regions'
           end if

! Loop for Modification of Trust Region
           do 300 i=1,ilist
              dx       = tdx(nlist(i))
              vdx(1,:) = xnew(:) - sampl(nlist(i),:)
              vg(:,1)  = gfun(nlist(i),:)
              vh(:,:)  = hfun(nlist(i),:,:)
              ycent    =  fun(nlist(i))
              if(inf(nlist(i))(3:6).eq.'F   '.and.dx.eq.0.d0)then
                write(*,*)'>> This Point does not have Trust Region'
                go to 300
              else if(inf(nlist(i))(3:6).eq.'FG  ')then
                d1y = matmul(vdx,vg)
                d2y = 0.d0
              else if(inf(nlist(i))(3:6).eq.'FGH ')then
                d1y = matmul(vdx,vg)
                d2y = matmul( matmul(vdx,vh),transpose(vdx) )
              else
                stop'trust region for Hvec?'
              end if
              yaprx = ycent + d1y(1,1) + 0.5d0*d2y(1,1)
              yreal = ynew(ifunc)
! Evaluation
              if(ycent.eq.yaprx)then
                pk1 = 0.d0
                if(yreal.ne.0.d0)then
                  pk2 = dabs((yaprx/yreal)-1.d0)
                else if(yaprx.ne.0.d0)then
                  pk2 = dabs((yreal/yaprx)-1.d0)
                else
                  pk2 = 0.d0
                end if
              else
                pk1 = (ycent-yreal)/(ycent-yaprx)
                pk2 = dabs(pk1-1.d0)
              end if
              if(pk1.lt.0.d0)then
                 Class = "Too Bad"
              else
                 if(     pk2.ge.0.d0.and.pk2.le.pkl)then
                   Class = "Very Good"
                 else if(pk2.gt.pkl .and.pk2.le.pku)then
                   Class = "Good"
                 else if(pk2.gt.pku                )then
                   Class = "Bad"
                 else
                   stop'range of pk2'
                 end if
              end if
              write(*,'(6x,a,2(f8.3,a),a)') &
              '>> (pk1,pk2) = (',pk1,' ,',pk2,' ), Category : ',Class

! Modification
              dd = 0.d0
              do k=1,ndim
                dd = dd + vdx(1,k)**2
              end do
              dd = dsqrt(dd)
              if(noone.eq.1)then ! The newpt is outside of any region
                if(dd.le.dx*mgfac)stop'why?'
                if(Class.eq.'Very Good')then
                    write(*,'(11x,a,2(f7.3,a),i5)') &
                    '>> Trust Region Augmentation to NewPt  ', &
                    dx,' ->',dd,' for Pt',nlist(i)
                    tdx(nlist(i)) = dd
                else if(Class.eq.'Good')then
                    write(*,'(11x,a,2(f7.3,a),i5)') &
                    '>> Trust Region Augmentation to NewPt/2', &
                    dx,' ->',dx+0.5*(dd-dx),' for Pt',nlist(i)
                    tdx(nlist(i)) = dx + 0.5d0*(dd-dx)
                end if
              else               ! The newpt is inside of trust regions
                if(dd.gt.dx*mgfac)stop'why??'
                if(Class.eq.'Too Bad'.or.Class.eq.'Bad')then
                    write(*,'(11x,a,f3.1,2(f7.3,a),i5)') &
                    '>> Trust Region Reduction to dd*',bdfac, &
                    dx,' ->',dd*bdfac,' for Pt',nlist(i)
                    tdx(nlist(i)) = dd*bdfac
                else if(Class.eq.'Good')then
                    write(*,'(11x,a,15x,f7.3,a,i5)') &
                    '>> Trust Region Preservation  ', &
                    dx,' for Pt',nlist(i)
                else if(Class.eq.'Very Good')then
                    write(*,'(11x,a,f3.1,2(f7.3,a),i5)') &
                    '>> Trust Region Expansion to dd*',vgfac, &
                    dx,' ->',dx*vgfac,' for Pt',nlist(i)
                    tdx(nlist(i)) = dx*vgfac
                else
                  stop'unknown Class'
                end if
              end if
300        continue ! ilist loop(i)

! Decision of Trust Region of New Pt
           dxnew = 1000.d0
           if(mode(3:6).eq.'F   ')then
              dxnew = 0.d0
           else if(mode(3:6).eq.'FG  '.or.mode(3:6).eq.'FGH ')then
             if(noone.eq.1)then ! outside
               do 400 i=1,rsample
                  if(inf(near(i))(3:6).eq.mode(3:6))then
                    dxnew = tdx(near(i))
                    go to 410
                  end if
400            continue
410            continue
             else               ! inside
               do 500 i=1,ilist
                  if(inf(nlist(i))(3:6).eq.mode(3:6))then
                    dxnew = min( dxnew,tdx(nlist(i)) )
                  end if
                  if(mode(3:6).eq.'FGH ' &
                .and.inf(nlist(i))(3:6).eq.'F   ' &
                .and.tdx(nlist(i)).ne.0.d0)then
                    dxnew = min( dxnew,tdx(nlist(i)) )
                    ! for grad/Hess eliminated point
                  end if
500            continue
             end if
           end if
           if(dxnew.gt.1.d0)then
             write(*,'(1x,a)') &
             '>> No corresponding trust region alongside'
             if(mode(3:6).eq.'FG  ')then
               dxnew = tdxinit(1)
             else if(mode(3:6).eq.'FGH ')then
               dxnew = tdxinit(2)
             end if
           end if
           write(*,'(1x,a,f12.7)') &
           '>> Trust Region Set for New Point = ',dxnew
           if(mode(3:6).ne.'F   ')then
             if(dxnew.lt.0.d0.or.dxnew.gt.1.d0)stop'the range of dxnew'
           end if

! return to tdxx
           do i=1,rsample
              j = icone(i)
              tdxx(j) = tdx(i)
           end do
! output the new set
           call write_trustdx(ifunc,dxnew)
100     continue ! ifunc loop

        end subroutine Check_Trust_Region
