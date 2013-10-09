        subroutine ANOVA(ifac,ndim,ANV)
        implicit none
        integer, intent(in) :: ifac,ndim
        double precision, dimension(ifac**ndim,1) :: ANV

        integer,          dimension(ndim) :: itp
        double precision, dimension(ndim) :: x
        integer :: ivol,idv,jdv
        integer :: i,j,k
        double precision :: dx,total,remain
        double precision :: mean,vari
        double precision, dimension(    ndim,0:ifac-1) :: mean_i
        double precision, dimension(0:ifac-1,0:ifac-1) :: mean_ij
        double precision, dimension(ndim)              :: vari_i
        double precision, dimension(ndim,ndim)         :: vari_ij

        ivol = 1  ! define the basic size of design space
        write(*,'(1x,a)')'>> Start to ANOVA'

        mean = 0.d0
        do 100 i=1,ifac**ndim
           call make_cartesian_higherD(i,ndim,ifac,itp,x)
           if(ivol.eq.0)then
             dx = (1.d0/dble(ifac-1))**(ndim)
           else
             dx = 1.d0
           end if
           do k=1,ndim
             if(itp(k).eq.0.or.itp(k).eq.ifac-1)then
               dx = dx * 0.5d0
             end if
           end do
           mean = mean + ANV(i,1)*dx
100     continue
        if(ivol.eq.0)then
        else
          mean = mean / (dble(ifac-1)**ndim)
        end if

        vari = 0.d0
        do 200 i=1,ifac**ndim
           call make_cartesian_higherD(i,ndim,ifac,itp,x)
           if(ivol.eq.0)then
             dx = (1.d0/dble(ifac-1))**(ndim)
           else
             dx = 1.d0
           end if
           do k=1,ndim
             if(itp(k).eq.0.or.itp(k).eq.ifac-1)then
               dx = dx * 0.5d0
             end if
           end do
           vari = vari + ((ANV(i,1)-mean)**2)*dx
200     continue
        if(ivol.eq.0)then
        else
          vari = vari / (dble(ifac-1)**ndim)
        end if
        write(*,'(6x,a,2e15.5)')'>> Mean and Variance = ',mean,vari

        mean_i = 0.d0
        vari_i = 0.d0
        do 300 idv=1,ndim
          do 310 i=1,ifac**ndim
           call make_cartesian_higherD(i,ndim,ifac,itp,x)
           if(ivol.eq.0)then
             dx = (1.d0/dble(ifac-1))**(ndim-1)
           else
             dx = 1.d0
           end if
           do k=1,ndim
             if(k.ne.idv)then
             if(itp(k).eq.0.or.itp(k).eq.ifac-1)then
               dx = dx * 0.5d0
             end if
             end if
           end do
           mean_i(idv,itp(idv)) = mean_i(idv,itp(idv)) + ANV(i,1)*dx
310       continue
          do 320 i=0,ifac-1
            if(ivol.eq.0)then
              mean_i(idv,i) = mean_i(idv,i) - mean
              dx = 1.d0 / dble(ifac-1)
            else
              mean_i(idv,i) = mean_i(idv,i)/(dble(ifac-1)**(ndim-1)) - mean
              dx = 1.d0
            end if
            if(i.eq.0.or.i.eq.ifac-1) dx = dx*0.5d0
            vari_i(idv) = vari_i(idv) + ((mean_i(idv,i))**2)*dx
320       continue
          if(ivol.eq.0)then
          else
            vari_i(idv) = vari_i(idv) / dble(ifac-1)
          end if
300     continue

        vari_ij = 0.d0
        do 400 idv =     1,ndim-1
        do 400 jdv = idv+1,ndim
           if(idv.eq.jdv)go to 400
           mean_ij = 0.d0
           do 410 i=1,ifac**ndim
             call make_cartesian_higherD(i,ndim,ifac,itp,x)
             if(ivol.eq.0)then
               dx = (1.d0/dble(ifac-1))**(ndim-2)
             else
               dx = 1.d0
             end if
             do k=1,ndim
               if(k.ne.idv.and.k.ne.jdv)then
               if(itp(k).eq.0.or.itp(k).eq.ifac-1)then
                 dx = dx * 0.5d0
               end if
               end if
             end do
             mean_ij(itp(idv),itp(jdv)) = &
             mean_ij(itp(idv),itp(jdv)) + ANV(i,1)*dx
410        continue ! main loop (i)
           do 420 i=0,ifac-1
           do 420 j=0,ifac-1
             if(ivol.eq.0)then
               dx = (1.d0/dble(ifac-1))**2
             else
               mean_ij(i,j) = mean_ij(i,j)/(dble(ifac-1)**(ndim-2))
               dx = 1.d0
             end if
             if(i.eq.0.or.i.eq.ifac-1) dx = dx*0.5d0
             if(j.eq.0.or.j.eq.ifac-1) dx = dx*0.5d0
             vari_ij(idv,jdv) = vari_ij(idv,jdv) + &
             ( (mean_ij(i,j)-mean_i(idv,i)-mean_i(jdv,j)-mean)**2 )*dx
420        continue
           if(ivol.eq.0)then
           else
             vari_ij(idv,jdv) = vari_ij(idv,jdv) / (dble(ifac-1)**2)
           end if
400     continue ! twin-dv loop (idv,jdv)
 
! output
        total = 0.d0
        do 500 i=1,ndim
           vari_i(i) = vari_i(i) / vari * 100.d0
           write(*,'(6x,a,i2,a,f10.3)') &
           '>> Significance of DV(',i,' )    = ',vari_i(i)
           total = total + vari_i(i)
500     continue
        do 510 i=1,ndim-1
        do 510 j=i+1,ndim
           vari_ij(i,j) = vari_ij(i,j) / vari * 100.d0
           write(*,'(6x,a,i2,a,i2,a,f10.3)') &
           '>> Significance of DV(',i,'-',j,' ) = ',vari_ij(i,j)
           total = total + vari_ij(i,j)
510     continue
        remain = 100.d0 - total
        write(*,'(6x,a,f10.3)') &
           '>> Higher Order Terms         = ',remain
        write(*,'(36x,a)')'+)__________'
        write(*,'(41x,a)')'100.000'

        end subroutine ANOVA
