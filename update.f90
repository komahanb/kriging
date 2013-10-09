        subroutine Update
        use dimKrig
        implicit none
        double precision, dimension(ndim)           :: xnew
        double precision, dimension(nfunc)          :: ynew
        double precision, dimension(nCOK,ndim)      :: gnew
        double precision, dimension(nCOK,ndim,ndim) :: hnew
        double precision, dimension(nCOK,ndim,2)    :: vnew
        double precision :: dummy
        integer :: i,j,k,l

        open(10,file='newsample1.dat',form='formatted',status='old')
        read(10,*) (xnew(i),i=1,ndim), dummy  
        close(10)

        call evalfunc(xnew,ndim,fct,0,hstat,ynew,gnew,hnew,vnew)

        if(id_proc.eq.0)then

          write(filenum,'(1x,2(a,i5),a)')'>> Update from',zsample,' to',zsample+1, ' sample points'
          open(11,file='sample.dat',form='formatted',status='unknown')
          write(11,'(3i6)')ndim,zsample+1,nfunc

          do i=1,zsample

             if (info(i)(3:6).ne.'FGHv') then
                write(11,101) info(i),(sample(i,j),j=1,ndim),(func(i,j),j=1,nfunc),((gfunc(i,nfCOK(j),k),k=1,ndim),j=1,nCOK),(((hfunc(i,nfCOK(j),k,l),l=1,ndim),k=1,ndim),j=1,nCOK)
             else
                write(11,101) info(i),(sample(i,j),j=1,ndim),(func(i,j),j=1,nfunc),((gfunc(i,nfCOK(j),k),k=1,ndim),j=1,nCOK),((hvect(i,nfCOK(j),k,1),k=1,ndim),j=1,nCOK),((hvect(i,nfCOK(j),k,2),k=1,ndim),j=1,nCOK)
             end if

          end do

          if( hstat.le.3 )then
           write(11,102) hstat,                    &
                      (xnew(i),i=1,ndim),                        &
                      (ynew(i),i=1,nfunc),                       &
                     ((gnew(i,j),j=1,ndim),i=1,nCOK),            &
                    (((hnew(i,j,k),k=1,ndim),j=1,ndim),i=1,nCOK)
          else 
           write(11,102)hstat,                     &
                      (xnew(i),i=1,ndim),                        &
                      (ynew(i),i=1,nfunc),                       &
                     ((gnew(i,j),j=1,ndim),i=1,nCOK),            &
                    (((vnew(i,j,k),k=1,2),j=1,ndim),i=1,nCOK),   &
                    (((hnew(i,j,k),k=1,ndim),j=1,ndim),i=1,nCOK)
          end if
          close(11)

        end if

101     format(a,10000e20.10)
102     format(i1,10000e20.10)

        end subroutine Update
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine Count_Data(Cline,ndat)
        implicit none
        character(len=1000000), intent(in) :: Cline
        integer, intent(out) :: ndat
        integer :: leng,lnew
        ndat = 0
        leng = index(Cline,'                ') - 1
100     continue
        lnew = index(Cline(:leng),'.',back=.true.)
        if(lnew.le.0)go to 200
        ndat = ndat + 1
        leng = lnew - 1
        go to 100
200     continue
        end subroutine Count_Data
