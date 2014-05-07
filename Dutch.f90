  subroutine Dutch(Ddib,fdib,gdib,hdib,orderex,Dutchorder,Dtoex,ftoex,DIMDI,DIM,NTOEX)

    implicit none

    integer :: ii,jj,jjj,kk,j,dim2,DIMDI,DIM,&
         info,ipvt(1+DIM+DIMDI),ipvt2(1+DIM+DIM**2+DIMDI),&
         NTOEX,orderex(0:DIMDI-1),Dutchorder(NTOEX)
    real*8 :: Ddib(DIM,0:DIMDI-1),fdib(0:DIMDI-1),gdib(DIM,0:DIMDI-1),&
         hdib(DIM,DIM,0:DIMDI-1),COOR(DIM)
    real*8 :: Dtoex(DIM,NTOEX),ftoex(NTOEX),fac1,fac2,fac3
    real*8 :: MTayl(0:DIMDI-1),term(DIM),&
         Css(DIM,DIM),alpha(DIM),x0(DIM),alpha0

    logical :: calclin,calcquad

    dim2=DIMDI-1

    x0(1:DIM)=Ddib(1:DIM,0)

    Css(:,:)=0.0
    do jj=1,DIM
       do kk=1,DIM
          Css(kk,jj)=Ddib(kk,jj)-x0(kk)
       end do
    end do

    call dgefa(Css,DIM,DIM,ipvt,info)

    if (info.ne.0) then
       write(*,*) 'Zero encountered in Dutch row',info
       stop
    end if
    
    do j=1,NTOEX

       if (Dutchorder(j).eq.1) then
          fac1=0.5
          fac2=2.0/3.0
          fac3=1.0/6.0
       else if (Dutchorder(j).eq.2) then
          fac1=1.0/3.0
          fac2=0.5
          fac3=1.0/12.0
       end if
  
       COOR(1:DIM)=Dtoex(1:DIM,j)
 
       do ii=0,dim2
             
          term(1:DIM)=COOR(1:DIM)-Ddib(1:DIM,ii)

          MTayl(ii)=fdib(ii)

          if (orderex(ii).eq.1) then
             do jj=1,DIM
                MTayl(ii)=MTayl(ii) + fac1 * term(jj)*gdib(jj,ii)
             end do
          else if (orderex(ii).eq.2) then
             do jj=1,DIM
                MTayl(ii)=MTayl(ii) + fac2 * term(jj)*gdib(jj,ii)
                do jjj=1,DIM
                  MTayl(ii)=MTayl(ii) + fac3 * hdib(jj,jjj,ii)*term(jj)*term(jjj)
               end do
             end do
          end if
          
       end do

       alpha(1:DIM)=COOR(1:DIM)-x0(1:DIM)
 
       call dgesl(Css,DIM,DIM,ipvt,alpha(1:DIM),0) 
       
       alpha0=1.0
       do jj=1,DIM
          alpha0=alpha0-alpha(jj)
       end do
 
       ftoex(j)=alpha0*MTayl(0)
       do jj=1,DIM
          ftoex(j)=ftoex(j)+alpha(jj)*MTayl(jj)
       end do
        
          
    end do
    


  end subroutine Dutch

