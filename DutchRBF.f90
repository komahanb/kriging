  subroutine DutchRBF(Ddib,fdib,gdib,hdib,orderex,Dutchorder,Dtoex,ftoex,DIMDI,DIM,NTOEX)

    implicit none

    integer :: ii,jj,jjj,j,dim2,DIMDI,DIM,info,ipvt(1+DIM+DIMDI),ipvt2(1+1.5*DIM+0.5*DIM**2+DIMDI),NTOEX,orderex(0:DIMDI-1),Dutchorder(NTOEX),nterms,mreg(DIMDI,DIM),dimquad,dimquad2
    real*8 :: Ddib(DIM,0:DIMDI-1),fdib(0:DIMDI-1),gdib(DIM,0:DIMDI-1),hdib(DIM,DIM,0:DIMDI-1),COOR(DIM)
    real*8 :: Dtoex(DIM,NTOEX),ftoex(NTOEX),phi,fac1,fac2,fac3
    real*8 :: MTayl(0:DIMDI-1),term(DIM),Css(1+DIM+DIMDI,1+DIM+DIMDI),Fs(1+DIM+DIMDI),Css2(1+1.5*DIM+0.5*DIM**2+DIMDI,1+1.5*DIM+0.5*DIM**2+DIMDI),Fs2(1+1.5*DIM+0.5*DIM**2+DIMDI)

    logical :: calclin,calcquad

    dim2=DIMDI-1 

    calclin=.false.
    calcquad=.false.

    do j=1,NTOEX

       if (Dutchorder(j).eq.1) then
          calclin=.true.
       else if (Dutchorder(j).eq.2) then
          calcquad=.true.
       else
          write (*,*) 'Dutchorder',Dutchorder(j),' not implemented'
          STOP
       end if

    end do

    if (calclin) then

       Css(:,:)=0.0
       
       do ii=2+DIM,1+DIM+DIMDI
          Css(1,ii)=1.0
          Css(ii,1)=1.0
          do jj=1,DIM
             Css(jj+1,ii)=Ddib(jj,ii-DIM-2)  
             Css(ii,jj+1)=Ddib(jj,ii-DIM-2)
          end do
          do jj=ii,1+DIM+DIMDI
             Css(jj,ii)=phi(DIM,Ddib(1:DIM,jj-DIM-2),Ddib(1:DIM,ii-DIM-2))
             Css(ii,jj)=Css(jj,ii)
          end do
       end do

       call dgefa(Css,1+DIM+DIMDI,1+DIM+DIMDI,ipvt,info)
       
       if (info.ne.0) then
          write(*,*) 'Zero encountered in LIN row',info
          stop
       end if

    end if

    if (calcquad) then

       dimquad=1+1.5*DIM+0.5*DIM**2+DIMDI
       dimquad2=1+1.5*DIM+0.5*DIM**2

       call basis(DIMDI,DIM,2,mreg,nterms)

       if (nterms.gt.DIMDI) then
          write(*,*) 'Need at least',nterms,'base points not', DIMDI
          STOP
       end if

!!$       write(*,'(a,3i6)')'DIM,DIMDI,nterms = ',DIM,DIMDI,nterms
!!$       do ii=1,nterms
!!$          write(*,'(i5,a,99i4)')ii,':',(mreg(ii,jj),jj=1,DIM)
!!$       end do
!!$       stop

       Css2(:,:)=0.0

       do jj=0,DIMDI-1
          do ii=1,nterms          
             fac1=1.0
             do jjj=1,DIM
                if (mreg(ii,jjj).eq.1) then
                   fac1=fac1*Ddib(jjj,jj)
                else if (mreg(ii,jjj).eq.2) then
                   fac1=fac1*Ddib(jjj,jj)**2
                end if
             end do
             Css2(ii,jj+dimquad2+1)=fac1
             Css2(jj+dimquad2+1,ii)=fac1
          end do
          do ii=0,jj
             fac1=phi(DIM,Ddib(1:DIM,jj),Ddib(1:DIM,ii))
             Css2(jj+dimquad2+1,ii+dimquad2+1)=fac1
             Css2(ii+dimquad2+1,jj+dimquad2+1)=fac1
          end do
       end do

!!$       do jj=1,dimquad
!!$          write(*,'(20e12.4)'), (Css2(jj,ii),ii=1,dimquad)
!!$       end do
       
       call dgefa(Css2,dimquad,dimquad,ipvt2,info)
       
       if (info.ne.0) then
          write(*,*) 'Zero encountered in QUAD row',info
          stop
       end if
       
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

       if (Dutchorder(j).eq.1) then
              
          Fs(1:1+DIM)=0.0
          Fs(2+DIM:2+DIM+dim2)=MTayl(0:dim2)
       
          call dgesl(Css,1+DIM+DIMDI,1+DIM+DIMDI,ipvt,Fs,0) 
          
          ftoex(j)=Fs(1)
          do jj=1,DIM
             ftoex(j)=ftoex(j)+Fs(jj+1)*COOR(jj)
          end do
          do ii=2+DIM,1+DIM+DIMDI
             ftoex(j)=ftoex(j)+Fs(ii)*phi(DIM,COOR(1:DIM),Ddib(1:DIM,ii-DIM-2))
          end do

       else if (Dutchorder(j).eq.2) then
   
          Fs2(1:dimquad2)=0.0
          Fs2(dimquad2+1:dimquad)=MTayl(0:dim2)
          
          call dgesl(Css2,dimquad,dimquad,ipvt2,Fs2,0) 

          ftoex(j)=0.0
          do ii=1,nterms          
             fac1=Fs2(ii)
             do jjj=1,DIM
                if (mreg(ii,jjj).eq.1) then
                   fac1=fac1*COOR(jjj)
                else if (mreg(ii,jjj).eq.2) then
                   fac1=fac1*COOR(jjj)**2
                end if
             end do
             ftoex(j)=ftoex(j)+fac1
          end do
          do ii=0,DIMDI-1
             fac1=Fs2(ii+dimquad2+1)*phi(DIM,COOR(1:DIM),Ddib(1:DIM,ii))
             ftoex(j)=ftoex(j)+fac1
          end do

       end if
          
    end do
    


  end subroutine DutchRBF



  FUNCTION phi(DIM,x1,x2)
    implicit none
    
    integer :: DIM,j
    real*8 :: x1(DIM),x2(DIM),phi,xnorm,tmp

    xnorm=0.0
    do j=1,DIM
       xnorm=xnorm+(x1(j)-x2(j))**2
    end do
    xnorm=SQRT(xnorm)

    !tmp=max(1.0-xnorm,0.0)
    !phi=tmp**6*(35.0*xnorm**2+18.0*xnorm+3.0)
    !phi=tmp**4*(4.0*xnorm+1.0)
    phi=xnorm**3
    !phi=xnorm


    RETURN
  END FUNCTION phi


  subroutine basis(maxdim,DIM,DIMPC,mregout,nterms)

      implicit none
      integer :: i,ii,j,k,isum,DIM,DIMPC,nterms,ent,maxdim
      integer :: mreg(maxdim,DIM),mregout(maxdim,DIM)
      
      call combination(DIM+DIMPC,DIM,nterms)    
                       
      mreg(:,:) = 0 

      if (DIMPC.ne.0) then

         do 100 i=1,nterms
            isum = 0
            do j=1,DIM
               isum = isum + mreg(i,j)
            end do
            if(isum.ne.DIMPC)then
               mreg(i+1,:) = mreg(i,:)
               mreg(i+1,1) = mreg(i+1,1) + 1
               go to 100
            else
               do j=1,DIM
                  if(mreg(i,j).ne.0)then
                     if(j.eq.DIM) go to 200
                     mreg(i+1,:) = mreg(i,:)
                     mreg(i+1,j) = 0
                     mreg(i+1,j+1) = mreg(i+1,j+1) + 1
                     go to 100
                  end if
               end do
               stop 'No target j in Make_Mreg'
            end if
 100     continue

         write(*,'(3i6)')DIMPC,nterms,DIM
         stop 'Error in Make_Mreg'
 200     continue
        
         if(i.ne.nterms)then
            write(*,*) i,nterms,DIMPC
            stop 'i.ne.nterms in Make_Mreg'
         end if

      end if


!     Resort

      mregout(:,:) = 0 
      ent=2
      do ii=1,DIMPC
         do i=2,nterms
         
            isum = 0
            do j=1,DIM
               isum = isum + mreg(i,j)
            end do
            
            if (isum.eq.ii) then
               mregout(ent,:)=mreg(i,:)
               mreg(i,:)=0
               ent=ent+1
            end if
        
         end do
      end do
   
    end subroutine
