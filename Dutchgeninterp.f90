  subroutine Dutchgeninterp(Ddib,fdib,gdib,hdib,orderex,Dutchorder,Dtoex,ftoex,DIMDI,DIM,NTOEX)

    implicit none

    integer :: ii,jj,jjj,j,DIMDI,DIM,info,ipvt(DIMDI),NTOEX,orderex(0:DIMDI-1),Dutchorder(NTOEX),nterms,mreg(DIMDI,DIM)
    real*8 :: Ddib(DIM,0:DIMDI-1),fdib(0:DIMDI-1),gdib(DIM,0:DIMDI-1),hdib(DIM,DIM,0:DIMDI-1),COOR(DIM)
    real*8 :: Dtoex(DIM,NTOEX),ftoex(NTOEX),phi,fac1,fac2,fac3
    real*8 :: MTayl(0:DIMDI-1),term(DIM),Css(DIMDI,DIMDI),Fs(DIMDI)

    logical :: calclin,calcquad

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

       do ii=1,DIMDI
          Css(ii,1)=1.0
          do jj=1,DIM
             Css(ii,jj+1)=Ddib(jj,ii-1)
          end do
       end do

       call dgefa(Css,DIMDI,DIMDI,ipvt,info)
       
       if (info.ne.0) then
          write(*,*) 'Zero encountered in LIN row',info
          stop
       end if

    end if

    if (calcquad) then

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

       do ii=1,DIMDI
          do jj=1,DIMDI          
             fac1=1.0
             do jjj=1,DIM
                   fac1=fac1*Ddib(jjj,ii-1)**mreg(jj,jjj)
             end do
             Css(ii,jj)=fac1
          end do
       end do

!!$       do jj=1,DIMDI
!!$          write(*,'(20e12.4)'), (Css(jj,ii),ii=1,DIMDI)
!!$       end do
!!$       stop
       
       call dgefa(Css,DIMDI,DIMDI,ipvt,info)
       
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
 
       do ii=0,DIMDI-1
             
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

       Fs(1:DIMDI)=MTayl(0:DIMDI-1)
        
       call dgesl(Css,DIMDI,DIMDI,ipvt,Fs,0)

       if (Dutchorder(j).eq.1) then
                                   
          ftoex(j)=Fs(1)
          do jj=1,DIM
             ftoex(j)=ftoex(j)+Fs(jj+1)*COOR(jj)
          end do

       else if (Dutchorder(j).eq.2) then
   
          ftoex(j)=0.0
          do ii=1,DIMDI         
             fac1=Fs(ii)
             do jj=1,DIM
                fac1=fac1*COOR(jj)**mreg(ii,jj)
             end do
             ftoex(j)=ftoex(j)+fac1
          end do

       end if
          
    end do
    


  end subroutine Dutchgeninterp


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
