subroutine Post_Process
  use dimKrig
  implicit none
  integer :: ifac
  character*60 :: export


  if(     Cmode(:7).eq.'Post_1D')then
     if(ndim.ne.1)stop'Post_1D, ndim != 1'
     ifac = 1001

  else if(Cmode(:7).eq.'Post_2D')then
     if(ndim.ne.2)stop'Post_2D, ndim != 2'

     if (fct.eq.20) then
        ifac = 51 ! I have database of only  51*51 for CFD
     else
        ifac=101 !101
     end if

  else if(Cmode(:9).eq.'Post_RMSE')then
     if(ndim.ne.2)stop'Post_RMSE, ndim != 2'
     ifac = 0

  else if(Cmode(:12).eq.'Post_HigherD')then
     ifac = int( exp(log(2.d6)/dble(ndim)))
     if(ifac.le.2)ifac = 2
     if (ndim.eq.5) ifac=10 !10^5
     if (ndim.eq.3)  ifac = 30 !4^9 2,48,882

  else if(Cmode(:10).eq.'Post_ANOVA')then
     ifac = int( exp(log(1.d6)/dble(ndim)))
     if(ifac.le.2)ifac = 2

  else if(Cmode(:15).eq.'Post_MonteCarlo')then
  else
     if(id_proc.eq.0)then
        write(*,*)'Unknown Cmode = ',Cmode
     end if
     call stop_all
  end if

  call find_Optimal
  call read_all_krig

  if(id_proc.eq.0.and.Cmode(:15).ne.'Post_MonteCarlo')&
       write(*,'(1x,a,i10)')'>> PostProcess on Metamodel, ifac=',ifac

  if(Cmode(:7).eq.'Post_1D'.or.Cmode(:7).eq.'Post_2D')then
     call Post_1or2D(ifac)
  else if(Cmode(:9).eq.'Post_RMSE')then
     call Post_RMSE
  else if(Cmode(:15).eq.'Post_MonteCarlo')then
     call MonteCarlo
  else
     call Post_Higher(ifac)
  end if

  call deallocate_all_krig
end subroutine Post_Process
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine Post_RMSE
  use dimKrig
  implicit none
  integer :: km,ka
  double precision, dimension(100) :: dm,da
  double precision, dimension(100,100,7) :: ad
  integer :: i,j,k,ict
  double precision :: diff,yhat,RMSE,EI
  double precision, dimension(ndim) :: x,yhatprime
  double precision, dimension(nfunc) :: err
  character*60 :: export

  open(10,file='info.dat',form='formatted',status='unknown')
  read(10,*)km,ka
  if(km.ge.100.or.ka.ge.100)stop'km,ka>100'
  read(10,*)(dm(k),k=1,km)
  read(10,*)(da(k),k=1,ka)
  close(10)
  open(10,file='load.dat',form='formatted',status='unknown')
  do i=1,km
     do j=1,ka
        read(10,*)(ad(i,j,k),k=1,2),(ad(i,j,k),k=5,7)
        ad(i,j,3) = (ad(i,j,1)-dm(1))/(dm(km)-dm(1))
        ad(i,j,4) = (ad(i,j,2)-da(1))/(da(ka)-da(1))
     end do
  end do
  close(10)

  do 100 k=1,nfunc-1
     ict  = 0
     diff = 0.d0
     do 200 i=1,km
        do 200 j=1,ka
           x(1) = ad(i,j,3)
           x(2) = ad(i,j,4)
           call meta_call(k,0,x,yhat,yhatprime,RMSE,EI)
           ict = ict + 1
           diff = diff + (yhat-ad(i,j,4+k))**2
200        continue ! data loop (i,j)
           diff = dsqrt(diff/dble(ict))
           if(id_proc.eq.0) &
                write(*,'(6x,a,i1,a,e15.5)')'>> RMSE of Func-',k,' = ',diff
           err(k) = diff
100        continue ! func loop (k)

           if(id_proc.eq.0)then
              open(10,file='error.dat',form='formatted',status='unknown')
              write(10,'(i5,99e20.10)')nsample,(err(k),k=1,nfunc-1)
              close(10)
           end if

         end subroutine Post_RMSE
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         subroutine Post_1or2D(ifac)
           use dimKrig
           implicit none
           integer, intent(in) :: ifac
           double precision, allocatable, dimension(:,:)   :: TEC1
           double precision, allocatable, dimension(:,:,:) :: TEC2
           double precision, dimension(ndim)      :: xin,xe,xy,yhatprime
           double precision, dimension(ndim)      :: df,v
           double precision, dimension(ndim,ndim) :: d2f
           double precision, dimension(ndim,rsample) :: xvis
           double precision, dimension(rsample,1)      :: fvis
           integer :: i,j,k,l,negr,igrad,ict,imax,jmax,readfromdatabase
           double precision :: yhat,RMSE,EI,diff,f,scal,xtmp1,xtmp2,maxerror
           double precision :: EImax,Ymin
           character(len=13) :: Csamp,Ckrig,Cexac
           character(len=14) :: Csampl

           double precision, DIMENSION(ifac**2) :: SIGV
           double precision, DIMENSION(ifac**2) :: SIGG
           double precision, DIMENSION(ifac**2,nfunc) :: FX
           double precision, DIMENSION(ifac**2) :: SIGMA
           double precision, dimension(ndim,rsample,1) :: fgvis
           double precision, DIMENSION(ndim,ifac**2) :: XSAMP
           integer :: Taylororder, IERR,iii,jjj
           double precision :: BETA, GAMM
           character*60 :: export


           imax=ifac
           !imax=41
           jmax=ifac


           if (fct.eq.20) then
              readfromdatabase=1 ! 0=call evalfunc, 1= read from tecex10.dat
           else
              readfromdatabase=0 ! 0=call evalfunc, 1= read from tecex10.dat
           end if


           if(     Cmode(:7).eq.'Post_1D')then
              allocate(TEC1(imax,5))
              TEC1 = 0.d0
           else if(Cmode(:7).eq.'Post_2D')then
              allocate(TEC2(imax,jmax,6))
              TEC2 = 0.d0
           end if

           do k=1,1!nfunc-1

              negr  =  0
              EImax = -1.d0
              Ymin  =  1.d10
              diff  =  0.d0
              maxerror = 0.0

              if(id_proc.eq.0)then

!!$          if (fct.eq.11) then
!!$             write(Cexac, '(a,i1,a)')'tecex0',k,'.dat'
!!$             open(10,file=Cexac,form='formatted',status='old')
!!$             read(10,*)
!!$             read(10,*)
!!$             read(10,*)
!!$             read(10,*)((TEC2(i,j,1),i=1,imax),j=1,jmax)
!!$             read(10,*)((TEC2(i,j,2),i=1,imax),j=1,jmax)
!!$             read(10,*)((TEC2(i,j,6),i=1,imax),j=1,jmax)
!!$             close(10)
!!$          end if




                 if (readfromdatabase.eq.1) then  
                    write(*,*) '>>  Reading from the existing database'
                    if (fct.LE.20) then
                       write(Cexac, '(a,i1,a)')'tecex20.dat' !!! CHANGED FROM K TO FCT
                       open(15,file='tec/'//Cexac,form='formatted',status='old')
                       read(15,*)
                       read(15,*)
                       read(15,*)
                       read(15,*)((TEC2(iii,jjj,1),iii=1,imax),jjj=1,jmax) !CHECK IF ITS IN DEGREES
                       read(15,*)((TEC2(iii,jjj,2),iii=1,imax),jjj=1,jmax)
                       read(15,*)((TEC2(iii,jjj,6),iii=1,imax),jjj=1,jmax)
                       close(15)                            
                    end if

                    do 210 i=1,imax
                       do 220 j=1,jmax

                          xin(1) = dble(i-1)/dble(imax-1)
                          if(ndim.eq.2) xin(2) = dble(j-1)/dble(jmax-1)
                          call meta_call(k,1,xin,yhat,yhatprime,RMSE,EI)

                          if(rmse.lt.0.d0)negr = negr + 1

                          if(EI.ge.EImax)then
                             EImax = EI
                             xe(:) = xin(:)
                          end if

                          if(yhat.le.Ymin)then
                             Ymin  = yhat
                             xy(:) = xin(:)
                          end if

                          if(Cmode(:7).eq.'Post_1D')then
                             TEC1(i,1) = xin(1)*(DS(2,1)-DS(1,1))+DS(1,1)
                             TEC1(i,2) = yhat
                             TEC1(i,3) = RMSE
                             TEC1(i,4) = EI
                             if (fct.lt.20) then
                                call evalfunc(xin,ndim,fct,0,0,f,df,d2f,v)
                                TEC1(i,5) = f
                                diff = diff + (f-yhat)**2
                                if (abs(f-yhat).gt.maxerror) maxerror=abs(f-yhat)
                             end if
                             go to 210

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          else if(Cmode(:7).eq.'Post_2D')then


                             ! TEC2(i,j,1) = xin(1)*(DS(2,1)-DS(1,1))+DS(1,1)
                             ! TEC2(i,j,2) = xin(2)*(DS(2,2)-DS(1,2))+DS(1,2)

                             TEC2(i,j,3) = yhat
                             TEC2(i,j,4) = RMSE
                             TEC2(i,j,5) = EI
                             if (fct.le.20) then                              
                                f=  TEC2(i,j,6)
                                diff = diff + (f-yhat)**2
                                if (abs(f-yhat).gt.maxerror) maxerror=abs(f-yhat)
                             end if
                             !                         if (fct.ge.10) then
                             !                             TEC2(i,j,1) = TEC2(i,j,1) * 180.0 / 4.0 /atan(1.0)   !RADINS TO degree
                             !                           f=TEC2(i,j,6)
                             !                           diff = diff + (f-yhat)**2
                             !                           if (abs(f-yhat).gt.maxerror) maxerror=abs(f-yhat)
                             !                        end if
                             !end if
                          end if

                          !   TEC2(i,j,6) = f
                          !   diff = diff + (f-yhat)**2

                          !   if (abs(f-yhat).gt.maxerror) maxerror=abs(f-yhat)

                          !  if (readfromdatabase.eq.0) then
                          !        if (k.eq.jmax-1) write(*,*) 'EVALFUNC is called',imax,'*',jmax,'times'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


220                       continue ! grid(j)
210                       continue ! grid(i)





                       else !calling exacfunc


                          do 110 i=1,imax
                             do 120 j=1,jmax

                                xin(1) = dble(i-1)/dble(imax-1)
                                if(ndim.eq.2) xin(2) = dble(j-1)/dble(jmax-1)

                                call meta_call(k,1,xin,yhat,yhatprime,RMSE,EI)

                                if(rmse.lt.0.d0)negr = negr + 1

                                if(EI.ge.EImax)then
                                   EImax = EI
                                   xe(:) = xin(:)
                                end if

                                if(yhat.le.Ymin)then
                                   Ymin  = yhat
                                   xy(:) = xin(:)
                                end if

                                if(Cmode(:7).eq.'Post_1D')then
                                   TEC1(i,1) = xin(1)*(DS(2,1)-DS(1,1))+DS(1,1)
                                   TEC1(i,2) = yhat
                                   TEC1(i,3) = RMSE
                                   TEC1(i,4) = EI
                                   if (fct.lt.20) then
                                      call evalfunc(xin,ndim,fct,0,0,f,df,d2f,v)
                                      TEC1(i,5) = f
                                      diff = diff + (f-yhat)**2
                                      if (abs(f-yhat).gt.maxerror) maxerror=abs(f-yhat)
                                   end if
                                   go to 110

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                else if(Cmode(:7).eq.'Post_2D')then
                                   !  if (readfromdatabase.eq.0) then                        
                                   TEC2(i,j,1) = xin(1)*(DS(2,1)-DS(1,1))+DS(1,1)
                                   TEC2(i,j,2) = xin(2)*(DS(2,2)-DS(1,2))+DS(1,2)
                                   ! end if
                                   TEC2(i,j,3) = yhat
                                   TEC2(i,j,4) = RMSE
                                   TEC2(i,j,5) = EI
                                   ! if (readfromdatabase.eq.0) then
                                   if (fct.lt.20) then 
                                      call evalfunc(xin,ndim,fct,0,0,f,df,d2f,v)
                                      TEC2(i,j,6) = f
                                      diff = diff + (f-yhat)**2
                                      if (abs(f-yhat).gt.maxerror) maxerror=abs(f-yhat)
                                   end if
                                   if (fct.ge.20) then
                                      TEC2(i,j,1) = TEC2(i,j,1) * 180.0 / 4.0 /atan(1.0)   !RADINS TO degree
                                      f=TEC2(i,j,6)
                                      diff = diff + (f-yhat)**2
                                      if (abs(f-yhat).gt.maxerror) maxerror=abs(f-yhat)
                                   end if
                                   !end if
                                end if

                                TEC2(i,j,6) = f
                                diff = diff + (f-yhat)**2
                                if (abs(f-yhat).gt.maxerror) maxerror=abs(f-yhat)

                                !  if (readfromdatabase.eq.0) then
                                !        if (k.eq.jmax-1) write(*,*) 'EVALFUNC is called',imax,'*',jmax,'times'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


120                             continue ! grid(j)
110                             continue ! grid(i)


                             end if




                             write(*,'(6x,a,i8,a)')'>> Negative RMSE on ',negr,' Pts'
                             write(*,'(6x,a,f12.7,a,99f8.3)') '>> EI   max = ',EImax,' @',(xe(i),i=1,ndim)
                             write(*,'(6x,a,f12.7,a,99f8.3)') '>> Yhat min = ', Ymin,' @',(xy(i),i=1,ndim)

                             if(ndim.eq.2)then
                                diff = dsqrt(diff/dble(imax*jmax))
                                write(*,*)
                                write(*,'(6x,a,e20.10)') '>> RMSE on 2D Function = ',diff
                                write(*,*)
                                write(93,'(i4,4e15.8)') nhs,diff,maxerror,difflocmax,diffloc2
                                rmsemat(runnum,loopcounter,1)=nhs
                                rmsemat(runnum,loopcounter,2)=diff
                             end if

!!$            xin(1)=0.25
!!$            xin(2)=0.25
!!$            call meta_call(1,2,xin,yhat,yhatprime,RMSE,EI)                
!!$            call evalfunc(xin,ndim,fct,0,1,f,df,d2f,v)
!!$            print *,'Derivative at',xin
!!$            print *,yhatprime
!!$            print *,df
!!$            print *




!!$         ict = 0
!!$         do i=1,imax
!!$            do j=1,jmax
!!$               ict = ict + 1
!!$               xin(1) = dble(i-1)/dble(imax-1)
!!$               if(ndim.eq.2) xin(2) = dble(j-1)/dble(jmax-1)
!!$               XSAMP(:,ict)=DS(1,:)+(DS(2,:)-DS(1,:))*xin(:)
!!$            end do
!!$         end do
!!$
!!$
!!$         ict = 0
!!$         do l=1,rsample
!!$            if(inf(l)(1:2).eq.'1_')then 
!!$               ict = ict + 1
!!$               xvis(:,ict) = DS(1,:)+(DS(2,:)-DS(1,:))*sampl(l,:)
!!$               fvis(ict,1)   = fun(l)                   
!!$               if(inf(l)(4:4).eq.'G')then 
!!$                  fgvis(:,ict,1)= gfun(l,:)/(DS(2,:)-DS(1,:))
!!$               end if
!!$            end if
!!$         end do
!!$
!!$         Taylororder=ict
!!$         
!!$         ! Calculate the best parameters beta and gamma
!!$         CALL MIR_BETA_GAMMA(nfunc-1, ndim, ict, xvis, fvis, SIGV, ict, xvis, fgvis, SIGG, Taylororder, 1, 1.0, BETA, GAMM, IERR)
!!$
!!$         WRITE (*,"('    BETA = ',F4.2,', GAMMA = 'F4.2)") BETA, GAMM
!!$         
!!$         ! Calculate Multivariate Interpolation and Regression approximation.
!!$         CALL MIR_EVALUATE(nfunc-1, ndim, imax*jmax, XSAMP, ict, xvis, fvis, SIGV, ict, xvis, fgvis, SIGG, BETA, GAMM, Taylororder, 1, FX, SIGMA, IERR)
!!$
!!$
!!$         !print *, i,j,f,yhat,FX(1,1),SIGMA(1)     
!!$         !stop
!!$
!!$         ict = 0
!!$         do i=1,imax
!!$            do j=1,jmax
!!$               ict = ict + 1               
!!$               TEC2(i,j,3) = FX(ict,1)
!!$            end do
!!$         end do















                          end if




                          ! just for output
                          call reduce_data(k,igrad,0.d0,0.d0)
                          if(igrad.eq.2)call indirect(k)

                          if(id_proc.eq.0)then

                             if(ndim.eq.1)then
                                write(*,'(6x,a)')'>> Output to gnusample/kriging.dat'
                                write(Csamp, '(a,i1,a)')'gnusamp0',k,'.dat'
                                write(Ckrig, '(a,i1,a)')'gnukrig0',k,'.dat'
                                write(Csampl,'(a,i1,a)')'gnusampl0',k,'.dat'

                                open(10,file=Csamp, form='formatted',status='unknown')
                                if(lmax.ne.1) open(11,file=Csampl,form='formatted',status='unknown')

                                do i=1,rsample
                                   if(inf(i)(1:2).eq.'1_')then
                                      write(10,'(2f15.8)')((DS(1,1)+(DS(2,1)-DS(1,1))*sampl(i,j)),j=1,ndim),fun(i)
                                   else
                                      write(11,'(2f15.8)')((DS(1,1)+(DS(2,1)-DS(1,1))*sampl(i,j)),j=1,ndim),fun(i)
                                   end if
                                end do
                                close(10)

                                if(lmax.ne.1) close(11)

                                open(10,file=Ckrig,form='formatted',status='unknown')
                                do i=1,imax
                                   write(10,'(5f15.8)')(TEC1(i,j),j=1,4)
                                end do
                                close(10)

                             else if(ndim.eq.2)then

                                write(*,'(6x,a)')'>> Output to tecsample/kriging.dat'
                                write(Csamp, '(a,i2.2,a)')'tecsamp',fct,'.dat'
                                write(Ckrig, '(a,i2.2,a)')'teckrig',fct,'.dat'
                                write(Csampl,'(a,i2.2,a)')'tecsampl',fct,'.dat'
                                write(Cexac, '(a,i2.2,a)')'tecex',fct,'.dat'

                                ict = 0
                                do i=1,rsample
                                   if(inf(i)(1:2).eq.'1_')then 
                                      ict = ict + 1
                                      xvis(:,ict) = DS(1,:)+(DS(2,:)-DS(1,:))*sampl(i,:)
                                      if (fct.ge.10) xvis(1,ict)=xvis(1,ict) * 180.0 / 4.0 /atan(1.0)   ! in degree
                                      fvis(ict,1)   = fun(i)
                                   end if
                                end do

                                open(10,file='./tec/'//Csamp, form='formatted',status='unknown')
                                write(10,'(a)')'TITLE = " "'
                                write(10,'(a)')'VARIABLES = "x" "y" "f"'
                                write(10,'(3(a,i5),a)')'ZONE T="hsample", I=',ict,', J=',1,', K=',1,', F=BLOCK'
                                write(10,'(9999f15.8)')(xvis(1,i),i=1,ict)
                                write(10,'(9999f15.8)')(xvis(2,i),i=1,ict)
                                write(10,'(9999f15.8)')(fvis(i,1),i=1,ict)
                                close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! added for exporting the samples to MIR
                                !       write(export, '(a,i3.3,a)')'KrigSamp',nhs,'.dat'
                                !      open(33,file='./KrigSamples/'//export, form='formatted',status='unknown')
                                !      write(33,'(a)')'TITLE = " "'
                                !      write(33,'(a)')'VARIABLES = "x" "y" "f"'
                                !      write(33,'(3(a,i5),a)')'ZONE T="hsample", I=',ict,', J=',1,', K=',1,', F=BLOCK'
                                !      write(33,'(9999f15.8)')(xvis(1,i),i=1,ict)
                                !      write(33,'(9999f15.8)')(xvis(2,i),i=1,ict)
                                !      write(33,'(9999f15.8)')(fvis(i,1),i=1,ict)
                                !      close(10)
                                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                if(lmax.ne.1)then
                                   ict = 0
                                   do i=1,rsample
                                      if(inf(i)(1:2).ne.'1_')then 
                                         ict = ict + 1
                                         xvis(:,ict) = DS(1,:)+(DS(2,:)-DS(1,:))*sampl(i,:)
                                         if (fct.ge.10) xvis(1,ict)=xvis(1,ict) * 180.0 / 4.0 /atan(1.0)   ! in degree
                                         fvis(ict,1)   = fun(i)
                                      end if
                                   end do
                                   open(11,file='./tec/'//Csampl, form='formatted',status='unknown')
                                   write(11,'(a)')'TITLE = " "'
                                   write(11,'(a)')'VARIABLES = "x" "y" "f"'
                                   write(11,'(3(a,i5),a)')'ZONE T="lsample", I=',ict,', J=',1,', K=',1,', F=BLOCK'
                                   write(11,'(9999f15.8)')(xvis(1,i),i=1,ict)
                                   write(11,'(9999f15.8)')(xvis(2,i),i=1,ict)
                                   write(11,'(9999f15.8)')(fvis(i,1),i=1,ict)
                                   close(11)
                                else
                                   open(11,file='./tec/'//Csampl, form='formatted',status='unknown')
                                   write(11,'(a)')'TITLE = " "'
                                   write(11,'(a)')'VARIABLES = "x" "y" "f"'
                                   write(11,'(3(a,i5),a)')'ZONE T="lsample", I=',1, ', J=',1,', K=',1,', F=BLOCK'
                                   write(11,'(9999f15.8)') xvis(1,1)
                                   write(11,'(9999f15.8)') xvis(2,1)
                                   write(11,'(9999f15.8)') fvis(1,1)
                                   close(11)
                                end if

                                open(10,file='./tec/'//Ckrig,form='formatted',status='unknown')
                                write(10,'(a)')'TITLE = " "'
                                write(10,'(a)')'VARIABLES = "x" "y" "f" "rmse" "EI"'
                                write(10,'(3(a,i5),a)')'ZONE T="Kriging", I=',imax,', J=',jmax,', K=',1,', F=BLOCK'
                                write(10,*)((TEC2(i,j,1),i=1,imax),j=1,jmax)
                                write(10,*)((TEC2(i,j,2),i=1,imax),j=1,jmax)
                                write(10,*)((TEC2(i,j,3),i=1,imax),j=1,jmax)
                                write(10,*)((TEC2(i,j,4),i=1,imax),j=1,jmax)
                                write(10,*)((TEC2(i,j,5),i=1,imax),j=1,jmax)
                                close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                !!Comment  this  cell out if reading samples from tecex10.dat!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                !!MAke sure you read the exact values above from this file instead of callinf exactfunction calls !!!!!!!!
                                if (readfromdatabase.eq.0) then
                                   ! write(*,*)'writing exact function values to tecex**.dat'
                                   if (fct.le.10) then
                                      open(10,file='./tec/'//Cexac,form='formatted',status='unknown')
                                      write(10,'(a)')'TITLE = " "'
                                      write(10,'(a)')'VARIABLES = "x" "y" "f"'
                                      write(10,'(3(a,i5),a)')'ZONE T="Exact", I=',imax,', J=',jmax,', K=',1,', F=BLOCK'
                                      write(10,*)((TEC2(i,j,1),i=1,imax),j=1,jmax)
                                      write(10,*)((TEC2(i,j,2),i=1,imax),j=1,jmax)
                                      write(10,*)((TEC2(i,j,6),i=1,imax),j=1,jmax)
                                      close(10)
                                   end if
                                end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$                open(10,file='rmse2D.dat',status='unknown',form='formatted')
!!$                write(10,'(2i8,e20.10)')rsample-psample,psample,diff
!!$                close(10)
                             end if

                          end if


                       end do  ! loop over functions

                       if(ndim.eq.1)then
                          deallocate( TEC1 )
                       else if(ndim.eq.2)then
                          deallocate( TEC2 )
                       end if


                     end subroutine Post_1or2D

