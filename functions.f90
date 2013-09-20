  subroutine evalfunc(x,DIM,fct,ifid,flag,f,df,d2f,v)
    use dimKrig, only: DS,fctindx,reusesamples,ndimt,xavgt,xstdt
    use omp_lib
    implicit none

    integer :: DIM,ifid,fct,flag,j,k,nseed(1)
    real*8 :: x(DIM),f,df(DIM),d2f(DIM,DIM),scal,scal2,prd,time,time2
    real*8 :: xtmp(ndimt),v(ndimt),dftmp(ndimt),d2ftmp(ndimt,ndimt)
    real*8 :: gtol,low(ndimt-DIM),up(ndimt-DIM)
    character(len=6)  :: Cstat
    character*2 :: fctindxnumber
    character*60 :: filename

    ! ifid      0: high fidelity  other:lower fidelity

    xtmp(1:ndimt-DIM)=xavgt(1:ndimt-DIM)
    do k=1,DIM
       scal=DS(2,k)-DS(1,k)
       xtmp(ndimt-DIM+k)=x(k)*scal+DS(1,k)
    end do

    if (fct.lt.10) then

       call calcf(xtmp,ndimt,fct,f)
       if (flag.ge.1) call calcdf(xtmp,ndimt,fct,dftmp)   
       if (flag.ge.2) call calcd2f(xtmp,ndimt,fct,d2ftmp)
          
       if (ifid.gt.0) then
          f=0.1*f
          dftmp(:)=0.1*dftmp(:)
          d2ftmp(:,:)=0.1*d2ftmp(:,:)
       end if

    else if (fct.eq.10) then
       CALL chdir('lowfid') ! Comment when using fine mesh

       call omp_set_num_threads(omp_get_max_threads())
       call Eulersolve(xtmp,ndimt,ifid,f,dftmp,d2ftmp,flag,v,fctindx)

       CALL chdir('../') !Comment when using fine mesh


    else if (fct.eq.11) then

       gtol=1e-6

       low(1:ndimt-DIM)=xtmp(1:ndimt-DIM)-xstdt(1:ndimt-DIM)
       up(1:ndimt-DIM)=xtmp(1:ndimt-DIM)+xstdt(1:ndimt-DIM)

       if (flag.ge.1) then
          write(*,*) 'Error in function call, optimization does not support gradient evalution'
          stop
       end if

       call omp_set_num_threads(omp_get_max_threads())
       call optimize(ndimt-DIM,xtmp,ndimt,f,df,low,up,gtol,.true.,.false.,fctindx+1)

    else if (fct.eq.12) then

       gtol=1e-6

       low(1:ndimt-DIM)=xtmp(1:ndimt-DIM)-xstdt(1:ndimt-DIM)
       up(1:ndimt-DIM)=xtmp(1:ndimt-DIM)+xstdt(1:ndimt-DIM)

       if (flag.ge.1) then
          write(*,*) 'Error in function call, optimization does not support gradient evalution'
          stop
       end if

       call omp_set_num_threads(omp_get_max_threads())

       if (fctindx.eq.0) then 

          call optimize(ndimt-DIM,xtmp,ndimt,f,dftmp,low,up,gtol,.true.,.false.,fctindx+10)

       else if (fctindx.eq.4) then

          call optimize(ndimt-DIM,xtmp,ndimt,f,dftmp,low,up,gtol,.false.,.false.,fctindx+10)

       else
          write(*,*) 'Wrong fctindx in function call'
          stop
       end if   

    end if

    df(1:DIM)=dftmp(ndimt-DIM+1:ndimt)
    d2f(1:DIM,1:DIM)=d2ftmp(ndimt-DIM+1:ndimt,ndimt-DIM+1:ndimt)
    do k=1,DIM
       scal=DS(2,k)-DS(1,k)
       df(k)=df(k)*scal
       do j=1,DIM
          scal2=DS(2,j)-DS(1,j)
          d2f(j,k)=d2f(j,k)*scal*scal2
       end do
    end do


    if (flag.eq.5) then   ! Hessian-vector product

       prd=0.0
       do k=1,DIM
          prd = prd + df(k)**2
       end do
       prd = dsqrt(prd)

       if (prd.gt.1.e-5) then ! Gradient not zero

          do k=1,DIM
             v(k) = df(k)/prd
          end do

       else
          
          ! Trick to get different random numbers for each call
          call get_seed(nseed)
          time2=secnds(time)
          time2=time2-int(time2)
          nseed(1)=nseed(1)+2*int(time2*1000000000.0)      
          call random_seed(put=nseed)
          call random_number(v)

          v(:)=2.0*v(:)-1.0

          prd=0.0 
          do k=1,DIM
              prd = prd + v(k)**2
          end do
          prd = dsqrt(prd)

          do k=1,DIM
             v(k) = v(k)/prd
          end do

       end if

       d2f(:,1) = matmul(d2f,v)

    end if

    if (reusesamples.eq.1) then

       ! Keep running tally of all function evaluations

       if (ifid.eq.0) then
          
          if(flag.eq.0)then
             Cstat = '1_F   '
          else if(flag.eq.1)then
             Cstat = '1_FG  '
          else if(flag.eq.3)then
             Cstat = '1_FGH '
          else if(flag.eq.5)then
             Cstat = '1_FGHv'
          end if
          
       end if
    
       if (ifid.gt.0) then
          
          ! Cstatusl
          if(     flag.eq.0)then
             Cstat = '2_F   '
          else if(flag.eq.1)then
             Cstat = '2_FG  '
          else if(flag.eq.3)then
             Cstat = '2_FGH '
          else if(flag.eq.5)then
             Cstat = '2_FGHv'
          end if
          
       end if

       filename='Function'
       call i_to_s(fctindx,fctindxnumber)
       filename(9:10)=fctindxnumber
       filename(11:14)='.dat'
       
       open(97,file=filename,form='formatted',status='unknown',access='append')
       if (flag.le.3) write(97,100) Cstat,(x(j)*(DS(2,j)-DS(1,j))+DS(1,j),j=1,dim),f,(df(j)/(DS(2,j)-DS(1,j)),j=1,dim),((d2f(k,j)/(DS(2,j)-DS(1,j))/(DS(2,k)-DS(1,k)),j=1,dim),k=1,dim) 
       if (flag.gt.3) write(97,100) Cstat,(x(j)*(DS(2,j)-DS(1,j))+DS(1,j),j=1,dim),f,(df(j)/(DS(2,j)-DS(1,j)),j=1,dim),(v(j),j=1,dim),(d2f(k,1),k=1,dim)
       close(97)

    end if

100    format(a,1x,10000e20.10)

  end subroutine evalfunc





  subroutine calcf(x,DIM,fct,f)
    use dimKrig, only: fctindx
    implicit none

    integer :: DIM,fct,k
    real*8 :: x(DIM),f,A,omeg

    real*8 :: rho, L, sigmay, pi, Fs, p, E
    
    if (fct.eq.1) then
       f=0.0
       do k=1,DIM
          f=f+x(k)
       end do
       f=cos(f)
    else if (fct.eq.2) then
       f=1.0
       do k=1,DIM
          f=f+x(k)**2
       end do
       f=1.0/f
    else if (fct.eq.3) then 
       f=0.0
       do k=1,DIM-1
          f=f+100.0*(x(k+1)-x(k)**2)**2+(1-x(k))**2
       end do
    else if (fct.eq.4) then

       A=10.0
       omeg=2.0*3.141592653589793

       f=A*real(DIM)
       do k=1,DIM
          f=f+x(k)**2-A*cos(omeg*x(k))
       end do
       f=f+5.0
    else if (fct.eq.5) then
       f=0.0
       do k=1,DIM
          f=f+x(k)
       end do
       f=cos(f)+0.01*cos(100.0*f)
    else if (fct.eq.6) then

       rho=0.2836
       sigmay=36260.0
       p=25000.0
       L=5.0
       E=30e6
       pi=4.0*atan(1.0)

       Fs=1.0

       if (fctindx.eq.0) then
          !---- OBJECTIVE FUNCTION
          f = rho*x(1)*L+rho*x(2)*sqrt(L**2+x(3)**2)
       else if (fctindx.eq.1) then
!---- INEQUALITY CONSTRAINTS
          f = p*Fs*sqrt(L**2+x(3)**2) / (x(2)*x(3)*sigmay) - 1.0
       else if (fctindx.eq.2) then  
          f = p*Fs*L / (x(1)*x(3)*sigmay) - 1.0
       else if (fctindx.eq.3) then
          f = 4.0*p*Fs*L**3 / (x(1)**2*x(3)*E*pi) - 1.0
       end if
       
    else if (fct.eq.7) then
       f=0.0
       do k=1,DIM
          f=f+x(k)**2
       end do
       f=f+5.0
     
    else if (fct.eq.8) then
       f=0.0
       do k=1,DIM
          f=f+x(k)**3
       end do
       f=f+5.0

    else if (fct.eq.9) then
  
!!$       if (x(1).le.1.0 .and. x(1).ge.-1.0 .and. x(2).le.1.0 .and. x(2).ge.-1.0) then
!!$          f=2.0
!!$       else
!!$          f=1.0
!!$       end if

       f=0.75*exp(-0.25*((9.0*x(1)-2.0)**2+(9.0*x(2)-2.0)**2))  &
        +0.75*exp(-1.0/49.0*(9.0*x(1)+1.0)**2-0.1*(9.0*x(2)+1.0)**2)  &
        +0.5*exp(-0.25*((9.0*x(1)-7.0)**2+(9.0*x(2)-3.0)**2))  &
        -0.2*exp(-(9.0*x(1)-4.0)**2-(9.0*x(2)-7.0)**2)

    end if

  end subroutine calcf

  
  subroutine calcdf(x,DIM,fct,df)
    use dimKrig, only: fctindx
    implicit none
    integer :: DIM,fct,k
    real*8 :: x(DIM),df(DIM),fac,A,omeg

    real*8 :: rho, L, sigmay, pi, Fs, p, E

    
    if (fct.eq.1) then
       !f=cos(x+y)

       fac=0.0
       do k=1,DIM
          fac=fac+x(k)
       end do
       
       fac=-sin(fac)
       do k=1,DIM
          df(k)=fac
       end do
       
    else if (fct.eq.2) then  
       !f=1.0/(1.0+x**2+y**2)

       fac=1.0
       do k=1,DIM
          fac=fac+x(k)**2
       end do
       fac=1.0/fac**2
       
       do k=1,DIM
          df(k)=-2.0*x(k)*fac
       end do

    else if (fct.eq.3) then 
       !f=(1.0-x)**2 + 100.0*(y-x**2)**2

  
       df(1)=-200.0*(x(2)-x(1)**2)*2.0*x(1)-2.d0*(1-x(1))
       do k=2,DIM-1
          df(k)=200.0*(x(k)-x(k-1)**2)-200.0*(x(k+1)-x(k)**2)*2.0*x(k)-2.0*(1-x(k))
       end do
       df(DIM)=200.0*(x(DIM)-x(DIM-1)**2)

    else if (fct.eq.4) then
       ! Rastrigin

       A=10.0
       omeg=2.0*3.141592653589793

       do k=1,DIM
          df(k)=2.0*x(k)+A*omeg*sin(omeg*x(k))
       end do 

    else if (fct.eq.5) then
       !f=cos(x+y)+0.01*cos(100*(x+y))

       fac=0.0
       do k=1,DIM
          fac=fac+x(k)
       end do
       
       fac=-sin(fac)-1.0*sin(100.0*fac)
       do k=1,DIM
          df(k)=fac
       end do
  
    else if (fct.eq.6) then 
 
       rho=0.2836
       sigmay=36260.0
       p=25000.0
       L=5.0
       E=30e6
       pi=4.0*atan(1.0)

       Fs=1.0

       if (fctindx.eq.0) then
          !---- OBJECTIVE FUNCTION
          df(1) = rho*L
          df(2) = rho*sqrt(L**2+x(3)**2)
          df(3) = rho*x(2)*x(3) / sqrt(L**2+x(3)**2)
       else if (fctindx.eq.1) then
!---- INEQUALITY CONSTRAINTS
          df(1) = 0.0
          df(2) =-p*Fs*sqrt(L**2+x(3)**2) / (x(2)**2*x(3)*sigmay)
          df(3) =-p*Fs*L**2 /sqrt(L**2/x(3)**2+1.0)/ (x(2)*x(3)**3*sigmay)
       else if (fctindx.eq.2) then  
          df(1) =-p*Fs*L / (x(1)**2*x(3)*sigmay)
          df(2) = 0.0
          df(3) =-p*Fs*L / (x(1)*x(3)**2*sigmay)
       else if (fctindx.eq.3) then
          df(1) =-8.0*p*Fs*L**3 / (pi*E*x(1)**3*x(3))
          df(2) = 0.0
          df(3) =-4.0*p*Fs*L**3 / (pi*E*x(1)**2*x(3)**2)
       end if

    else if (fct.eq.7) then  
       !f=x**2+y**2
       
       do k=1,DIM
          df(k)=2.0*x(k)
       end do
   
    else if (fct.eq.8) then
       !f=x**3+y**3

       do k=1,DIM
          df(k)=3.0*x(k)**2
       end do

    else if (fct.eq.9) then

       !df(:)=0.0

       df(1)=0.75*exp(-0.25*((9.0*x(1)-2.0)**2+(9.0*x(2)-2.0)**2))  *2.0*(-0.25)*(9.0*x(1)-2.0)*9.0  &
        +0.75*exp(-1.0/49.0*(9.0*x(1)+1.0)**2-0.1*(9.0*x(2)+1.0)**2) * 2.0*(-1.0/49.0)*(9.0*x(1)+1.0)*9.0 &
        +0.5*exp(-0.25*((9.0*x(1)-7.0)**2+(9.0*x(2)-3.0)**2)) *2.0* (-0.25)*(9.0*x(1)-7.0)*9.0 &
        -0.2*exp(-(9.0*x(1)-4.0)**2-(9.0*x(2)-7.0)**2) *2.0*(9.0*x(1)-4.0)*(-9.0)

       df(2)=0.75*exp(-0.25*((9.0*x(1)-2.0)**2+(9.0*x(2)-2.0)**2)) *2.0*(-0.25)*(9.0*x(2)-2.0)*9.0  &
        +0.75*exp(-1.0/49.0*(9.0*x(1)+1.0)**2-0.1*(9.0*x(2)+1.0)**2) *2.0*(-0.1)*(9.0*x(2)+1.0)*9.0 &
        +0.5*exp(-0.25*((9.0*x(1)-7.0)**2+(9.0*x(2)-3.0)**2)) *2.0*(-0.25)*(9.0*x(2)-3.0)*9.0  &
        -0.2*exp(-(9.0*x(1)-4.0)**2-(9.0*x(2)-7.0)**2) *2.0*(9.0*x(2)-7.0)*(-9.0)
       
    end if
      
  end subroutine calcdf

  subroutine calcd2f(x,DIM,fct,d2f)
    implicit none
    integer :: DIM,fct,j,k
    real*8 :: x(DIM),d2f(DIM,DIM),fac,A,omeg


    if (fct.eq.1) then
       
       !f=cos(x+y)

       fac=0.0
       do k=1,DIM
          fac=fac+x(k)
       end do
       
       fac=-cos(fac)
       do j=1,DIM
          do k=1,DIM
             d2f(j,k)=fac
          end do
       end do
       
    else if (fct.eq.2) then
       
       !f=1.0/(1.0+x**2+y**2)

       fac=1.0
       do k=1,DIM
          fac=fac+x(k)**2
       end do
       fac=1.0/fac
       
       do j=1,DIM
          do k=1,DIM
             if (j.ne.k) then
                d2f(j,k)=8.0*x(j)*x(k)*fac**3
             else
                d2f(j,k)=(-2.0+8.0*x(k)**2*fac)*fac**2
             end if
          end do
       end do
  
    else if (fct.eq.3) then 
       !f=(1.0-x)**2 + 100.0*(y-x**2)**2

       d2f(:,:)=0.0

       d2f(1,1)=(-400.0*x(2)+1200.0*x(1)**2+2.0)
       do k=2,DIM-1
          d2f(k,k)=200.0-400.0*x(k+1)+1200.0*x(k)**2+2.0
       end do
       d2f(DIM,DIM)=200.0
     
       do k=2,DIM
          d2f(k,k-1)=-400.0*x(k-1)
          d2f(k-1,k)=d2f(k,k-1)
       end do

    else if (fct.eq.4) then
       ! Rastrigin

       A=10.0
       omeg=2.0*3.141592653589793

       d2f(:,:)=0.0
       do k=1,DIM
          d2f(k,k)=2.0+A*omeg**2*cos(omeg*x(k))
       end do  

    else if (fct.eq.5) then

       !f=cos(x+y)+0.01*cos(100*(x+y))

       fac=0.0
       do k=1,DIM
          fac=fac+x(k)
       end do
       
       fac=-cos(fac)-100.0*cos(100.0*fac)
       do j=1,DIM
          do k=1,DIM
             d2f(j,k)=fac
          end do
       end do

    else if (fct.eq.6) then
       
       d2f(:,:)=0.0

    else if (fct.eq.7) then

       !f=x**2+y**2

       d2f(:,:)=0.0
       do k=1,DIM
          d2f(k,k)=2.0
       end do 
      
 
    else if (fct.eq.8) then

       !f=x**3+y**3
       
       d2f(:,:)=0.0
       do k=1,DIM
          d2f(k,k)=6.0*x(k)
       end do 

    else if (fct.eq.9) then

       d2f(:,:)=0.0
                
    end if
      
  end subroutine calcd2f




