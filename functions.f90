
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

  
    xtmp(1:ndimt)=xavgt(1:ndimt)

 
    do k=1,DIM
       scal=DS(2,k)-DS(1,k)
       xtmp(ndimt-DIM+k)=x(k)*scal+DS(1,k)
    end do


  if (fct.lt.20) then

     call calcf(xtmp,ndimt,fct,f)
     if (flag.ge.1) call calcdf(xtmp,ndimt,fct,dftmp)   
     if (flag.ge.2) call calcd2f(xtmp,ndimt,fct,d2ftmp)

     if (ifid.gt.0) then
        f=0.1*f
        dftmp(:)=0.1*dftmp(:)
        d2ftmp(:,:)=0.1*d2ftmp(:,:)
     end if

  else if (fct.eq.20) then

     CALL chdir('lowfid') ! Comment when using fine mesh

     call omp_set_num_threads(omp_get_max_threads())
     call Eulersolve(xtmp,ndimt,ifid,f,dftmp,d2ftmp,flag,v,fctindx)

     CALL chdir('../') !Comment when using fine mesh

     ! Meant for optimizations

  else if (fct.eq.21) then ! Two bar truss 
     
     gtol=1e-6

     low(1:ndimt-DIM)=xtmp(1:ndimt-DIM)-xstdt(1:ndimt-DIM)
     up(1:ndimt-DIM)=xtmp(1:ndimt-DIM)+xstdt(1:ndimt-DIM)

     if (flag.ge.1) then
        write(*,*) 'Error in function call, optimization does not support gradient evalution'
        stop
     end if


  if (fctindx.eq.0) then ! what is the worst objective function 

         call optimize(ndimt-DIM,xtmp,ndimt,f,dftmp,low,up,gtol,.false.,.false.,fctindx)

      else

        call optimize(ndimt-DIM,xtmp,ndimt,f,dftmp,low,up,gtol,.true.,.false.,fctindx)
        
      end if

!     call optimize(ndimt-DIM,xtmp,ndimt,f,dftmp,low,up,gtol,.true.,.false.,fctindx)

  else if (fct.eq.22) then !CFD

     gtol=1e-6

     low(1:ndimt-DIM)=xtmp(1:ndimt-DIM)-xstdt(1:ndimt-DIM)
     up(1:ndimt-DIM)=xtmp(1:ndimt-DIM)+xstdt(1:ndimt-DIM)

     if (flag.ge.1) then
        write(*,*) 'Error in function call, optimization does not support gradient evalution'
        stop
     end if

!     print*,'before optimize',fctindx,fct,xtmp,ndimt,low(1:ndimt-dim),up(1:ndimt-dim)

     call omp_set_num_threads(omp_get_max_threads())

     if (fctindx.eq.0) then !what is the max possible drag? (objective function)

        call optimize(ndimt-DIM,xtmp,ndimt,f,dftmp,low,up,gtol,.false.,.false.,fctindx)

     else if (fctindx.eq.4) then !what is the least lift possible? (lift constraint)

      call optimize(ndimt-DIM,xtmp,ndimt,f,dftmp,low,up,gtol,.true.,.false.,fctindx)


     else
        write(*,*) 'Wrong fctindx in function call'
        stop
     end if

!   print*,'after optimize',xtmp,f,dftmp

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

100 format(a,1x,10000e20.10)

end subroutine evalfunc





subroutine calcf(x,DIM,fct,f)
  use dimKrig, only: fctindx, DAT,mainprog
  implicit none

  integer :: DIM,fct,k
  real*8 :: x(DIM),f,A,omeg
  real*8:: tau_allow,M,V,B,D
  real*8 :: rho, L, sigmay, pi, Fs, p, E,R, T,sigma_allow !Truss design parameters
  real*8::gamma
  real*8::tensile_sigma1_max,tensile_sigma2_max,tensile_sigma3_max
  real*8::comp_sigma1_max,comp_sigma2_max,comp_sigma3_max
  real*8::max_u_disp,max_v_disp,theta,pu,pv,u,sigma(dim)

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
!     if (dim.ne.3) stop'Wrong dimension'

     pi=4.0*atan(1.0)

!!$     ! Storing variables into DAT array to passthrough
!!$     ! Use this in your calling subprogam
!!$     
!!$     dat(1)=rho
!!$     dat(2)=L
!!$     dat(3)=E
!!$     dat(4)=P
!!$     dat(5)=sigmay
!!$     dat(6)=Fs

     if (mainprog) then
        rho=0.2836
        sigmay=36260.0
        p=25000.0
        L=5.0
        E=30.0e6
        Fs=1.0
     else

        !Read variables from DAT array to use in expressions

        rho=dat(1)
        L=dat(2)
        E=dat(3)
        P=dat(4)
        sigmay=dat(5)
        Fs=dat(6)

     end if

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
     else
        stop'Wrong fct index'
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
  else if (fct.eq.9) then ! Tubular column

!     if (dim.ne.3) stop'Wrong dimension'

     !Thanks: Arora Section 3.7

     pi=4.0*atan(1.0)

     if (mainprog) then

        rho=7.833e-6           !kg/m3
        p=10.0e6               !10 MN
        E=207000.0               !N/mm2
        sigma_allow=248.0        !N/mm2
        Fs=1.0

     else

        !Read variables from DAT array to use in expressions

        rho=dat(1)
        E=dat(2)
        P=dat(3)
        sigma_allow=dat(4)
        Fs=dat(5)

     end if

!!$
!!$     ! Store variables into DAT array to passthrough
!!$     ! Use this in your calling program
!!$     
!!$       dat(1)=rho
!!$       dat(2)=E
!!$       dat(3)=P
!!$       dat(4)=sigma_allow
!!$       dat(5)=Fs


     !define R,T

     R=x(1)
     T=x(2)
     L=x(3)

     if (fctindx.eq.0) then

        !---- OBJECTIVE FUNCTION
        f = 2.0*rho*L*pi*R*T

     else if (fctindx.eq.1) then

        !---- INEQUALITY CONSTRAINTS
        f = P*Fs / (2.0*pi*R*T*sigma_allow) - 1.0

     else if (fctindx.eq.2) then

        f = 4.0*p*Fs*L**2 / (t*E*(R*pi)**3) - 1.0

        !     else if (fctindx.eq.3) then

        !        f= x(1)*Fs/(32.0*x(2)) - 1.0

     else

        print*, 'Wrong function index for this test case',fctindx
        stop
     end if


  else if (fct.eq.10) then ! Cantilever beam 
!     if (dim.ne.2) stop'Wrong dimension'

     ! Thanks: Section 3.8 Arora 
     if (mainprog) then

        M=40.0d6 !Nm
        V=150000.0 !N
        sigma_allow= 10.0d6 !N/m2
        tau_allow= 2.0d6 !N/m2     
        Fs=1.0

     else

        M=dat(1)
        V=dat(2)
        sigma_allow=dat(3)
        tau_allow=dat(4)
        Fs=dat(5)

     end if

!!$       !Store variables into DAT array to passthrough
!!$       ! Use this in your calling program
!!$       dat(1)=M
!!$       dat(2)=V
!!$       dat(3)=sigma_allow
!!$       dat(4)=tau_allow
!!$       dat(5)=Fs

     !define R,T

     B=x(1)
     D=x(2)
     !     L=x(3)

     if (fctindx.eq.0) then

        !---- OBJECTIVE FUNCTION
        f = B*D

     else if (fctindx.eq.1) then

        !---- INEQUALITY CONSTRAINTS
        !bending stress constraint

        f=(6.0*M*fs)/(b*(d**2)*sigma_allow)-1.0       

     else if (fctindx.eq.2) then

        ! Inequality constraint 2
        ! Shear Stress constraint

        f=(3.0*V*fs)/(2.0*b*d*tau_allow) -1.0

     else if (fctindx.eq.3) then

        f= d*FS/(2.0*b) - 1.0

     else 

        print*, 'Wrong function index for this test case',fctindx
        stop
     end if

  else if (fct.eq.11) then ! Three bar truss 


!     if (dim.ne.3) stop'wrong dim for this problem'

     pi=4.0*atan(1.0)

     if (mainprog) then

        !Problem data and other constants
        L=10.0 !height ref
        E=1.0e7 !E
        gamma=0.1 !gamma
        theta=135.0*pi/180.0
        P=20000.0

        ! Max constraint values

        !Tensile
        tensile_sigma1_max=5000.0    ! psi tensile_sigma1_max=dat(6)      
        tensile_sigma2_max=20000.0    ! psi tensile_sigma2_max=dat(7)
        tensile_sigma3_max=5000.0    ! psi tensile_sigma3_max=dat(8)
        !Compressive
        comp_sigma1_max=5000.0    ! psi comp_sigma1_max=dat(9)
        comp_sigma2_max=20000.0   ! psi comp_sigma2_max=dat(10)
        comp_sigma3_max=5000.0   ! psi comp_sigma3_max=dat(11)
        !Displacement
        max_u_disp=0.005    ! in  max_u_disp=dat(12)
        max_v_disp=0.005    ! in  max_v_disp=dat(12)

        Fs=1.0      ! Factor of safety

     else

        L=dat(1)
        E=dat(2)
        gamma=dat(3)
        theta=dat(4)
        P=dat(5)

        ! Max constraint values

        !Tensile
        tensile_sigma1_max=dat(6)
        tensile_sigma2_max= dat(7)
        tensile_sigma3_max= dat(8)

        !Compressive
        comp_sigma1_max= dat(9)
        comp_sigma2_max= dat(10)
        comp_sigma3_max=dat(11)

        !Displacement
        max_u_disp= dat(12)
        max_v_disp= dat(13)
        Fs= dat(14)
     end if

!!$     ! Use this in caling program and pass to pcestimate subtoutine
!!$     !Problem data and other constants
!!$     dat(1)=10.0 !height ref
!!$     dat(2)=1.0e7 !E
!!$     dat(3)=0.1 !gamma
!!$     dat(4)=45.0*pi/180.0
!!$     dat(5)=20000.0
!!$
!!$     ! Max constraint values
!!$
!!$     !Tensile
!!$     dat(6)=5000.0    ! psi tensile_sigma1_max=dat(6)      
!!$     dat(7)=20000.0    ! psi tensile_sigma2_max=dat(7)
!!$     dat(8)=5000.0    ! psi tensile_sigma3_max=dat(8)
!!$     !Compressive
!!$     dat(9)=5000.0    ! psi comp_sigma1_max=dat(9)
!!$     dat(10)=20000.0   ! psi comp_sigma2_max=dat(10)
!!$     dat(11)=5000.0   ! psi comp_sigma3_max=dat(11)
!!$     !Displacement
!!$     dat(12)=0.005    ! in  max_u_disp=dat(12)
!!$     dat(13)=0.005    ! in  max_v_disp=dat(12)
!!$     dat(14)=1.0      ! Factor of safety
!!$

     pu=P*cos(theta)
     pv=P*sin(theta)

     if (fctindx.gt.0) then
        u=(L/E)*(x(1)*pu + 2*dsqrt(2.0)*x(2)*pu + x(3)*pu + x(3)*pv - x(1)*pv)/(x(1)*x(2) + dsqrt(2.0)*x(1)*x(3) + x(2)*x(3))

        v=(L/E)*(-x(1)*pu + x(3)*pu + x(1)*pv + x(3)*pv)/(x(1)*x(2) + dsqrt(2.0)*x(1)*x(3) + x(2)*x(3))

!!$if (loadcase.eq.2) then

        sigma(1)=-1.0*(sqrt(2.0)*x(2)*pu + x(3)*pu  +x(3)*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))

!!$else

!!$   sigma(1)=(sqrt(2.0)*x(2)*pu + x(3)*pu  +x(3)*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))

!!$ end if

        sigma(2)= (-(x(1)-x(3))*pu+(x(1)+x(3))*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))

        sigma(3)=(-sqrt(2.0)*x(2)*pu -x(1)*pu  +x(1)*pv)/(x(1)*x(2)+sqrt(2.0)*x(1)*x(3)+x(2)*x(3))

     end if

     if (fctindx.eq.0) then

        F= x(1)*gamma*L*sqrt(2.0) + x(2)*gamma*L +  x(3)*gamma*L*sqrt(2.0)  

     else if (fctindx.eq.1) then

        F = (sigma(1) - tensile_sigma1_max)/tensile_sigma1_max    !tensile 1

     else if (fctindx.eq.2) then

        F = (sigma(2) - tensile_sigma2_max)/tensile_sigma2_max   !tensile 2

     else if (fctindx.eq.3) then

        F = (sigma(3) - tensile_sigma3_max)/tensile_sigma3_max    ! tensile 3

     else if (fctindx.eq.4) then

        F = (-1.0*sigma(1)/comp_sigma1_max) -1.0     !compressive 1

     else if (fctindx.eq.5) then

        F = (-1.0*sigma(2)/comp_sigma2_max) -1.0  !compressive 2

     else if (fctindx.eq.6) then

        F = (-1.0*sigma(3) / comp_sigma3_max) -1.0 !compressive 3

     else if (fctindx.eq.7) then

        F = (-1.0*u -max_u_disp)/max_u_disp

     else if (fctindx.eq.8) then

        F = (v -max_v_disp)/max_v_disp

     else 

        print*,fctindx
        stop'Unsupported function index'

     end if

  else if (fct.eq.12) then ! Threebar 6d problem

!     if (dim.ne.6) stop'wrong dim for this problem'

     pi=4.0*atan(1.0)

     if (mainprog) then
        ! Use these settings if the program is called from main.f90.
        ! If PC is used as library DAT is passed as an input vector from calling program such as IPOPT

        !Problem data and other constants
        dat(1)=10.0 !height ref
        dat(2)=1.0e7 !E
        dat(3)=0.1 !gamma
        dat(4)=45.0*pi/180.0
        dat(5)=20000.0

        ! Max constraint values

        !Tensile
        dat(6)=5000.0    ! psi tensile_sigma1_max=dat(6)      
        dat(7)=20000.0    ! psi tensile_sigma2_max=dat(7)
        dat(8)=5000.0    ! psi tensile_sigma3_max=dat(8)
        !Compressive
        dat(9)=5000.0    ! psi comp_sigma1_max=dat(9)
        dat(10)=20000.0   ! psi comp_sigma2_max=dat(10)
        dat(11)=5000.0   ! psi comp_sigma3_max=dat(11)
        !Displacement
        dat(12)=0.005    ! in  max_u_disp=dat(12)
        dat(13)=0.005    ! in  max_v_disp=dat(12)
        dat(14)=1.0      ! Factor of safety
     end if

     if (fctindx.eq.0) then

        call threebarf(0,dat,x,F)

     else if (fctindx.eq.1) then

        call threebarf(1,dat,x,F)

     else if (fctindx.eq.2) then

        call threebarf(2,dat,x,F)

     else if (fctindx.eq.3) then

        call threebarf(3,dat,x,F)

     else if (fctindx.eq.4) then

        call threebarf(4,dat,x,F)

     else if (fctindx.eq.5) then

        call threebarf(5,dat,x,F)

     else if (fctindx.eq.6) then

        call threebarf(6,dat,x,F)

     else if (fctindx.eq.7) then

        call threebarf(7,dat,x,F)

     else if (fctindx.eq.8) then

        call threebarf(8,dat,x,F)

     else

        stop'wrong fctindx'

     end if

  else
     stop'Wrong Fn Number'

  end if

end subroutine calcf

  
subroutine calcdf(x,DIM,fct,df)
  use dimKrig, only: fctindx,DAT,mainprog
  implicit none
  integer :: DIM,fct,k
  real*8 :: x(DIM),df(DIM),fac,A,omeg

  real*8 :: rho, L, sigmay, pi, Fs, p, E,R, T,sigma_allow !Truss design parameters
  real*8:: tau_allow,M,V,B,D
  real*8::gamma
  real*8::tensile_sigma1_max,tensile_sigma2_max,tensile_sigma3_max
  real*8::comp_sigma1_max,comp_sigma2_max,comp_sigma3_max
  real*8::max_u_disp,max_v_disp,theta,pu,pv,u,sigma(dim)



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
!     if (dim.ne.3) stop'Wrong dimension'

     pi=4.0*atan(1.0)

!!$     ! Storing variables into DAT array to passthrough
!!$     ! Use this in your calling subprogam
!!$     
!!$     dat(1)=rho
!!$     dat(2)=L
!!$     dat(3)=E
!!$     dat(4)=P
!!$     dat(5)=sigmay
!!$     dat(6)=Fs

     if (mainprog) then
        rho=0.2836
        sigmay=36260.0
        p=25000.0
        L=5.0
        E=30e6
        Fs=1.0
     else

        !Read variables from DAT array to use in expressions

        rho=dat(1)
        L=dat(2)
        E=dat(3)
        P=dat(4)
        sigmay=dat(5)
        Fs=dat(6)

     end if

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
     !Tubular column
     !Thanks: Arora Section 3.7

!     if (dim.ne.3) stop'Wrong dimension'

     !Thanks: Arora Section 3.7

     pi=4.0*atan(1.0)

     if (mainprog) then

        rho=7.833e-6           !kg/m3
        p=10.0e6               !10 MN
        E=207000.0               !N/mm2
        sigma_allow=248.0        !N/mm2
        Fs=1.0

     else

        !Read variables from DAT array to use in expressions

        rho=dat(1)
        E=dat(2)
        P=dat(3)
        sigma_allow=dat(4)
        Fs=dat(5)

     end if

!!$
!!$     ! Store variables into DAT array to passthrough
!!$     ! Use this in your calling program
!!$     
!!$       dat(1)=rho
!!$       dat(2)=E
!!$       dat(3)=P
!!$       dat(4)=sigma_allow
!!$       dat(5)=Fs

     !define R,T

     R=x(1)
     T=x(2)
     L=x(3)

     if (fctindx.eq.0) then

        !---- OBJECTIVE FUNCTION
        df(1)=2.0*pi*rho*T*L
        df(2)=2.0*pi*rho*R*L
        df(3)=2.0*pi*rho*R*T

     else if (fctindx.eq.1) then

        !---- INEQUALITY CONSTRAINT 1

        df(1)= -P*FS / (2.0*pi*sigma_allow*(R**2)*T)
        df(2)= -P*FS / (2.0*pi*sigma_allow*R*T**2)
        df(3)= 0.0

     else if (fctindx.eq.2) then

        !---- INEQUALITY CONSTRAINT 2

        df(1)=-12.0*fs*P*L**2 / (E*T*(pi**3)*(R**4))
        df(2)=-4.0*P*FS*L**2 / (E*(T**2)*(pi**3)*(R**3))
        df(3)= 8.0*P*L*FS / (E*T*(pi**3)*(R**3))


     else
        print*, 'Wrong function index for this test case',fctindx
        stop

     end if


  else if (fct.eq.10) then ! Cantilever beam 
!     if (dim.ne.2) stop'Wrong dimension'

     ! Thanks: Section 3.8 Arora 
     if (mainprog) then

        M=40.0d6 !Nm
        V=150000.0 !N
        sigma_allow= 10.0d6 !N/m2
        tau_allow= 2.0d6 !N/m2     
        Fs=1.0

     else

        M=dat(1)
        V=dat(2)
        sigma_allow=dat(3)
        tau_allow=dat(4)
        Fs=dat(5)

     end if

!!$       !Store variables into DAT array to passthrough
!!$       ! Use this in your calling program
!!$       dat(1)=M
!!$       dat(2)=V
!!$       dat(3)=sigma_allow
!!$       dat(4)=tau_allow
!!$       dat(5)=Fs

     !define R,T

     B=x(1)
     D=x(2)

     if (fctindx.eq.0) then

        !---- OBJECTIVE FUNCTION
        df(1)=d
        df(2)=b

     else if (fctindx.eq.1) then

        !---- INEQUALITY CONSTRAINT 1

        df(1)= -(6.0*fs*M)/((b*d)**2*sigma_allow)
        df(2)= -(12.0*M*fs)/(b*(d**3)*sigma_allow)


     else if (fctindx.eq.2) then

        !---- INEQUALITY CONSTRAINT 2

        df(1)=-(3.0*V*fs)/(2.0*(b**2)*d*tau_allow)
        df(2)=-(3.0*V*FS)/(2.0*b*(d**2)*tau_allow)

     else if (fctindx.eq.3) then

        df(1) = -1.0*d*FS/(2.0*b**2)
        df(2) = 1.0*FS/(2.0*b) 
     else

        print*, 'Wrong function index for this test case',fctindx
        stop

     end if


  else if (fct.eq.11) then

!     if (dim.ne.3) stop'wrong dim for this problem'

     pi=4.0*atan(1.0)

     if (mainprog) then

        !Problem data and other constants
        L=10.0 !height ref
        E=1.0e7 !E
        gamma=0.1 !gamma
        theta=135.0*pi/180.0
        P=20000.0

        ! Max constraint values

        !Tensile
        tensile_sigma1_max=5000.0    ! psi tensile_sigma1_max=dat(6)      
        tensile_sigma2_max=20000.0    ! psi tensile_sigma2_max=dat(7)
        tensile_sigma3_max=5000.0    ! psi tensile_sigma3_max=dat(8)
        !Compressive
        comp_sigma1_max=5000.0    ! psi comp_sigma1_max=dat(9)
        comp_sigma2_max=20000.0   ! psi comp_sigma2_max=dat(10)
        comp_sigma3_max=5000.0   ! psi comp_sigma3_max=dat(11)
        !Displacement
        max_u_disp=0.005    ! in  max_u_disp=dat(12)
        max_v_disp=0.005    ! in  max_v_disp=dat(12)

        Fs=1.0      ! Factor of safety

     else

        L=dat(1)
        E=dat(2)
        gamma=dat(3)
        theta=dat(4)
        P=dat(5)

        ! Max constraint values

        !Tensile
        tensile_sigma1_max=dat(6)
        tensile_sigma2_max= dat(7)
        tensile_sigma3_max= dat(8)

        !Compressive
        comp_sigma1_max= dat(9)
        comp_sigma2_max= dat(10)
        comp_sigma3_max=dat(11)

        !Displacement
        max_u_disp= dat(12)
        max_v_disp= dat(13)
        Fs= dat(14)
     end if

!!$     ! Use this in caling program and pass to pcestimate subtoutine
!!$     !Problem data and other constants
!!$     dat(1)=10.0 !height ref
!!$     dat(2)=1.0e7 !E
!!$     dat(3)=0.1 !gamma
!!$     dat(4)=45.0*pi/180.0
!!$     dat(5)=20000.0
!!$
!!$     ! Max constraint values
!!$
!!$     !Tensile
!!$     dat(6)=5000.0    ! psi tensile_sigma1_max=dat(6)      
!!$     dat(7)=20000.0    ! psi tensile_sigma2_max=dat(7)
!!$     dat(8)=5000.0    ! psi tensile_sigma3_max=dat(8)
!!$     !Compressive
!!$     dat(9)=5000.0    ! psi comp_sigma1_max=dat(9)
!!$     dat(10)=20000.0   ! psi comp_sigma2_max=dat(10)
!!$     dat(11)=5000.0   ! psi comp_sigma3_max=dat(11)
!!$     !Displacement
!!$     dat(12)=0.005    ! in  max_u_disp=dat(12)
!!$     dat(13)=0.005    ! in  max_v_disp=dat(12)
!!$     dat(14)=1.0      ! Factor of safety
!!$

     pu=P*cos(theta)
     pv=P*sin(theta)

     if (fctindx.eq.0) then

        df(1)= L*sqrt(2.0)*gamma
        df(2) = L*gamma
        df(3) = L*sqrt(2.0)*gamma

     else if (fctindx.eq.1) then


        df(1)=((x(2) + 2.0**(1.0/2.0)*x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(2)=((x(1) + x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*pu)/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

        df(3)=((x(2) + 2.0**(1.0/2.0)*x(1))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu + pv)/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))    


     else if (fctindx.eq.2) then
        df(1)=((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(3)))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu - pv)/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

        df(2)=((x(1) + x(3))*(pu*(x(1) - x(3)) - pv*(x(1) + x(3))))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(3)=(pu + pv)/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) + ((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(1)))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)


     else if (fctindx.eq.3) then

        df(1)=((x(2) + 2.0**(1.0/2.0)*x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu - pv)/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

        df(2)=((x(1) + x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*pu)/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

        df(3)=((x(2) + 2.0**(1.0/2.0)*x(1))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

     else if (fctindx.eq.4) then

        df(1)=-((x(2) + 2.0**(1.0/2.0)*x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(2)=(2.0**(1.0/2.0)*pu)/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(1) + x(3))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(3)=(pu + pv)/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(2) + 2.0**(1.0/2.0)*x(1))*(x(3)*pu + x(3)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

     else if (fctindx.eq.5) then


        df(1)=(pu - pv)/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(3)))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(2)=-((x(1) + x(3))*(pu*(x(1) - x(3)) - pv*(x(1) + x(3))))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(3)=- (pu + pv)/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((pu*(x(1) - x(3)) - pv*(x(1) + x(3)))*(x(2) + 2.0**(1.0/2.0)*x(1)))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)


     else if (fctindx.eq.6) then

        df(1)=(pu - pv)/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(2) + 2.0**(1.0/2.0)*x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(2)=(2.0**(1.0/2.0)*pu)/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((x(1) + x(3))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(3)=-((x(2) + 2.0**(1.0/2.0)*x(1))*(x(1)*pu - x(1)*pv + 2.0**(1.0/2.0)*x(2)*pu))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

     else if (fctindx.eq.7) then
        df(1) = ((1.0/max_u_disp)*L*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1)*pu + x(3)*pu - x(1)*pv + x(3)*pv + 2.0*2.0**(1.0/2.0)*x(2)*pu))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - ((1.0/max_u_disp)*L*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

        df(2) = ((1.0/max_u_disp)*L*(x(1) + x(3))*(x(1)*pu + x(3)*pu - x(1)*pv + x(3)*pv + 2.0*2.0**(1.0/2.0)*x(2)*pu))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (400.0*2.0**(1.0/2.0)*L*pu)/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

        df(3) = ((1.0/max_u_disp)*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1)*pu + x(3)*pu - x(1)*pv + x(3)*pv + 2.0*2.0**(1.0/2.0)*x(2)*pu))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - ((1.0/max_u_disp)*L*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3)))

     else if (fctindx.eq.8) then

        df(1)= -((1.0/max_v_disp)*L*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((1.0/max_v_disp)*L*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(3)*pu - x(1)*pu + x(1)*pv + x(3)*pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(2)= -((1.0/max_v_disp)*L*(x(1) + x(3))*(x(3)*pu - x(1)*pu + x(1)*pv + x(3)*pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        df(3)= ((1.0/max_v_disp)*L*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))) - ((1.0/max_v_disp)*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(3)*pu - x(1)*pu + x(1)*pv + x(3)*pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)


     else

        stop'wrong function index for this case'


     endif


  else if (fct.eq.12) then

!     if (dim.ne.6) stop'wrong dim for this problem'

     pi=4.0*atan(1.0)
     if (mainprog) then 
        ! Use these settings if the program is called from main.f90.
        ! If PC is used as library DAT is passed as an input vector from calling program such as IPOPT

        !Problem data and other constants
        dat(1)=10.0 !height ref
        dat(2)=1.0e7 !E
        dat(3)=0.1 !gamma
        dat(4)=45.0*pi/180.0
        dat(5)=20000.0

        ! Max constraint values

        !Tensile
        dat(6)=5000.0    ! psi tensile_sigma1_max=dat(6)      
        dat(7)=20000.0    ! psi tensile_sigma2_max=dat(7)
        dat(8)=5000.0    ! psi tensile_sigma3_max=dat(8)
        !Compressive
        dat(9)=5000.0    ! psi comp_sigma1_max=dat(9)
        dat(10)=20000.0   ! psi comp_sigma2_max=dat(10)
        dat(11)=5000.0   ! psi comp_sigma3_max=dat(11)
        !Displacement
        dat(12)=0.005    ! in  max_u_disp=dat(12)
        dat(13)=0.005    ! in  max_v_disp=dat(12)
        dat(14)=1.0      ! Factor of safety
     end if

     if (fctindx.eq.0) then
        call threebardf(0,dat,x,df)

     else if (fctindx.eq.1) then
        call threebardf(1,dat,x,df)

     else if (fctindx.eq.2) then
        call threebardf(2,dat,x,df)

     else if (fctindx.eq.3) then
        call threebardf(3,dat,x,df)

     else if (fctindx.eq.4) then
        call threebardf(4,dat,x,df)

     else if (fctindx.eq.5) then
        call threebardf(5,dat,x,df)

     else if (fctindx.eq.6) then
        call threebardf(6,dat,x,df)

     else if (fctindx.eq.7) then
        call threebardf(7,dat,x,df)

     else if (fctindx.eq.8) then
        call threebardf(8,dat,x,df)
     else

        stop'wrong fctindx'

     end if


  else

     stop'Unsupported function number'

  end if

end subroutine calcdf

subroutine calcd2f(x,DIM,fct,d2f)
  use dimKrig, only: fctindx,DAT,mainprog
  implicit none
  integer :: DIM,fct,j,k
  real*8 :: x(DIM),d2f(DIM,DIM),fac,A,omeg,P,rho,Fs,sigmay
  real*8::gamma,L,E,pi,R, T,sigma_allow !Truss design parameters
  real*8::tensile_sigma1_max,tensile_sigma2_max,tensile_sigma3_max
  real*8::comp_sigma1_max,comp_sigma2_max,comp_sigma3_max
  real*8::max_u_disp,max_v_disp,theta,pu,pv,u,sigma(dim)
  real*8:: tau_allow,M,V,B,D

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
!     if (dim.ne.3) stop'Wrong dimension'

     pi=4.0*atan(1.0)

!!$     ! Storing variables into DAT array to passthrough
!!$     ! Use this in your calling subprogam
!!$     
!!$     dat(1)=rho
!!$     dat(2)=L
!!$     dat(3)=E
!!$     dat(4)=P
!!$     dat(5)=sigmay
!!$     dat(6)=Fs

     if (mainprog) then
        rho=0.2836
        sigmay=36260.0
        p=25000.0
        L=5.0
        E=30e6
        Fs=1.0
     else

        !Read variables from DAT array to use in expressions

        rho=dat(1)
        L=dat(2)
        E=dat(3)
        P=dat(4)
        sigmay=dat(5)
        Fs=dat(6)

     end if
     d2f(:,:)=0.0d0
     if (fctindx.eq.0) then
        !---- OBJECTIVE FUNCTION
        d2f(1,1) = 0.0d0
        d2f(2,1) = 0.0d0
        d2f(3,1) = rho*x(2)*x(3) / sqrt(L**2+x(3)**2)

     else if (fctindx.eq.1) then
        !---- INEQUALITY CONSTRAINTS

        d2f(2,3)= p*Fs*L**2/ sqrt(L**2/x(3)**2+1.0) / (x(2)**2*x(3)**3*sigmay)
        d2f(3,3)= -p*Fs*L**2/(x(2)*sigmay)*(L**2/sqrt(L**2/x(3)**2+1.0)**3/x(3)**6 - 3.0/x(3)**4/sqrt(L**2/x(3)**2+1.0))

        !copy to the upper triangle
        d2f(3,2)=d2f(2,3)

     else if (fctindx.eq.2) then  

        d2f(1,3)= p*Fs*L / (x(1)**2*x(3)**2*sigmay)
        d2f(3,3)= 2.0*p*Fs*L / (x(1)*x(3)**3*sigmay)

        !copy to upper

        d2f(3,1)=d2f(1,3)

     else if (fctindx.eq.3) then

        d2f(1,3)= 8.0*p*Fs*L**3 / (pi*x(3)**2*E*x(1)**3)
        d2f(3,3)= 8.0*p*Fs*L**3 / (pi*x(3)**3*E*x(1)**2)

        !copy to upper Triangle

        d2f(3,1)=d2f(1,3)

     end if

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

!     if (dim.ne.3) stop'Wrong dimension'

     !Thanks: Arora Section 3.7

     pi=4.0*atan(1.0)

     if (mainprog) then

        rho=7.833e-6           !kg/m3
        p=10.0e6               !10 MN
        E=207000.0               !N/mm2
        sigma_allow=248.0        !N/mm2
        Fs=1.0

     else

        !Read variables from DAT array to use in expressions

        rho=dat(1)
        E=dat(2)
        P=dat(3)
        sigma_allow=dat(4)
        Fs=dat(5)

     end if

!!$
!!$     ! Store variables into DAT array to passthrough
!!$     ! Use this in your calling program
!!$     
!!$       dat(1)=rho
!!$       dat(2)=E
!!$       dat(3)=P
!!$       dat(4)=sigma_allow
!!$       dat(5)=Fs

     !define R,T

     R=x(1)
     T=x(2)
     L=x(3)

     if (fctindx.eq.0) then

        !---- OBJECTIVE FUNCTION
        !lower triangle
        d2f(1,1)=0.0
        d2f(2,2)=0.0
        d2f(3,3)=0.0

        d2f(2,1)=2.0*pi*rho*L
        d2f(3,2)=2.0*pi*rho*R

        d2f(3,1)=2.0*pi*rho*T

        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)


     else if (fctindx.eq.1) then

        !---- INEQUALITY CONSTRAINT 1

        d2f(1,1)= P*FS/(pi*(R**3)*T*sigma_allow)
        d2f(2,2)= P*FS/(pi*R*(T**3)*sigma_allow)
        d2f(3,3)= 0.0

        d2f(2,1)= P*FS/(2.0*pi*(R**2)*(T**2)*sigma_allow)
        d2f(3,2)= 0.0

        d2f(3,1)= 0.0

        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)

     else if (fctindx.eq.2) then

        !---- INEQUALITY CONSTRAINT 2

        !lower triangle

        d2f(1,1)= 48.0*FS*(L**2)*P/(E*T*(pi**3)*(R**5))
        d2f(2,2)= 8.0*FS*(L**2)*P/(E*(pi*R*T)**3)
        d2f(3,3)= 8.0*P*FS/(E*T*(pi*R)**3)

        d2f(2,1)= (12.0*L**2*P*FS)/(E*pi**3*R**4*T**2)
        d2f(3,2)= -8.0*L*P*FS/((pi**3)*E*(R**3)*(T**2))

        d2f(3,1)= -(24.0*L*P*FS)/(E*pi**3*R**4*T)

        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)

     else

        print*, 'Wrong function index for this test case',fctindx
        stop

     end if

  else if (fct.eq.10) then !cantilever beam
!     if (dim.ne.2) stop'Wrong dimension'

     ! Thanks: Section 3.8 Arora 
     if (mainprog) then

        M=40.0d6 !Nm
        V=150000.0 !N
        sigma_allow= 10.0d6 !N/m2
        tau_allow= 2.0d6 !N/m2     
        Fs=1.0

     else

        M=dat(1)
        V=dat(2)
        sigma_allow=dat(3)
        tau_allow=dat(4)
        Fs=dat(5)

     end if

!!$       !Store variables into DAT array to passthrough
!!$       ! Use this in your calling program
!!$       dat(1)=M
!!$       dat(2)=V
!!$       dat(3)=sigma_allow
!!$       dat(4)=tau_allow
!!$       dat(5)=Fs

     !define R,T

     B=x(1)
     D=x(2)


     if (fctindx.eq.0) then

        !---- OBJECTIVE FUNCTION

        !lower triangle
        d2f(1,1)=0.0
        d2f(2,2)=0.0

        d2f(2,1)=1.0

        !copy to upper triangle

        d2f(1,2)=d2f(2,1)

     else if (fctindx.eq.1) then

        !---- INEQUALITY CONSTRAINT 1

        d2f(1,1)= (12.0*M*fs)/((b**3)*(d**2)*sigma_allow)
        d2f(2,2)= (36.0*M*fs)/(b*(d**4)*sigma_allow)

        d2f(2,1)= (12.0*M*Fs)/((b**2)*(d**3)*sigma_allow)

        !copy to upper triangle

        d2f(1,2)=d2f(2,1)

     else if (fctindx.eq.2) then

        !---- INEQUALITY CONSTRAINT 2

        !lower triangle

        d2f(1,1)= (3.0*V*fs)/((b**3)*d*tau_allow)
        d2f(2,2)= (3.0*V*fs)/(b*(d**3)*tau_allow)

        d2f(2,1)=(3.0*V*fs)/(2.0*(b**2)*(d**2)*tau_allow)

        !copy to upper triangle

        d2f(1,2)=d2f(2,1)

     else if (fctindx.eq.3) then

        !lower triangle

        d2f(1,1)= d*FS/(B**3)
        d2f(2,2)= 0.0             
        d2f(2,1)=-1.0*FS/(2.0*b**2)

        !copy 

        d2f(1,2)=d2f(2,1)

     else

        print*, 'Wrong function index for this test case',fctindx
        stop

     end if !fctindx

  else if (fct.eq.11) then

!     if (dim.ne.3) stop'wrong dim for this problem'

     pi=4.0*atan(1.0)

     if (mainprog) then

        !Problem data and other constants
        L=10.0 !height ref
        E=1.0e7 !E
        gamma=0.1 !gamma
        theta=135.0*pi/180.0
        P=20000.0

        ! Max constraint values

        !Tensile
        tensile_sigma1_max=5000.0    ! psi tensile_sigma1_max=dat(6)      
        tensile_sigma2_max=20000.0    ! psi tensile_sigma2_max=dat(7)
        tensile_sigma3_max=5000.0    ! psi tensile_sigma3_max=dat(8)
        !Compressive
        comp_sigma1_max=5000.0    ! psi comp_sigma1_max=dat(9)
        comp_sigma2_max=20000.0   ! psi comp_sigma2_max=dat(10)
        comp_sigma3_max=5000.0   ! psi comp_sigma3_max=dat(11)
        !Displacement
        max_u_disp=0.005    ! in  max_u_disp=dat(12)
        max_v_disp=0.005    ! in  max_v_disp=dat(12)

        Fs=1.0      ! Factor of safety

     else

        L=dat(1)
        E=dat(2)
        gamma=dat(3)
        theta=dat(4)
        P=dat(5)

        ! Max constraint values

        !Tensile
        tensile_sigma1_max=dat(6)
        tensile_sigma2_max= dat(7)
        tensile_sigma3_max= dat(8)

        !Compressive
        comp_sigma1_max= dat(9)
        comp_sigma2_max= dat(10)
        comp_sigma3_max=dat(11)

        !Displacement
        max_u_disp= dat(12)
        max_v_disp= dat(13)
        Fs= dat(14)
     end if

!!$     ! Use this in caling program and pass to pcestimate subtoutine
!!$     !Problem data and other constants
!!$     dat(1)=10.0 !height ref
!!$     dat(2)=1.0e7 !E
!!$     dat(3)=0.1 !gamma
!!$     dat(4)=45.0*pi/180.0
!!$     dat(5)=20000.0
!!$
!!$     ! Max constraint values
!!$
!!$     !Tensile
!!$     dat(6)=5000.0    ! psi tensile_sigma1_max=dat(6)      
!!$     dat(7)=20000.0    ! psi tensile_sigma2_max=dat(7)
!!$     dat(8)=5000.0    ! psi tensile_sigma3_max=dat(8)
!!$     !Compressive
!!$     dat(9)=5000.0    ! psi comp_sigma1_max=dat(9)
!!$     dat(10)=20000.0   ! psi comp_sigma2_max=dat(10)
!!$     dat(11)=5000.0   ! psi comp_sigma3_max=dat(11)
!!$     !Displacement
!!$     dat(12)=0.005    ! in  max_u_disp=dat(12)
!!$     dat(13)=0.005    ! in  max_v_disp=dat(12)
!!$     dat(14)=1.0      ! Factor of safety
!!$

     pu=P*cos(theta)
     pv=P*sin(theta)

     d2f(:,:)=0.0

     if (fctindx.eq.0) then

        d2f(:,:)=0.0

     else  if (fctindx.eq.1) then

        !lower triangle

        d2f(1,1)= -(2.0*(x(2) + 2.0**(1.0/2.0)*x(3))**2.0*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(tensile_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)

        d2f(2,2)=-((2.0*(x(1) + x(3))**2.0*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3 - (2.0*2.0**(1.0/2.0)*pu*(x(1) + x(3)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)/tensile_sigma1_max

        d2f(3,3)=-((2.0*(x(2) + 2.0**(1.0/2.0)*x(1))**2.0*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3 - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(pu + pv))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)/tensile_sigma1_max


        d2f(2,1)=((pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 + (2.0**(1.0/2.0)*pu*(x(2) + 2.0**(1.0/2.0)*x(3)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 - (2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1) + x(3))*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)/tensile_sigma1_max

        d2f(3,2)=((pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 + ((pu + pv)*(x(1) + x(3)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 + (2.0**(1.0/2.0)*pu*(x(2) + 2.0**(1.0/2.0)*x(1)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1) + x(3))*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)/tensile_sigma1_max


        d2f(3,1)=((2.0**(1.0/2.0)*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 + ((x(2) + 2.0**(1.0/2.0)*x(3))*(pu + pv))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)/tensile_sigma1_max


        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)

     else  if (fctindx.eq.2) then

        !lower triangle
        d2f(1,1)=((2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu - pv))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 + (2.0*(x(2) + 2.0**(1.0/2.0)*x(3))**2.0*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)/tensile_sigma2_max

        d2f(2,2)=(2.0*(x(1) + x(3))**2.0*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)

        d2f(3,3)=-((2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(pu + pv))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))**2.0*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)/tensile_sigma2_max


        d2f(2,1)=((x(1) + x(3))*(pu - pv))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pv*(x(1) + x(3)) - pu*(x(1) - x(3)))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1) + x(3))*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(tensile_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)

        d2f(3,2)=-((pv*(x(1) + x(3)) - pu*(x(1) - x(3)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 + ((pu + pv)*(x(1) + x(3)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1) + x(3))*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)/tensile_sigma2_max


        d2f(3,1)=(((x(2) + 2.0**(1.0/2.0)*x(1))*(pu - pv))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 - ((x(2) + 2.0**(1.0/2.0)*x(3))*(pu + pv))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 - (2.0**(1.0/2.0)*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 + (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(2) + 2.0**(1.0/2.0)*x(3))*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)/tensile_sigma2_max


        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)
     else  if (fctindx.eq.3) then

        !lower triangle
        d2f(1,1)=((2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu - pv))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 - (2.0*(x(2) + 2.0**(1.0/2.0)*x(3))**2.0*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)/tensile_sigma3_max

        d2f(2,2)=-((2.0*(x(1) + x(3))**2.0*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3 - (2.0*2.0**(1.0/2.0)*pu*(x(1) + x(3)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)/tensile_sigma3_max

        d2f(3,3)=-(2.0*(x(2) + 2.0**(1.0/2.0)*x(1))**2.0*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)


        d2f(2,1)=((pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 + ((x(1) + x(3))*(pu - pv))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 + (2.0**(1.0/2.0)*pu*(x(2) + 2.0**(1.0/2.0)*x(3)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2 - (2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1) + x(3))*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)/tensile_sigma3_max

        d2f(3,2)=(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (2.0**(1.0/2.0)*pu*(x(2) + 2.0**(1.0/2.0)*x(1)))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1) + x(3))*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)


        d2f(3,1)=(2.0**(1.0/2.0)*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + ((x(2) + 2.0**(1.0/2.0)*x(1))*(pu - pv))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(tensile_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)


        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)
     else  if (fctindx.eq.4) then

        !lower triangle
        d2f(1,1)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(3))**2.0*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)

        d2f(2,2)=(2.0*(x(1) + x(3))**2.0*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - (2.0*2.0**(1.0/2.0)*pu*(x(1) + x(3)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        d2f(3,3)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(1))**2.0*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(pu + pv))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)


        d2f(2,1)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1) + x(3))*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - (2.0**(1.0/2.0)*pu*(x(2) + 2.0**(1.0/2.0)*x(3)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        d2f(3,2)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1) + x(3))*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - ((pu + pv)*(x(1) + x(3)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*pu*(x(2) + 2.0**(1.0/2.0)*x(1)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        d2f(3,1)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - ((x(2) + 2.0**(1.0/2.0)*x(3))*(pu + pv))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*(pu*x(3) + pv*x(3) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma1_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)


        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)
     else  if (fctindx.eq.5) then

        !lower triangle
        d2f(1,1)=- (2.0*(x(2) + 2.0**(1.0/2.0)*x(3))**2.0*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - (2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu - pv))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        d2f(2,2)=-(2.0*(x(1) + x(3))**2.0*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)

        d2f(3,3)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(pu + pv))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))**2.0*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)


        d2f(2,1)=(pv*(x(1) + x(3)) - pu*(x(1) - x(3)))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - ((x(1) + x(3))*(pu - pv))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1) + x(3))*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)

        d2f(3,2)=(pv*(x(1) + x(3)) - pu*(x(1) - x(3)))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + ((pu + pv)*(x(1) + x(3)))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1) + x(3))*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)


        d2f(3,1)=((x(2) + 2.0**(1.0/2.0)*x(3))*(pu + pv))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (2.0**(1.0/2.0)*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - ((x(2) + 2.0**(1.0/2.0)*x(1))*(pu - pv))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(2) + 2.0**(1.0/2.0)*x(3))*(pv*(x(1) + x(3)) - pu*(x(1) - x(3))))/(comp_sigma2_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)


        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)
     else  if (fctindx.eq.6) then

        !lower triangle
        d2f(1,1)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(3))**2.0*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - (2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu - pv))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        d2f(2,2)=(2.0*(x(1) + x(3))**2.0*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - (2.0*2.0**(1.0/2.0)*pu*(x(1) + x(3)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        d2f(3,3)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(1))**2.0*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)


        d2f(2,1)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1) + x(3))*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - ((x(1) + x(3))*(pu - pv))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*pu*(x(2) + 2.0**(1.0/2.0)*x(3)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)

        d2f(3,2)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1) + x(3))*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - (2.0**(1.0/2.0)*pu*(x(2) + 2.0**(1.0/2.0)*x(1)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)


        d2f(3,1)=(2.0*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - ((x(2) + 2.0**(1.0/2.0)*x(1))*(pu - pv))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*(pu*x(1) - pv*x(1) + 2.0**(1.0/2.0)*pu*x(2)))/(comp_sigma3_max*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2)
        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)
     else  if (fctindx.eq.7) then

        !lower triangle
        d2f(1,1)=((2.0*L*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(3))**2.0*(pu*x(1) + pu*x(3) - pv*x(1) + pv*x(3) + 2.0*2.0**(1.0/2.0)*pu*x(2)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3))/max_u_disp

        d2f(2,2)=-((2.0*L*(x(1) + x(3))**2.0*(pu*x(1) + pu*x(3) - pv*x(1) + pv*x(3) + 2.0*2.0**(1.0/2.0)*pu*x(2)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - (4*2.0**(1.0/2.0)*L*pu*(x(1) + x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2))/max_u_disp

        d2f(3,3)=((2.0*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(1))**2.0*(pu*x(1) + pu*x(3) - pv*x(1) + pv*x(3) + 2.0*2.0**(1.0/2.0)*pu*x(2)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3))/max_u_disp


        d2f(2,1)=((L*(pu*x(1) + pu*x(3) - pv*x(1) + pv*x(3) + 2.0*2.0**(1.0/2.0)*pu*x(2)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (L*(x(1) + x(3))*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (2.0*2.0**(1.0/2.0)*L*pu*(x(2) + 2.0**(1.0/2.0)*x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1) + x(3))*(pu*x(1) + pu*x(3) - pv*x(1) + pv*x(3) + 2.0*2.0**(1.0/2.0)*pu*x(2)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3))/max_u_disp

        d2f(3,2)=((L*(pu*x(1) + pu*x(3) - pv*x(1) + pv*x(3) + 2.0*2.0**(1.0/2.0)*pu*x(2)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (L*(pu + pv)*(x(1) + x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (2.0*2.0**(1.0/2.0)*L*pu*(x(2) + 2.0**(1.0/2.0)*x(1)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1) + x(3))*(pu*x(1) + pu*x(3) - pv*x(1) + pv*x(3) + 2.0*2.0**(1.0/2.0)*pu*x(2)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3))/max_u_disp


        d2f(3,1)=((L*(x(2) + 2.0**(1.0/2.0)*x(1))*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (L*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (2.0**(1.0/2.0)*L*(pu*x(1) + pu*x(3) - pv*x(1) + pv*x(3) + 2.0*2.0**(1.0/2.0)*pu*x(2)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu*x(1) + pu*x(3) - pv*x(1) + pv*x(3) + 2.0*2.0**(1.0/2.0)*pu*x(2)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3))/max_u_disp


        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)
     else  if (fctindx.eq.8) then

        !lower triangle
        d2f(1,1)=((2.0*L*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(3))**2.0*(pu*x(3) - pu*x(1) + pv*x(1) + pv*x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3))/max_v_disp

        d2f(2,2)=(2.0*L*(x(1) + x(3))**2.0*(pu*x(3) - pu*x(1) + pv*x(1) + pv*x(3)))/(E*max_v_disp*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)

        d2f(3,3)=((2.0*L*(x(2) + 2.0**(1.0/2.0)*x(1))**2.0*(pu*x(3) - pu*x(1) + pv*x(1) + pv*x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3) - (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2))/max_v_disp


        d2f(2,1)=(L*(x(1) + x(3))*(pu - pv))/(E*max_v_disp*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (L*(pu*x(3) - pu*x(1) + pv*x(1) + pv*x(3)))/(E*max_v_disp*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(3))*(x(1) + x(3))*(pu*x(3) - pu*x(1) + pv*x(1) + pv*x(3)))/(E*max_v_disp*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3)

        d2f(3,2)=-((L*(pu*x(3) - pu*x(1) + pv*x(1) + pv*x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (L*(pu + pv)*(x(1) + x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(1) + x(3))*(pu*x(3) - pu*x(1) + pv*x(1) + pv*x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3))/max_v_disp


        d2f(3,1)=((L*(x(2) + 2.0**(1.0/2.0)*x(1))*(pu - pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (2.0**(1.0/2.0)*L*(pu*x(3) - pu*x(1) + pv*x(1) + pv*x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) - (L*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu + pv))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**2) + (2.0*L*(x(2) + 2.0**(1.0/2.0)*x(1))*(x(2) + 2.0**(1.0/2.0)*x(3))*(pu*x(3) - pu*x(1) + pv*x(1) + pv*x(3)))/(E*(x(1)*x(2) + x(2)*x(3) + 2.0**(1.0/2.0)*x(1)*x(3))**3))/max_v_disp

        !copy to upper triangle

        d2f(1,2)=d2f(2,1)
        d2f(2,3)=d2f(3,2)

        d2f(1,3)=d2f(3,1)

     else

        print *,'Warning',fct
        stop'Invalid fct for Hessian'

     end if

  else 

     stop'Wrong Fn number'

  end if

end subroutine calcd2f




