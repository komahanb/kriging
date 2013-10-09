        subroutine SCF(iscf,hstat,xf,xb,tt,pp,v)
        implicit none

        integer ,         intent(in)  :: iscf,hstat
        double precision, intent(in)  :: xf,xb,tt,pp
        double precision, intent(out) :: v
        double precision :: xx

        xx = dabs(xf-xb)

        if(iscf.eq.0)then      ! Gaussian

          if(pp.ne.2.d0)stop'pp!=2 in gauss_scf'
          v   = exp(-1.d0*tt*xx**2.d0)

        else if(iscf.eq.2)then ! Cubic Spline

          if(      xx.lt.(1.d0/(2.d0*tt))                       )then
             v   = 1.d0 - 6.d0*(xx*tt)**2 + 6.d0*(xx*tt)**3
          else if( xx.ge.(1.d0/(2.d0*tt)) .and. xx.lt.(1.d0/tt) )then
             v   = 2.d0*(1.d0-xx*tt)**3 
          else if( xx.ge.(1.d0/tt)                              )then
             v   = 0.d0
          end if

        else if(iscf.eq.3)then ! Wendland C2

          if(xx.le.1.d0/tt)then
             v = ((1.d0-tt*xx)**4)*(4.d0*tt*xx+1.d0)
          else
             v = 0.d0
          end if

        else if(iscf.eq.4)then ! Wendland C4

          if(xx.le.1.d0/tt)then
             v = ((1.d0-tt*xx)**6)*(35.d0*(tt*xx)**2 + 18.d0*tt*xx + 3.d0)/3.d0
          else
             v = 0.d0
          end if

        else if(iscf.eq.5)then ! Matern function

          if(hstat.eq.0)then !nu=1/2
             v = exp(-tt*xx) 
          else if(hstat.eq.1)then !nu=3/2
             v = (1.0+sqrt(3.0)*tt*xx)*exp(-sqrt(3.0)*tt*xx)
          else !nu=5/2
             v = (1.0+sqrt(5.0)*tt*xx+5.0/3.0*(tt*xx)**2)*exp(-sqrt(5.0)*tt*xx)
          end if

        else
          stop'unknown iSCF'
        end if

        end subroutine SCF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine SCF_OUT(iscf,xf,xb,tt,pp,v)
        use dimKrig, only:hstat
        implicit none
       
        integer ,         intent(in)  :: iscf
        double precision, intent(in)  :: xf,xb,tt,pp
        double precision, intent(out) :: v
        double precision :: xx
        ! Gaussian is different from above, power of that
        ! The above is only for Automatic Differentiation

        xx = dabs(xf-xb)
        if(iscf.eq.0.or.iscf.eq.1)then      ! Gaussian
          v   = exp(-1.d0*tt*((xf-xb)**pp))

        else if(iscf.eq.2)then ! Cubic Spline
          if(      xx.lt.(1.d0/(2.d0*tt))                       )then
             v   = 1.d0 - 6.d0*(xx*tt)**2 + 6.d0*(xx*tt)**3
          else if( xx.ge.(1.d0/(2.d0*tt)) .and. xx.lt.(1.d0/tt) )then
             v   = 2.d0*(1.d0-xx*tt)**3 
          else if( xx.ge.(1.d0/tt)                              )then
             v   = 0.d0
          end if

        else if(iscf.eq.3)then ! Wendland C2
          if(xx.le.1.d0/tt)then
             v = ((1.d0-tt*xx)**4)*(4.d0*tt*xx+1.d0)
          else
             v = 0.d0
          end if

        else if(iscf.eq.4)then ! Wendland C4
          if(xx.le.1.d0/tt)then
             v = ((1.d0-tt*xx)**6)*(35.d0*(tt*xx)**2 + 18.d0*tt*xx + 3.d0)/3.d0
          else
             v = 0.d0
          end if

        else if(iscf.eq.5)then ! Matern function

          if(hstat.eq.0)then !nu=1/2
             !v = exp(-tt*xx)
             if(xx.le.1.d0/tt)then
                v = ((1.d0-tt*xx)**6)*(35.d0*(tt*xx)**2 + 18.d0*tt*xx + 3.d0)/3.d0
             else
                v = 0.d0
             end if
          else if(hstat.eq.1)then !nu=3/2
             v = (1.0+sqrt(3.0)*tt*xx)*exp(-sqrt(3.0)*tt*xx)
          else !nu=5/2
             v = (1.0+sqrt(5.0)*tt*xx+5.0/3.0*(tt*xx)**2)*exp(-sqrt(5.0)*tt*xx)
          end if

        else
          stop'unknown iSCF'
        end if

        end subroutine SCF_OUT
