        subroutine makeMatrixPetitR(mode,idif,jdif, &
                                    ndim,nsize,rsample,iscf,mode_dck, &
                                    lmax,tdim,ict_sample,             &
                                    train,xin,nparent,dxadd,inf,hvec, &
                                    theta,power,ccrf,r)
          use dimKrig, only:hstat
! not use dimKrig only here for meta
! mode = 0 : r
!        1 : dr/dX_idif
!        2 : d2r/dX_idif/dX_jdif
        implicit none
        integer, intent(in) :: mode,idif,jdif
        integer, intent(in) :: ndim,nsize,rsample,iscf,mode_dck
        integer, intent(in) :: tdim,lmax
        double precision, intent(in) :: ccrf
        double precision, dimension(rsample,ndim),   intent(in) :: train
        double precision, dimension(ndim),           intent(in) :: xin
        double precision, dimension(ndim,tdim),      intent(in) :: theta,power
        integer,          dimension(rsample),        intent(in) :: nparent
        double precision, dimension(rsample),        intent(in) :: dxadd
        character(len=6), dimension(rsample),        intent(in) :: inf
        double precision, dimension(rsample,ndim,2), intent(in) :: hvec
        integer,          dimension(0:10,0:10),      intent(in) :: ict_sample
        double precision, dimension(nsize,1), intent(out) :: r

        integer :: i,k,l,m,t
        integer :: l1,idf,idg,idh,idv,nhes
        integer :: ict
        double precision :: xf,xb,tt,pp,value
        double precision :: scf,dscf_f1,dscf_b1
        double precision :: d2scf_f2,d2scf_f1b1,d2scf_b2
        double precision :: d3scf_f2b1,d3scf_f1b2,d4scf_f2b2
        character(len=2) :: cl
        double precision, dimension(ndim,ndim) :: H,Hfv
        double precision, dimension(ndim,1)    :: HV,V

        if(mode.ge.1)then
          if(idif.le.0.or.idif.gt.ndim)stop'idif in petitR'
        end if
        if(mode.eq.2)then
          if(jdif.le.0.or.jdif.gt.ndim)stop'jdif in petitR'
        end if

        nhes = ndim*(ndim+1)/2
        r = 1.d0

        do 100 i=1,rsample
          ! l1 is the fidelity level
          do l=1,lmax
           write(cl,105)l
105        format(i1,'_')
           if(cl.eq.inf(i)(1:2))then
             l1 = l
             go to 101
           end if
          end do
          stop'unknown level of i in PetitR'
101       continue
          ! which theta/power
          ict = 0
          do 110 l=1,1
          do 110 m=l,lmax
            ict = ict + 1
            if(l1.eq.m)then
              t = ict
              go to 111
            end if
110       continue
          stop'unknown tdim in PetitR'
111       continue
          ! which index on PetitR
          call find_index(i,l1,idf,idg,idh,idv, &
                          ndim,nhes,nsize,mode_dck,rsample, &
                          ict_sample,inf)
          Hfv = 1.d0
          do 120 k=1,ndim
            xf = train(i,k)
            xb = xin(k)
            tt = theta(k,t)
            pp = power(k,t)
            if(iscf.eq.0.or.iscf.eq.1)then
              if(pp.lt.1.d0)pp = 1.d0
              if(pp.gt.2.d0)pp = 2.d0
            end if
            ! AD
            if(mode.eq.0)then
              call SCF_DF_DF(      iscf,hstat,xf,1.d0,1.d0,xb,tt,pp, &
                                   scf,dscf_f1,d2scf_f2)
            else if(mode.eq.1)then
             if(     inf(i)(3:6).eq.'FH  '.or.inf(i)(3:6).eq.'FGH '.or. &
                     inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv')then
              call SCF_DF_DF_DB(   iscf,hstat,xf,1.d0,1.d0,xb,1.d0,tt,pp, &
                                   scf,dscf_f1,d2scf_f2,d3scf_f2b1)
              call SCF_DF_DB(      iscf,hstat,xf,1.d0,xb,1.d0,tt,pp, &
                                   scf,dscf_f1,d2scf_f1b1)
              call SCF_DB(         iscf,hstat,xf,xb,1.d0,tt,pp, &
                                   scf,dscf_b1)
             else if(inf(i)(3:6).eq.'FG  ')then
              call SCF_DF_DB(      iscf,hstat,xf,1.d0,xb,1.d0,tt,pp, &
                                   scf,dscf_f1,d2scf_f1b1)
              call SCF_DB(         iscf,hstat,xf,xb,1.d0,tt,pp, &
                                   scf,dscf_b1)
             else if(inf(i)(3:6).eq.'F   ')then
              call SCF_DB(         iscf,hstat,xf,xb,1.d0,tt,pp, &
                                   scf,dscf_b1)
             end if
            else if(mode.eq.2)then
             if(     inf(i)(3:6).eq.'FH  '.or.inf(i)(3:6).eq.'FGH '.or. &
                     inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv')then
              call SCF_DF_DF_DB_DB(iscf,hstat,xf,1.d0,1.d0,xb,1.d0,1.d0,tt,pp, &
                                   scf,dscf_f1,d2scf_f2,d3scf_f2b1,d4scf_f2b2)
              call SCF_DB_DB_DF(   iscf,hstat,xf,1.d0,xb,1.d0,1.d0,tt,pp, &
                                   scf,dscf_b1,d2scf_b2,d3scf_f1b2)
              call SCF_DF_DB(      iscf,hstat,xf,1.d0,xb,1.d0,tt,pp, &
                                   scf,dscf_f1,d2scf_f1b1)
             else if(inf(i)(3:6).eq.'FG  ')then
              call SCF_DB_DB_DF(   iscf,hstat,xf,1.d0,xb,1.d0,1.d0,tt,pp, &
                                   scf,dscf_b1,d2scf_b2,d3scf_f1b2)
              call SCF_DF_DB(      iscf,hstat,xf,1.d0,xb,1.d0,tt,pp, &
                                   scf,dscf_f1,d2scf_f1b1)
             else if(inf(i)(3:6).eq.'F   ')then
              call SCF_DB_DB(      iscf,hstat,xf,xb,1.d0,1.d0,tt,pp, &
                                   scf,dscf_b1,d2scf_b2)
             end if
            else
             stop'unknown mode in PetitR'
            end if

            ! F
            if(mode.eq.0)then
              r(idf,1) = r(idf,1) * scf
            else if(mode.eq.1)then
              if(k.eq.idif)then
                r(idf,1) = r(idf,1) * dscf_b1
              else
                r(idf,1) = r(idf,1) * scf
              end if
            else if(mode.eq.2)then
              if(      k.eq.idif.and.k.eq.jdif )then
                r(idf,1) = r(idf,1) * d2scf_b2
              else if((k.NE.idif.and.k.eq.jdif).or. &
                      (k.eq.idif.and.k.NE.jdif))then
                r(idf,1) = r(idf,1) * dscf_b1
              else if( k.NE.idif.and.k.NE.jdif )then
                r(idf,1) = r(idf,1) * scf
              end if
            end if
            ! G
            if(inf(i)(4:4).eq.'G')then
              do 130 l=1,ndim
                if(mode.eq.0)then
                  if(k.eq.l)then
                    r(idg+l,1) = r(idg+l,1) * dscf_f1
                  else
                    r(idg+l,1) = r(idg+l,1) * scf
                  end if
                else if(mode.eq.1)then
                  if(       k.eq.l.and.k.eq.idif )then
                    r(idg+l,1) = r(idg+l,1) * d2scf_f1b1
                  else if( (k.NE.l.and.k.eq.idif))then
                    r(idg+l,1) = r(idg+l,1) * dscf_b1
                  else if( (k.eq.l.and.k.NE.idif))then
                    r(idg+l,1) = r(idg+l,1) * dscf_f1
                  else if( (k.NE.l.and.k.NE.idif))then
                    r(idg+l,1) = r(idg+l,1) * scf
                  else
                    stop'1-G in PetitR'
                  end if
                else if(mode.eq.2)then
                  if(       k.eq.idif.and.k.eq.jdif.and.k.eq.l )then
                    r(idg+l,1) = r(idg+l,1) * d3scf_f1b2
                  else if( (k.NE.idif.and.k.eq.jdif.and.k.eq.l).or. &
                           (k.eq.idif.and.k.NE.jdif.and.k.eq.l))then
                    r(idg+l,1) = r(idg+l,1) * d2scf_f1b1
                  else if(  k.eq.idif.and.k.eq.jdif.and.k.NE.l )then
                    r(idg+l,1) = r(idg+l,1) * d2scf_b2
                  else if( (k.NE.idif.and.k.eq.jdif.and.k.NE.l).or. &
                           (k.eq.idif.and.k.NE.jdif.and.k.NE.l))then
                    r(idg+l,1) = r(idg+l,1) * dscf_b1
                  else if(  k.NE.idif.and.k.NE.jdif.and.k.eq.l )then
                    r(idg+l,1) = r(idg+l,1) * dscf_f1
                  else if(  k.NE.idif.and.k.NE.jdif.and.k.NE.l )then
                    r(idg+l,1) = r(idg+l,1) * scf
                  else
                    stop'2-G in PetitR'
                  end if
                end if
130           continue
            end if
            ! H
            if(inf(i)(3:6).eq.'FH  '.or.inf(i)(3:6).eq.'FGH ')then
              ict = 0
              do 140 l=1,ndim
              do 140 m=l,ndim
                if(mode_dck.eq.1.and.l.ne.m)go to 140
                ict = ict + 1
                if(mode.eq.0)then
                  if(      k.eq.l.and.k.eq.m )then
                    r(idh+ict,1) = r(idh+ict,1) * d2scf_f2
                  else if((k.NE.l.and.k.eq.m).or. &
                          (k.eq.l.and.k.NE.m))then
                    r(idh+ict,1) = r(idh+ict,1) * dscf_f1
                  else if((k.NE.l.and.k.NE.m))then
                    r(idh+ict,1) = r(idh+ict,1) * scf
                  end if 
                else if(mode.eq.1)then
                  if(      k.eq.l.and.k.eq.m.and.k.eq.idif )then
                    r(idh+ict,1) = r(idh+ict,1) * d3scf_f2b1
                  else if((k.NE.l.and.k.eq.m.and.k.eq.idif).or. &
                          (k.eq.l.and.k.NE.m.and.k.eq.idif))then
                    r(idh+ict,1) = r(idh+ict,1) * d2scf_f1b1
                  else if((k.eq.l.and.k.eq.m.and.k.NE.idif))then
                    r(idh+ict,1) = r(idh+ict,1) * d2scf_f2
                  else if((k.NE.l.and.k.eq.m.and.k.NE.idif).or. &
                          (k.eq.l.and.k.NE.m.and.k.NE.idif))then
                    r(idh+ict,1) = r(idh+ict,1) * dscf_f1
                  else if( k.NE.l.and.k.NE.m.and.k.eq.idif )then
                    r(idh+ict,1) = r(idh+ict,1) * dscf_b1
                  else if( k.NE.l.and.k.NE.m.and.k.NE.idif )then
                    r(idh+ict,1) = r(idh+ict,1) * scf
                  else
                    stop'1-H in PetitR'
                  end if
                else if(mode.eq.2)then
                  if(      k.eq.l.and.k.eq.m.and.k.eq.idif.and.k.eq.jdif )then
                    r(idh+ict,1) = r(idh+ict,1) * d4scf_f2b2
                  else if((k.NE.l.and.k.eq.m.and.k.eq.idif.and.k.eq.jdif).or. &
                          (k.eq.l.and.k.NE.m.and.k.eq.idif.and.k.eq.jdif))then
                    r(idh+ict,1) = r(idh+ict,1) * d3scf_f1b2
                  else if((k.eq.l.and.k.eq.m.and.k.NE.idif.and.k.eq.jdif).or. &
                          (k.eq.l.and.k.eq.m.and.k.eq.idif.and.k.NE.jdif))then
                    r(idh+ict,1) = r(idh+ict,1) * d3scf_f2b1
                  else if( k.NE.l.and.k.NE.m.and.k.eq.idif.and.k.eq.jdif )then
                    r(idh+ict,1) = r(idh+ict,1) * d2scf_f2
                  else if( k.eq.l.and.k.eq.m.and.k.NE.idif.and.k.NE.jdif )then
                    r(idh+ict,1) = r(idh+ict,1) * d2scf_b2
                  else if((k.NE.l.and.k.eq.m.and.k.NE.idif.and.k.eq.jdif).or. &
                          (k.NE.l.and.k.eq.m.and.k.eq.idif.and.k.NE.jdif).or. &
                          (k.eq.l.and.k.NE.m.and.k.NE.idif.and.k.eq.jdif).or. &
                          (k.eq.l.and.k.NE.m.and.k.eq.idif.and.k.NE.jdif))then
                    r(idh+ict,1) = r(idh+ict,1) * d2scf_f1b1
                  else if((k.NE.l.and.k.NE.m.and.k.NE.idif.and.k.eq.jdif).or. &
                          (k.NE.l.and.k.NE.m.and.k.eq.idif.and.k.NE.jdif))then
                    r(idh+ict,1) = r(idh+ict,1) * dscf_b1
                  else if((k.NE.l.and.k.eq.m.and.k.NE.idif.and.k.NE.jdif).or. &
                          (k.eq.l.and.k.NE.m.and.k.NE.idif.and.k.NE.jdif))then
                    r(idh+ict,1) = r(idh+ict,1) * dscf_f1
                  else if((k.NE.l.and.k.NE.m.and.k.NE.idif.and.k.NE.jdif))then
                    r(idh+ict,1) = r(idh+ict,1) * scf
                  else
                    stop'2-H in PetitR'
                  end if
                end if
140           continue
            end if
            ! Hv
            if(inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv')then
              do 150 l=1,ndim
              do 150 m=l,ndim
                if(mode.eq.0)then
                  if(      k.eq.l.and.k.eq.m )then
                    Hfv(l,m) = Hfv(l,m) * d2scf_f2
                  else if((k.NE.l.and.k.eq.m).or. &
                          (k.eq.l.and.k.NE.m))then
                    Hfv(l,m) = Hfv(l,m) * dscf_f1
                  else if( k.NE.l.and.k.NE.m )then
                    Hfv(l,m) = Hfv(l,m) * scf
                  end if
                else if(mode.eq.1)then
                  if(      k.eq.l.and.k.eq.m.and.k.eq.idif )then
                    Hfv(l,m) = Hfv(l,m) * d3scf_f2b1
                  else if((k.NE.l.and.k.eq.m.and.k.eq.idif).or. &
                          (k.eq.l.and.k.NE.m.and.k.eq.idif))then
                    Hfv(l,m) = Hfv(l,m) * d2scf_f1b1
                  else if( k.eq.l.and.k.eq.m.and.k.NE.idif )then
                    Hfv(l,m) = Hfv(l,m) * d2scf_f2
                  else if((k.NE.l.and.k.eq.m.and.k.NE.idif).or. &
                          (k.eq.l.and.k.NE.m.and.k.NE.idif))then
                    Hfv(l,m) = Hfv(l,m) * dscf_f1
                  else if( k.NE.l.and.k.NE.m.and.k.eq.idif )then
                    Hfv(l,m) = Hfv(l,m) * dscf_b1
                  else if( k.NE.l.and.k.NE.m.and.k.NE.idif )then
                    Hfv(l,m) = Hfv(l,m) * scf
                  else
                    stop'1-Hv in PetitR'
                  end if
                else if(mode.eq.2)then
                  if(      k.eq.l.and.k.eq.m.and.k.eq.idif.and.k.eq.jdif )then
                    Hfv(l,m) = Hfv(l,m) * d4scf_f2b2
                  else if((k.NE.l.and.k.eq.m.and.k.eq.idif.and.k.eq.jdif).or. &
                          (k.eq.l.and.k.NE.m.and.k.eq.idif.and.k.eq.jdif))then
                    Hfv(l,m) = Hfv(l,m) * d3scf_f1b2
                  else if((k.eq.l.and.k.eq.m.and.k.NE.idif.and.k.eq.jdif).or. &
                          (k.eq.l.and.k.eq.m.and.k.eq.idif.and.k.NE.jdif))then
                    Hfv(l,m) = Hfv(l,m) * d3scf_f2b1
                  else if( k.NE.l.and.k.NE.m.and.k.eq.idif.and.k.eq.jdif )then
                    Hfv(l,m) = Hfv(l,m) * d2scf_b2
                  else if( k.eq.l.and.k.eq.m.and.k.NE.idif.and.k.NE.jdif )then
                    Hfv(l,m) = Hfv(l,m) * d2scf_f2
                  else if((k.NE.l.and.k.eq.m.and.k.NE.idif.and.k.eq.jdif).or. &
                          (k.NE.l.and.k.eq.m.and.k.eq.idif.and.k.NE.jdif).or. &
                          (k.eq.l.and.k.NE.m.and.k.NE.idif.and.k.eq.jdif).or. &
                          (k.eq.l.and.k.NE.m.and.k.eq.idif.and.k.NE.jdif))then
                    Hfv(l,m) = Hfv(l,m) * d2scf_f1b1
                  else if((k.NE.l.and.k.NE.m.and.k.NE.idif.and.k.eq.jdif).or. &
                          (k.NE.l.and.k.NE.m.and.k.eq.idif.and.k.NE.jdif))then
                    Hfv(l,m) = Hfv(l,m) * dscf_b1
                  else if((k.NE.l.and.k.eq.m.and.k.NE.idif.and.k.NE.jdif).or. &
                          (k.eq.l.and.k.NE.m.and.k.NE.idif.and.k.NE.jdif))then
                    Hfv(l,m) = Hfv(l,m) * dscf_f1
                  else if( k.NE.l.and.k.NE.m.and.k.NE.idif.and.k.NE.jdif )then
                    Hfv(l,m) = Hfv(l,m) * scf
                  else
                    stop'2-Hv in PetitR'
                  end if
                end if
150           continue
            end if
120       continue ! ndim(k) loop

          ! Nugget Effect for Indirect
          if(nsize.eq.rsample)then
            call calc_nugget(rsample,nparent,dxadd,i,1,value) ! 1 is dummy
            r(idf,1) = r(idf,1) * value
          end if
          ! Hessian Vector
          if(inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv')then
            do 200 l=1,ndim
            do 200 m=l,ndim
               H(l,m) = Hfv(l,m)
               H(m,l) = Hfv(l,m)
               V(l,1) = hvec(i,l,1)
200         continue
            HV = matmul(H,V)
            do 210 k=1,ndim
               r(idv+k,1) = HV(k,1)
210         continue
          end if
100     continue ! rsample(i) loop

! Cross Correlation Relaxation Factor
        if(lmax.ne.1)then
          call find_index(1,2,idf,idg,idh,idv, &
                          ndim,nhes,nsize,mode_dck,rsample, &
                          ict_sample,inf)
          do i=idf,nsize
             r(i,1) = r(i,1) * ccrf
          end do
        end if

        end subroutine makeMatrixPetitR
