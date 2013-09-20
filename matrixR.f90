        subroutine makeMatrixRij(theta,power)
        use dimKrig
        implicit none
! power is active only for iscf = 0/1 (gaussian)
        double precision, dimension(ndim,tdim) :: theta,power
        integer :: i,j,k,l,m,n,o
        integer :: l1,l2,t
        integer :: ict,ict1,ict2
        integer :: ideb
        integer :: idf,idg,idh,idv,jdf,jdg,jdh,jdv
        double precision :: xf,xb,tt,pp
        double precision :: value,prd,gam
        double precision ::  dscfdf, dscfdb,dscf,scf
        double precision :: d2scfdf,d2scfdb,d2scfdfdb
        double precision :: d3scfdf,d3scfdb
        double precision :: d4scf
        character(len=2) :: cl

        double precision, dimension(ndim,ndim,ndim,ndim) :: Hhv
        double precision, dimension(ndim,ndim,ndim) :: Hgv
        double precision, dimension(ndim,ndim) :: H,Hfv
        double precision, dimension(ndim,1)    :: HV,V

        ideb = 0
        if(nsize.ne.rsample.and.iscf.eq.1)then
          stop'Gaussian Power for direct Cokriging'
        end if

        Rij = 1.d0
        do 100 i=1,rsample
          do l=1,lmax
             write(cl,105)l
105          format(i1,'_')
             if(cl.eq.inf(i)(1:2))then
               l1 = l
               go to 101
             end if
          end do
          stop'unknown level of i in Rij'
101       continue
          do 110 j=1,rsample
            do l=1,lmax
               write(cl,105)l
               if(cl.eq.inf(j)(1:2))then
                 l2 = l
                 go to 111
               end if
            end do
            stop'unknown level of j in Rij'
111         continue

            ! l1,l2 are the fidelity levels of i,j
            ! which theta/power
            ict = 0
            do 120 l=1,lmax
            do 120 m=l,lmax
               ict = ict + 1
               if((l1.eq.l.and.l2.eq.m).or.(l1.eq.m.and.l2.eq.l))then
                 t = ict
                 go to 121
               end if
120         continue
            stop'unknown tdim in Rij'
121         continue

            ! which index on Rij
!           call find_index(i,l1,idf,idg,idh,idv, &
!                           ndim,nhes,nsize,mode_dck,rsample, &
!                           ict_sample,inf)
!           call find_index(j,l2,jdf,jdg,jdh,jdv, &
!                           ndim,nhes,nsize,mode_dck,rsample, &
!                           ict_sample,inf)
            idf = iRij(i,1)
            idg = iRij(i,2)
            idh = iRij(i,3)
            idv = iRij(i,4)
            jdf = iRij(j,1)
            jdg = iRij(j,2)
            jdh = iRij(j,3)
            jdv = iRij(j,4)

            if( inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv'.or. &
                inf(j)(3:6).eq.'FHv '.or.inf(j)(3:6).eq.'FGHv' ) then
              Hfv = 1.d0
              Hgv = 1.d0
              Hhv = 1.d0
            end if
            do 150 k=1,ndim !---------------------------------ndim loop
             xf = sampl(i,k)
             xb = sampl(j,k)
             tt = theta(k,t)
             pp = power(k,t)
             if(iscf.eq.0.or.iscf.eq.1)then
               if(pp.gt.0.9999d0.and.pp.lt.1.d0)pp = 1.d0
               if(pp.lt.2.0001d0.and.pp.gt.2.d0)pp = 2.d0
               if(pp.lt.1.d0.or.pp.gt.2.d0)stop'power in Rij'
             end if

             ! AD
             if(      (inf(i)(3:6).eq.'F   ').and. &
                      (inf(j)(3:6).eq.'F   ') ) then
               call SCF_OUT(        iscf,xf,xb,tt,pp,scf)
             else if( (inf(i)(3:6).eq.'FH  '.or.inf(i)(3:6).eq.'FGH '.or. &
                       inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv').and. &
                      (inf(j)(3:6).eq.'FH  '.or.inf(j)(3:6).eq.'FGH '.or. &
                       inf(j)(3:6).eq.'FHv '.or.inf(j)(3:6).eq.'FGHv') ) then
               call SCF_DF_DF_DB_DB(iscf,hstat,xf,1.d0,1.d0,xb,1.d0,1.d0,tt,pp, &
                    scf,dscfdf,d2scfdf,d3scfdf,d4scf)
               call SCF_DB_DB_DF(   iscf,hstat,xf,1.d0,xb,1.d0,1.d0,tt,pp, &
                    scf,dscfdb,d2scfdb,d3scfdb)
               call SCF_DF_DB(      iscf,hstat,xf,1.d0,xb,1.d0,tt,pp, &
                    scf,dscf  ,d2scfdfdb)
             else if( (inf(i)(3:6).eq.'FG  ').and. &
                      (inf(j)(3:6).eq.'FH  '.or.inf(j)(3:6).eq.'FGH '.or. &
                       inf(j)(3:6).eq.'FHv '.or.inf(j)(3:6).eq.'FGHv') ) then
               call SCF_DB_DB_DF(   iscf,hstat,xf,1.d0,xb,1.d0,1.d0,tt,pp, &
                    scf,dscfdb,d2scfdb,d3scfdb)
               call SCF_DF_DB(      iscf,hstat,xf,1.d0,xb,1.d0,tt,pp, &
                    scf,dscfdf,d2scfdfdb)
             else if( (inf(i)(3:6).eq.'FH  '.or.inf(i)(3:6).eq.'FGH '.or. &
                       inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv').and. &
                      (inf(j)(3:6).eq.'FG  ') ) then
               call SCF_DF_DF_DB(   iscf,hstat,xf,1.d0,1.d0,xb,1.d0,tt,pp, &
                    scf,dscfdf,d2scfdf,d3scfdf)
               call SCF_DB_DF(      iscf,hstat,xf,1.d0,xb,1.d0,tt,pp, &
                    scf,dscfdb,d2scfdfdb)
             else if( (inf(i)(3:6).eq.'FG  ').and. &
                      (inf(j)(3:6).eq.'FG  ') ) then
               call SCF_DF_DB(      iscf,hstat,xf,1.d0,xb,1.d0,tt,pp, &
                    scf,dscfdf,d2scfdfdb)
               call SCF_DB(         iscf,hstat,xf,xb,1.d0,tt,pp, &
                    scf,dscfdb)
             else if( (inf(i)(3:6).eq.'F   ').and. &
                      (inf(j)(3:6).eq.'FH  '.or.inf(j)(3:6).eq.'FGH '.or. &
                       inf(j)(3:6).eq.'FHv '.or.inf(j)(3:6).eq.'FGHv') ) then
               call SCF_DB_DB(      iscf,hstat,xf,xb,1.d0,1.d0,tt,pp, &
                    scf,dscfdb,d2scfdb)
             else if( (inf(i)(3:6).eq.'FH  '.or.inf(i)(3:6).eq.'FGH '.or. &
                       inf(i)(3:6).eq.'FHv '.or.inf(i)(3:6).eq.'FGHv').and. &
                      (inf(j)(3:6).eq.'F   ') ) then
               call SCF_DF_DF(      iscf,hstat,xf,1.d0,1.d0,xb,tt,pp, &
                    scf,dscfdf,d2scfdf)
             else if( (inf(i)(3:6).eq.'F   ').and. &
                      (inf(j)(3:6).eq.'FG  ') ) then
               call SCF_DB(         iscf,hstat,xf,xb,1.d0,tt,pp, &
                    scf,dscfdb)
             else if( (inf(i)(3:6).eq.'FG  ').and. &
                      (inf(j)(3:6).eq.'F   ') ) then
               call SCF_DF(         iscf,hstat,xf,1.d0,xb,tt,pp, &
                    scf,dscfdf)
             else
               write(*,*)inf(i),inf(j)
               stop'unknown set for AD'
             end if


             if((l1.lt.l2).or.(l1.eq.l2.and.i.le.j))then
               if(ideb.eq.1.or.ideb.eq.10) &
               write(*,'(a,1x,2a,2i5)')'F-F',inf(i),inf(j),idf,jdf
               Rij(idf,jdf) = Rij(idf,jdf) * scf                    ! F-F
             end if
             if( (inf(i)(3:6).eq.'F   ').and.(inf(j)(3:6).eq.'F   ') )go to 150

             if((inf(j)(4:4).eq.'G'))then                           ! F-G
              if(l1.le.l2)then
               do 200 l=1,ndim
                 if(ideb.eq.2.or.ideb.eq.10) &
                 write(*,'(a,1x,2a,2i5)')'F-G',inf(i),inf(j),idf,jdg+l
                 if(k.eq.l)then
                   Rij(idf,jdg+l) = Rij(idf,jdg+l) * dscfdb
                 else
                   Rij(idf,jdg+l) = Rij(idf,jdg+l) * scf
                 end if
200            continue
              else
               do 205 l=1,ndim
                 if(ideb.eq.2.or.ideb.eq.10) &
                 write(*,'(a,1x,2a,2i5)')'G-F',inf(i),inf(j),jdg+l,idf
                 if(k.eq.l)then
                   Rij(jdg+l,idf) = Rij(jdg+l,idf) * dscfdb
                 else
                   Rij(jdg+l,idf) = Rij(jdg+l,idf) * scf
                 end if
205            continue
              end if
             end if

             if((inf(j)(3:6).eq.'FH  '.or.inf(j)(3:6).eq.'FGH '))then ! F-H
              if(l1.le.l2)then
               ict = 0
               do 210 l=1,ndim
               do 210 m=l,ndim
                  if(mode_dck.eq.1.and.l.ne.m)go to 210
                  ict = ict + 1
                  if(ideb.eq.3.or.ideb.eq.10) &
                  write(*,'(a,1x,2a,2i5,f15.8)')'F-H', &
                  inf(i),inf(j),idf,jdh+ict,d2scfdb
                  if(       k.eq.l.and.k.eq.m  )then
                    Rij(idf,jdh+ict) = Rij(idf,jdh+ict) * d2scfdb
                  else if( (k.eq.l.and.k.NE.m).or. &
                           (k.NE.l.and.k.eq.m) )then
                    Rij(idf,jdh+ict) = Rij(idf,jdh+ict) * dscfdb
                  else if(  k.NE.l.and.k.NE.m  )then
                    Rij(idf,jdh+ict) = Rij(idf,jdh+ict) * scf
                  else
                    stop'Rij F-H'
                  end if
210            continue
              else
               ict = 0
               do 215 l=1,ndim
               do 215 m=l,ndim
                  if(mode_dck.eq.1.and.l.ne.m)go to 215
                  ict = ict + 1
                  if(ideb.eq.3.or.ideb.eq.10) &
                  write(*,'(a,1x,2a,2i5,f15.8)')'H-F', &
                  inf(i),inf(j),jdh+ict,idf,d2scfdb
                  if(       k.eq.l.and.k.eq.m  )then
                    Rij(jdh+ict,idf) = Rij(jdh+ict,idf) * d2scfdb
                  else if( (k.eq.l.and.k.NE.m).or. &
                           (k.NE.l.and.k.eq.m) )then
                    Rij(jdh+ict,idf) = Rij(jdh+ict,idf) * dscfdb
                  else if(  k.NE.l.and.k.NE.m  )then
                    Rij(jdh+ict,idf) = Rij(jdh+ict,idf) * scf
                  else
                    stop'Rij H-F'
                  end if
215            continue
              end if
             end if

             if((inf(j)(3:6).eq.'FHv '.or.inf(j)(3:6).eq.'FGHv'))then ! F-Hv
               do 220 l=1,ndim
               do 220 m=l,ndim
                  if(       k.eq.l.and.k.eq.m  )then
                    Hfv(l,m) = Hfv(l,m) * d2scfdb
                  else if( (k.eq.l.and.k.NE.m).or. &
                           (k.NE.l.and.k.eq.m) )then
                    Hfv(l,m) = Hfv(l,m) * dscfdb
                  else if(  k.NE.l.and.k.NE.m  )then
                    Hfv(l,m) = Hfv(l,m) * scf
                  else
                    stop'Rij F-Hv'
                  end if
220            continue
             end if

             if(inf(i)(4:4).eq.'G'.and.inf(j)(4:4).eq.'G')then          ! G-G
              if((l1.lt.l2).or.(l1.eq.l2.and.i.le.j))then
                if(ideb.eq.5.or.ideb.eq.10) &
                write(*,'(a,1x,2a,2i5)')'G-G',inf(i),inf(j),idg+l,jdg+m
                do 230 l=1,ndim
                do 230 m=1,ndim
                   if(     k.eq.l.and.k.eq.m)then
                     Rij(idg+l,jdg+m) = Rij(idg+l,jdg+m) * d2scfdfdb
                   else if(k.eq.l.and.k.NE.m)then
                     Rij(idg+l,jdg+m) = Rij(idg+l,jdg+m) * dscfdf
                   else if(k.NE.l.and.k.eq.m)then
                     Rij(idg+l,jdg+m) = Rij(idg+l,jdg+m) * dscfdb
                   else if(k.NE.l.and.k.NE.m)then
                     Rij(idg+l,jdg+m) = Rij(idg+l,jdg+m) * scf
                   else
                     stop'Rij G-G'
                   end if
230             continue
              end if
             end if

             if( (inf(i)(4:4).eq.'G').and. &
                 (inf(j)(3:6).eq.'FGH '.or.inf(j)(3:6).eq.'FH  ') )then ! G-H
              if(l1.le.l2)then
               do 240 l=1,ndim
                 ict = 0
                 do 241 m=1,ndim
                 do 241 n=m,ndim
                   if(mode_dck.eq.1.and.m.ne.n)go to 241
                   ict = ict + 1
                   if(ideb.eq.6.or.ideb.eq.10) &
                   write(*,'(a,1x,2a,5i5)')'G-H',inf(i),inf(j),idg+l,jdh+ict
                   if(       k.eq.l.and.k.eq.m.and.k.eq.n  )then
                     Rij(idg+l,jdh+ict) = Rij(idg+l,jdh+ict) * d3scfdb
                   else if(  k.NE.l.and.k.eq.m.and.k.eq.n  )then
                     Rij(idg+l,jdh+ict) = Rij(idg+l,jdh+ict) * d2scfdb
                   else if( (k.eq.l.and.k.eq.m.and.k.NE.n).or. &
                            (k.eq.l.and.k.NE.m.and.k.eq.n) )then
                     Rij(idg+l,jdh+ict) = Rij(idg+l,jdh+ict) * d2scfdfdb
                   else if( (k.NE.l.and.k.eq.m.and.k.NE.n).or. &
                            (k.NE.l.and.k.NE.m.and.k.eq.n) )then
                     Rij(idg+l,jdh+ict) = Rij(idg+l,jdh+ict) * dscfdb
                   else if(  k.eq.l.and.k.NE.m.and.k.NE.n  )then
                     Rij(idg+l,jdh+ict) = Rij(idg+l,jdh+ict) * dscfdf
                   else if(  k.NE.l.and.k.NE.m.and.k.NE.n  )then
                     Rij(idg+l,jdh+ict) = Rij(idg+l,jdh+ict) * scf
                   else
                     stop'Rij G-H'
                   end if
241              continue
240            continue
              else
               do 245 l=1,ndim
                 ict = 0
                 do 246 m=1,ndim
                 do 246 n=m,ndim
                   if(mode_dck.eq.1.and.m.ne.n)go to 246
                   ict = ict + 1
                   if(ideb.eq.6.or.ideb.eq.10) &
                   write(*,'(a,1x,2a,2i5)')'H-G',inf(i),inf(j),jdh+ict,idg+l
                   if(       k.eq.l.and.k.eq.m.and.k.eq.n  )then
                     Rij(jdh+ict,idg+l) = Rij(jdh+ict,idg+l) * d3scfdb
                   else if(  k.NE.l.and.k.eq.m.and.k.eq.n  )then
                     Rij(jdh+ict,idg+l) = Rij(jdh+ict,idg+l) * d2scfdb
                   else if( (k.eq.l.and.k.eq.m.and.k.NE.n).or. &
                            (k.eq.l.and.k.NE.m.and.k.eq.n) )then
                     Rij(jdh+ict,idg+l) = Rij(jdh+ict,idg+l) * d2scfdfdb
                   else if( (k.NE.l.and.k.eq.m.and.k.NE.n).or. &
                            (k.NE.l.and.k.NE.m.and.k.eq.n) )then
                     Rij(jdh+ict,idg+l) = Rij(jdh+ict,idg+l) * dscfdb
                   else if(  k.eq.l.and.k.NE.m.and.k.NE.n  )then
                     Rij(jdh+ict,idg+l) = Rij(jdh+ict,idg+l) * dscfdf
                   else if(  k.NE.l.and.k.NE.m.and.k.NE.n  )then
                     Rij(jdh+ict,idg+l) = Rij(jdh+ict,idg+l) * scf
                   else
                     stop'Rij H-G'
                   end if
246              continue
245            continue
              end if
             end if

             if( (inf(i)(4:4).eq.'G').and. &
                 (inf(j)(3:6).eq.'FGHv'.or.inf(j)(3:6).eq.'FHv ') )then ! G-Hv
               do 250 l=1,ndim
                 do 251 m=1,ndim
                 do 251 n=m,ndim
                   if(       k.eq.l.and.k.eq.m.and.k.eq.n  )then
                     Hgv(l,m,n) = Hgv(l,m,n) * d3scfdb
                   else if(  k.NE.l.and.k.eq.m.and.k.eq.n  )then
                     Hgv(l,m,n) = Hgv(l,m,n) * d2scfdb
                   else if( (k.eq.l.and.k.eq.m.and.k.NE.n).or. &
                            (k.eq.l.and.k.NE.m.and.k.eq.n) )then
                     Hgv(l,m,n) = Hgv(l,m,n) * d2scfdfdb
                   else if( (k.NE.l.and.k.eq.m.and.k.NE.n).or. &
                            (k.NE.l.and.k.NE.m.and.k.eq.n) )then
                     Hgv(l,m,n) = Hgv(l,m,n) * dscfdb
                   else if(  k.eq.l.and.k.NE.m.and.k.NE.n  )then
                     Hgv(l,m,n) = Hgv(l,m,n) * dscfdf
                   else if(  k.NE.l.and.k.NE.m.and.k.NE.n  )then
                     Hgv(l,m,n) = Hgv(l,m,n) * scf
                   else
                     stop'Rij G-Hv'
                   end if
251              continue
250            continue
             end if

             if( (inf(i)(3:6).eq.'FGH '.or.inf(i)(3:6).eq.'FH  ').and. &
                 (inf(j)(3:6).eq.'FGH '.or.inf(j)(3:6).eq.'FH  ') )then ! H-H
              if((l1.lt.l2).or.(l1.eq.l2.and.i.le.j))then
                ict1 = 0
                do 260 l=1,ndim
                do 260 m=l,ndim
                  if(mode_dck.eq.1.and.l.ne.m)go to 260
                  ict1 = ict1 + 1
                  ict2 = 0
                  do 261 n=1,ndim
                  do 261 o=n,ndim
                    if(mode_dck.eq.1.and.n.ne.o)go to 261
                    ict2 = ict2 + 1
                    if(ideb.eq.8.or.ideb.eq.10) &
                    write(*,'(a,1x,2a,2i5,f15.8)')'H-H', &
                    inf(i),inf(j),idh+ict1,jdh+ict2,d4scf
                    if(       k.eq.l.and.k.eq.m.and.k.eq.n.and.k.eq.o  )then
                      Rij(idh+ict1,jdh+ict2) = Rij(idh+ict1,jdh+ict2)*d4scf
                    else if( (k.NE.l.and.k.eq.m.and.k.eq.n.and.k.eq.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.eq.n.and.k.eq.o) )then
                      Rij(idh+ict1,jdh+ict2) = Rij(idh+ict1,jdh+ict2)*d3scfdb
                    else if( (k.eq.l.and.k.eq.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.eq.l.and.k.eq.m.and.k.eq.n.and.k.NE.o) )then
                      Rij(idh+ict1,jdh+ict2) = Rij(idh+ict1,jdh+ict2)*d3scfdf
                    else if(  k.NE.l.and.k.NE.m.and.k.eq.n.and.k.eq.o  )then
                      Rij(idh+ict1,jdh+ict2) = Rij(idh+ict1,jdh+ict2)*d2scfdb
                    else if(  k.eq.l.and.k.eq.m.and.k.NE.n.and.k.NE.o  )then
                      Rij(idh+ict1,jdh+ict2) = Rij(idh+ict1,jdh+ict2)*d2scfdf
                    else if( (k.NE.l.and.k.eq.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.NE.l.and.k.eq.m.and.k.eq.n.and.k.NE.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.eq.n.and.k.NE.o) )then
                      Rij(idh+ict1,jdh+ict2) = Rij(idh+ict1,jdh+ict2)*d2scfdfdb
                    else if( (k.NE.l.and.k.NE.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.NE.l.and.k.NE.m.and.k.eq.n.and.k.NE.o) )then
                      Rij(idh+ict1,jdh+ict2) = Rij(idh+ict1,jdh+ict2)*dscfdb
                    else if( (k.NE.l.and.k.eq.m.and.k.NE.n.and.k.NE.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.NE.n.and.k.NE.o) )then
                      Rij(idh+ict1,jdh+ict2) = Rij(idh+ict1,jdh+ict2)*dscfdf
                    else if(  k.NE.l.and.k.NE.m.and.k.NE.n.and.k.NE.o  )then
                      Rij(idh+ict1,jdh+ict2) = Rij(idh+ict1,jdh+ict2)*scf
                    else
                      stop'Rij H-H'
                    end if
261               continue
260             continue
              end if
             end if
             if( (inf(i)(3:6).eq.'FGH '.or.inf(i)(3:6).eq.'FH  ').and. &
                 (inf(j)(3:6).eq.'FGHv'.or.inf(j)(3:6).eq.'FHv ') ) then ! H-Hv
               do 270 l=1,ndim
               do 270 m=l,ndim
                 do 271 n=1,ndim
                 do 271 o=n,ndim
                    if(       k.eq.l.and.k.eq.m.and.k.eq.n.and.k.eq.o  )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d4scf
                    else if( (k.NE.l.and.k.eq.m.and.k.eq.n.and.k.eq.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.eq.n.and.k.eq.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d3scfdb
                    else if( (k.eq.l.and.k.eq.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.eq.l.and.k.eq.m.and.k.eq.n.and.k.NE.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d3scfdf
                    else if(  k.NE.l.and.k.NE.m.and.k.eq.n.and.k.eq.o  )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d2scfdb
                    else if(  k.eq.l.and.k.eq.m.and.k.NE.n.and.k.NE.o  )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d2scfdf
                    else if( (k.NE.l.and.k.eq.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.NE.l.and.k.eq.m.and.k.eq.n.and.k.NE.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.eq.n.and.k.NE.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d2scfdfdb
                    else if( (k.NE.l.and.k.NE.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.NE.l.and.k.NE.m.and.k.eq.n.and.k.NE.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * dscfdb
                    else if( (k.NE.l.and.k.eq.m.and.k.NE.n.and.k.NE.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.NE.n.and.k.NE.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * dscfdf
                    else if(  k.NE.l.and.k.NE.m.and.k.NE.n.and.k.NE.o  )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * scf
                    else
                      stop'Rij H-Hv'
                    end if
271              continue
270            continue
             end if
             if( (inf(i)(3:6).eq.'FGHv'.or.inf(i)(3:6).eq.'FHv ').and. &
                 (inf(j)(3:6).eq.'FGHv'.or.inf(j)(3:6).eq.'FHv ') ) then ! Hv-Hv
              if((l1.lt.l2).or.(l1.eq.l2.and.i.le.j))then
               do 280 l=1,ndim
               do 280 m=1,ndim
                 do 281 n=1,ndim
                 do 281 o=1,ndim
                    if(       k.eq.l.and.k.eq.m.and.k.eq.n.and.k.eq.o  )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d4scf
                    else if( (k.NE.l.and.k.eq.m.and.k.eq.n.and.k.eq.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.eq.n.and.k.eq.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d3scfdb
                    else if( (k.eq.l.and.k.eq.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.eq.l.and.k.eq.m.and.k.eq.n.and.k.NE.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d3scfdf
                    else if(  k.NE.l.and.k.NE.m.and.k.eq.n.and.k.eq.o  )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d2scfdb
                    else if(  k.eq.l.and.k.eq.m.and.k.NE.n.and.k.NE.o  )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d2scfdf
                    else if( (k.NE.l.and.k.eq.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.NE.l.and.k.eq.m.and.k.eq.n.and.k.NE.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.eq.n.and.k.NE.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * d2scfdfdb
                    else if( (k.NE.l.and.k.NE.m.and.k.NE.n.and.k.eq.o).or. &
                             (k.NE.l.and.k.NE.m.and.k.eq.n.and.k.NE.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * dscfdb
                    else if( (k.NE.l.and.k.eq.m.and.k.NE.n.and.k.NE.o).or. &
                             (k.eq.l.and.k.NE.m.and.k.NE.n.and.k.NE.o) )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * dscfdf
                    else if(  k.NE.l.and.k.NE.m.and.k.NE.n.and.k.NE.o  )then
                      Hhv(l,m,n,o) = Hhv(l,m,n,o) * scf
                    else
                      stop'Rij H-Hv'
                    end if
281              continue
280            continue
              end if
             end if
150         continue !--------------------------------------- ndim(k) loop

            ! Nugget Effect for Indirect
            if(nsize.eq.rsample)then
              call calc_nugget(rsample,nparent,dxadd,i,j,value)
              Rij(idf,jdf) = Rij(idf,jdf) * value
            end if

            ! Combination of Hessian Vector
            if(inf(j)(3:6).eq.'FHv '.or.inf(j)(3:6).eq.'FGHv')then ! F-Hv
              do 310 l=1,ndim
               do 311 m=l,ndim
                H(l,m) = Hfv(l,m)
                H(m,l) = Hfv(l,m)
                V(l,1) = hvec(j,l,1)
311            continue
310           continue
              HV = matmul(H,V)
              if(l1.le.l2)then
                do 312 k=1,ndim
                   Rij(idf,jdv+k) = HV(k,1)
                   if(ideb.eq.4.or.ideb.eq.10) &
                   write(*,'(a,1x,2a,2i5,f15.8)')'F-V', &
                   inf(i),inf(j),idf,jdv+k
312             continue
              else
                do 313 k=1,ndim
                   Rij(jdv+k,idf) = HV(k,1)
                   if(ideb.eq.4.or.ideb.eq.10) &
                   write(*,'(a,1x,2a,2i5,f15.8)')'V-F', &
                   inf(i),inf(j),jdv+k,idf
313             continue
              end if
            end if
            if( (inf(i)(4:4).eq.'G').and. &
                (inf(j)(3:6).eq.'FGHv'.or.inf(j)(3:6).eq.'FHv ') )then ! G-Hv
              do 320 l=1,ndim
               do 321 m=1,ndim
               do 321 n=m,ndim
                 H(m,n) = Hgv(l,m,n)
                 H(n,m) = Hgv(l,m,n)
                 V(m,1) = hvec(j,m,1)
321            continue
               HV = matmul(H,V)
               if(l1.le.l2)then
                 do 322 k=1,ndim
                   Rij(idg+l,jdv+k) = HV(k,1)
                   if(ideb.eq.7.or.ideb.eq.10) &
                   write(*,'(a,1x,2a,2i5,f15.8)')'G-V', &
                   inf(i),inf(j),idg+l,jdv+k
322              continue
               else
                 do 323 k=1,ndim
                   Rij(jdv+k,idg+l) = HV(k,1)
                   if(ideb.eq.7.or.ideb.eq.10) &
                   write(*,'(a,1x,2a,2i5,f15.8)')'V-G', &
                   inf(i),inf(j),jdv+k,idg+l
323              continue
               end if
320           continue
            end if
            if( (inf(i)(3:6).eq.'FGH '.or.inf(i)(3:6).eq.'FH  ').and. &
                (inf(j)(3:6).eq.'FGHv'.or.inf(j)(3:6).eq.'FHv ') ) then ! H-Hv
              ict = 0
              do 330 l=1,ndim
              do 330 m=l,ndim
               if(mode_dck.eq.1.and.l.ne.m)go to 330
               ict = ict + 1
               do 331 n=1,ndim
               do 331 o=n,ndim
                 H(n,o) = Hhv(l,m,n,o)
                 H(o,n) = Hhv(l,m,n,o)
                 V(n,1) = hvec(j,n,1)
331            continue
               HV = matmul(H,V)
               if(l1.le.l2)then
                 do 332 k=1,ndim
                   Rij(idh+ict,jdv+k) = HV(k,1)
                   if(ideb.eq.9.or.ideb.eq.10) &
                   write(*,'(a,1x,2a,2i5,f15.8)')'H-V', &
                   inf(i),inf(j),idh+ict,jdv+k
332              continue
               else
                 do 333 k=1,ndim
                   Rij(jdv+k,idh+ict) = HV(k,1)
                   if(ideb.eq.9.or.ideb.eq.10) &
                   write(*,'(a,1x,2a,2i5,f15.8)')'V-H', &
                   inf(i),inf(j),jdv+k,idh+ict
333              continue
               end if
330           continue
            end if
            if( (inf(i)(3:6).eq.'FGHv'.or.inf(i)(3:6).eq.'FHv ').and. &
                (inf(j)(3:6).eq.'FGHv'.or.inf(j)(3:6).eq.'FHv ') ) then ! Hv-Hv
              do 340 l=1,ndim
              do 340 n=1,ndim
               prd = 0.d0
               do 341 m=1,ndim
               do 341 o=1,ndim
                 prd = prd + Hhv(l,m,n,o)*hvec(i,m,1)*hvec(j,o,1)
341            continue
               if((l1.lt.l2).or.(l1.eq.l2.and.i.le.j))then
                 Rij(idv+l,jdv+n) = prd
                 if(ideb.eq.9.or.ideb.eq.10) &
                 write(*,'(a,1x,2a,2i5,f15.8)')'V-V', &
                 inf(i),inf(j),idv+l,jdv+n
               end if
340           continue
            end if

110       continue ! j loop
100     continue ! i loop

! Cross Correlation Relaxation Factor
        if(lmax.ne.1)then
         do 400 l1=1,lmax-1
          call find_index(1,l1,idf,idg,idh,idv, &
                          ndim,nhes,nsize,mode_dck,rsample, &
                          ict_sample,inf)
          l2=l1+1
          call find_index(1,l2,jdf,jdg,jdh,jdv, &
                          ndim,nhes,nsize,mode_dck,rsample, &
                          ict_sample,inf)
          do 410 i=idf,jdf-1
          do 410 j=jdf,nsize
             Rij(i,j) = Rij(i,j) * ccrf
410       continue
400      continue
        end if
! Symmetric Correlation Matrix / Diagonal Terms
        gam = dble(1000+nsize)*Reps
        do 500 i=1,nsize
        do 500 j=i,nsize
          Rij(j,i) = Rij(i,j)
          if(i.eq.j)then
            if(Reps.ne.0.d0)then
              Rij(i,j) = Rij(i,j) + gam
            end if
!           if(i.gt.2)Rij(i,j) = Rij(i,j) * (-1.d0)
          else
!           if(i.le.2.and.j.gt.2)then
!             Rij(j,i) = Rij(i,j)*(-1.d0)
!           end if
          end if
500     continue

        return
999     continue
!       write(*,'(4i5)')id_proc,nsize,rsample,ndim
        do i=1,nsize
          if(id_proc.eq.0) &
          write(*,'(i3,99f9.4)')i,(Rij(i,j),j=1,nsize)
        end do
        return
        call stop_all
        end subroutine makeMatrixRij
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine find_index(inp,linp,idf,idg,idh,idv, &
                              ndim,nhes,nsize,mode_dck,rsample,ict,inf)
        implicit none
! inp  : index of sample
! linp : level of fidelity
        integer, intent(in) :: inp,linp
        integer, intent(in) :: ndim,nhes,nsize,mode_dck,rsample
        integer, dimension(0:10,0:10), intent(in) :: ict
        character(len=6), dimension(rsample), intent(in) :: inf
        integer, intent(out) :: idf,idg,idh,idv

        integer :: i,l
        integer :: idx,nh
        character(len=2) :: cl

        nh = nhes
        if(mode_dck.eq.1)nh = ndim

        ! for higher fidelity levels
        idx = 0
        do 100 l=1,linp-1
           idx = idx + ict(l,1)*(1) &
                     + ict(l,2)*(1+ndim) &
                     + ict(l,4)*(1+ndim+nh) &
                     + ict(l,6)*(1+ndim   +ndim)
100     continue

        ! for target fidelity level
        write(cl,101)linp
101     format(i1,'_')

        idf = idx
        idg = idx + ict(linp,0)
        idh = idx + ict(linp,0)+ndim*(ict(linp,2)+ict(linp,4)+ict(linp,6))
        idv = idx + ict(linp,0)+ndim*(ict(linp,2)+ict(linp,4)+ict(linp,6)) &
                               +nh  *(ict(linp,4))
        do 200 i=1,inp-1
          if(inf(i)(1:2).ne.cl)go to 200
          idf = idf + 1
          if(inf(i)(3:6).ne.'F   ')idg = idg + ndim
          if(inf(i)(3:6).eq.'FGH ')idh = idh + nh
          if(inf(i)(3:6).eq.'FGHv')idv = idv + ndim
200     continue

        idf = idf + 1 !! idf is on the place !!
!       write(*,'(7i3)')inp,linp,idf,idg,idh,idv,nsize
        if(idf.gt.nsize)stop'idf.gt.nsize in find_index'
        if(idg.gt.nsize)stop'idg.gt.nsize in find_index'
        if(idh.gt.nsize)stop'idh.gt.nsize in find_index'
        if(idv.gt.nsize)stop'idv.gt.nsize in find_index'
        end subroutine find_index
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine make_index_Rij
        use dimKrig
        implicit none
        integer :: i,l,l1
        character(len=2) :: cl
        integer :: idf,idg,idh,idv

        iRij = 0
        do 100 i=1,rsample
          do l=1,lmax
             write(cl,105)l
105          format(i1,'_')
             if(cl.eq.inf(i)(1:2))then
               l1 = l
               go to 101
             end if
          end do
          stop'unknown level of i in Rij in make_index'
101       continue

          ! which index on Rij
          call find_index(i,l1,idf,idg,idh,idv, &
                          ndim,nhes,nsize,mode_dck,rsample, &
                          ict_sample,inf)
          iRij(i,1) = idf
          iRij(i,2) = idg
          iRij(i,3) = idh
          iRij(i,4) = idv
100     continue

        end subroutine make_index_Rij
