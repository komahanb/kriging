        subroutine read_set
        use dimKrig
        implicit none
        character(len=10) :: Cdum
        integer :: i,id,ii,jd,ji

        open(10,file='sample.dat',form='formatted',status='old')
        read(10,*)ndim,nsample,nfunc !cheat
        close(10)

        open(10,file='KrigSet.inp',form='formatted',status='old')
        read(10,*)Cdum
        read(10,*)nOPT
        if(nOPT.gt.10)stop'so match nOPT'
        read(10,*)(nfOPT(i),i=1,nOPT)

        read(10,*)Cdum
        read(10,*)nCON
        if(nCON.gt.10)stop'so match nCON'
        do i=1,nCON
        read(10,*)nfCON(i),cCON1(i),vCON(i),cCON2(i)
        end do

        read(10,*)Cdum
        read(10,*)nKRI,(nfKRI(i),i=1,nKRI)
        read(10,*)nDCK,(nfDCK(i),i=1,nDCK)
        read(10,*)nICK,(nfICK(i),i=1,nICK)
        read(10,*)nDMF,(nfDMF(i),i=1,nDMF)
        read(10,*)nVFM,(nfVFM(i),i=1,nVFM)
        if(nKRI.gt.10)stop'so match nKRI'
        if(nDCK.gt.10)stop'so match nDCK'
        if(nICK.gt.10)stop'so match nICK'
        if(nDMF.gt.10)stop'so match nDMF'
        if(nVFM.gt.10)stop'so match nVFM'

        read(10,*)Cdum
        read(10,*)nMXS,nMXG,nMXH,nMXV

        read(10,*)Cdum
        read(10,*)kord
        if(     kord.eq.0)then
           kreg_orig = 1
        else if(kord.lt.0)then
           kreg_orig = abs(kord)*ndim+1
        else if(kord.eq.1)then
           kreg_orig = ndim+1
        else if(kord.eq.2)then
           kreg_orig = (ndim**2-ndim)/2+ndim + ndim+1
        else
           call combination(ndim+kord,ndim,kreg_orig)
        end if

        read(10,*)Cdum
        read(10,*)iscf
        if(iscf.lt.0.or.iscf.gt.5)stop'iscf'

        read(10,*)Cdum
        read(10,*)devmin,devmax
        if(devmin.lt.0.d0)stop'negative devmin'
        if(devmax.lt.devmin)stop'devmax<devmin'

        read(10,*)Cdum
        read(10,*)ndebug
        read(10,*)Cdum
        read(10,*)newflag
        read(10,*)Cdum
        read(10,*)mEI
        read(10,*)Cdum
        read(10,*)mLC
        read(10,*)Cdum
        read(10,*)Reps,ccrf
        Reps = abs(Reps)
        if(mEI.lt.0.or.mEI.gt.4)stop'mode of EI'
        if(mLC.lt.0)stop'mode of Local Correction'
        if(Reps.ge.1.d-5)stop'too big Reps'
        if(ccrf.le.0.99d0.or.ccrf.ge.1.01d0)stop'unreasonable ccrf'
        close(10)

        if(nDCK+nICK.gt.10)then
          stop'so match nCOK'
        else
          nCOK = nDCK+nICK
          id = 1
          ii = 1
          do i=1,nCOK
             jd = nfDCK(id)
             ji = nfICK(ii)
             if(id.gt.nDCK)jd = 100
             if(ii.gt.nICK)ji = 100
             if(jd.lt.ji)then
               nfCOK(i) = jd
               id = id + 1
             else if(jd.gt.ji)then
               nfCOK(i) = ji
               ii = ii + 1
             else
               stop'nfCOK-1'
             end if
          end do
          if(jd.ne.nDCK.and.nDCK.ne.0)stop'nfCOK-2'
          if(ji.ne.nICK.and.nICK.ne.0)stop'nfCOK-3'
        end if

        mode_dck = 0
        if(nDCK.ne.0)then
          open(10,file='direct.inp',form='formatted',status='old')
          read(10,*) Cdum
          if(Cdum.eq.'Diagonal'.or.Cdum.eq.'diagonal')mode_dck = 1
          close(10)
        end if
        if(nICK.ne.0)then
          open(10,file='indirect.inp',form='formatted',status='old')
          read(10,*) itrust,ngput,nhput
          read(10,*) Cnug,Cdum,Vnug
          read(10,*)(tdxinit(i),i=1,2)
          read(10,*) tdxmin,tdxmax
          read(10,*)(tdd_thres(i),i=1,4)
          read(10,*) dpar
          read(10,*) it4opt
          read(10,*) pkl,pku
          read(10,*) vgfac
          read(10,*) bdfac
          read(10,*) mgfac
          close(10)
          if(ngput.lt.-2.or.ngput.gt.-1)stop'ngput'
          if(nhput.lt.-4)stop'nhput'
          if(Cnug.ne.'C'.and.Cnug.ne.'N')stop'Cnug'
          if(Cnug.eq.'C')then
          if(Vnug.lt.0.d0.or.Vnug.gt.1.d0)stop'Vnug'
          else
          if(Vnug.lt.0.d0)stop'Vnug'
          end if
          if(tdxinit(1).lt.tdxmin.or.tdxinit(1).gt.tdxmax)stop'tdxinit1'
          if(tdxinit(2).lt.tdxmin.or.tdxinit(2).gt.tdxmax)stop'tdxinit2'
          if(dpar.lt.0.d0)stop'dpar<0'
          if(pkl.ge.pku)stop'pkl>=pku'
          if(vgfac.lt.1.d0)stop'verygoodfactor<1'
          if(bdfac.ge.1.d0)stop'badfactor>1'
          if(mgfac.lt.1.d0)stop'marginfactor<1'
          ! tdd_thres(1) is final barrier distance to reduce addpts
          ! tdd_thres(2) is a dist ratio per dx to eliminate addpts near realpt
          ! tdd_thres(3) is a dist ratio per dx to eliminate addpts for overlap
          ! tdd_thres(4) is a dist ratio per dx to eliminate addpts for onemore
          do i=1,1
            tdd_thres(i) = dsqrt( dble(ndim)*(tdd_thres(i))**2 )
          end do
        end if

        if(id_proc.eq.0)then
          write(filenum,'(1x,a)')'>> Reading KrigSet.inp'
          write(filenum,'(6x,a,1i4)')'>> # of Optimizing Functions = ',nOPT
          write(filenum,'(6x,a,1i4)')'>> # of Constraint Functions = ',nCON
          write(filenum,'(6x,a,5i4)')'>> # of K/DCOK/ICOK/DUM/VFM  = ', &
                               nKRI,nDCK,nICK,nDMF,nVFM
          write(filenum,'(6x,a,4i4)')'>> # of Max Sample/G/H/Hv    = ', &
                               nMXS,nMXG,nMXH,nMXV
          if(     kord.eq.0.and.kreg_orig.eq.1)then
          write(filenum,'(6x,a)')'>> Ordinary Kriging'
          else if(kord.lt.0)then
          write(filenum,'(6x,a,i2,a,i4,a)') &
          '>> Universal Kriging (order=',abs(kord), &
          ', nterm=',kreg_orig,') (diag)'
          else
          write(filenum,'(6x,a,i2,a,i4,a)') &
          '>> Universal Kriging (order=',kord,', nterm=',kreg_orig,')'
          end if
          if(iscf.eq.0)then
          write(filenum,'(6x,a)')'>> Use Gaussian SCF'
          else if(iscf.eq.1)then
          write(filenum,'(6x,a)')'>> Use Gaussian SCF with Variable Power'
          else if(iscf.eq.2)then
          write(filenum,'(6x,a)')'>> Use Cubic Spline SCF'
          else if(iscf.eq.3)then
          write(filenum,'(6x,a)')'>> Use Wendland C2 SCF'
          else if(iscf.eq.4)then
          write(filenum,'(6x,a)')'>> Use Wendland C4 SCF'
          else if(iscf.eq.5)then
          write(filenum,'(6x,a)')'>> Use Matern SCF'
          end if
          if(itrust.eq.1.and.nICK.ne.0)then
          write(filenum,'(6x,a)')'>> Use Trust Region Approach'
          end if
          if(nDCK.ne.0)then
          if(mode_dck.eq.1)then
          write(filenum,'(6x,a)')'>> Use Diagonal Hess for Direct Cokriging'
          else
          write(filenum,'(6x,a)')'>> Use All Hess for Direct Cokriging'
          end if
          end if
          if(mLC.eq.1)then
          write(filenum,'(6x,a)')'>> Use 2nd Order Additive Correction'
          else if(mLC.eq.2)then
          write(filenum,'(6x,a)')'>> Use 2nd Order Multiplicative Correction'
          else if(mLC.ne.0)then
          write(filenum,'(6x,a,i3)')'>> Use Combined Correction',mLC
          end if
          if(Reps.ne.0.d0) &
          write(filenum,'(6x,a,e10.3)') &
          '>> Add epsilon on Diagonal Rij, epsilon =  ',Reps
          if(ccrf.ne.1.d0) &
          write(filenum,'(6x,a,f12.7)') &
          '>> Cross Correlation Relaxation Factor  =',ccrf
        end if
        end subroutine read_set
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
