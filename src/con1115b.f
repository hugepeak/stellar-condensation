C     "con1115B.for"    11/16/94
      subroutine readcontfile(filename)
      IMPLICIT REAL*8(A-H,O-Z)
C $DEFINE mswindows
C $UNDEFINE mswindows
C $IF DEFINED (mswindows)
C      INCLUDE 'FLIB.FD'
C $ENDIF
      INCLUDE 'con0621i.f'
C
      CHARACTER*2  we2
      CHARACTER*80 A80
      CHARACTER*20 filename
      dimension np(islimit),mflag(mlimit)
C $IF DEFINED (mswindows)
      INTEGER*2 IUNIT,IREQ
C      RECORD /QWINFO/ WINFO
C $ENDIF
C
      OPEN (9,FILE=filename)
      idebug=0
      isummary=0
      itable=0
      imode=0
      do i=1,mlimit
         wtb(i)='  '
      enddo
C     read flags:
C     idebug=0 normal; =1 or 2 print more info; =3 print initial values & stop
C     isummary=0 no summary file; =1 make summary file
C     itable=0 no table file; >=1 make table file --- cf. subroutine tables()
      read(9,'(A80)') a80
      read(a80,'(3I2,50A2)')idebug,isummary,itable,(wtb(i),i=1,mlimit)
c      print *,(wtb(i),',',i=1,mlimit)
      i=index(filename,'.')
      j=index(filename,' ')
C     Output file name is "????????.OUT"
      if (i.eq.0) then
         if ((j.eq.0).or.(j.gt.9)) then
            filename(9:12)='.OUT'
         else
            filename(j:j+3)='.OUT'
         endif
      else
         filename(i:i+3)='.OUT'
      endif
      OPEN (8,FILE=filename)
C $IF DEFINED (mswindows)
C         if (iwindow.ge.1) then
C             iunit=17
C             ireq=qwin$framemax
C             ii=getwsizeqq(iunit,ireq,winfo)
C             iwidth=winfo.W
C             iheight=winfo.H
C         endif
C         if (iwindow.ge.1) then
C             open (17,file='USER',title='Output')
C             iunit=17
C             winfo.type=qwin$set
C             winfo.H=7
C             winfo.W=iwidth
C             if (iwidth.gt.86) winfo.W=82
C             winfo.X=4
C             winfo.Y=iheight-17
C             if (iheight-17.lt.3) winfo.Y=3
C             ii=setwsizeqq(iunit,winfo)
C         endif
C         if (iwindow.ge.1) then
C             open (18,file='USER',title='Status')
C             iunit=18
C             winfo.type=qwin$set
C             winfo.H=6
C             winfo.W=iwidth
C             if (iwidth.gt.45) winfo.W=45
C             winfo.X=iwidth-45
C             winfo.Y=0
C             if (iwidth-45.lt.1) winfo.X=0
C             ii=setwsizeqq(iunit,winfo)
C         endif
C $ENDIF
C     Summary file name is "????????.SUM"
      if (isummary.ge.1) then
         if (i.eq.0) then
            if ((j.eq.0).or.(j.gt.9)) then
               filename(9:12)='.SUM'
            else
               filename(j:j+3)='.SUM'
            endif
         else
            filename(i:i+3)='.SUM'
         endif
         open (7,file=filename)
        if (isummary.ge.2) then
         if (i.eq.0) then
            if ((j.eq.0).or.(j.gt.9)) then
               filename(9:12)='.NON'
            else
               filename(j:j+3)='.NON'
            endif
         else
            filename(i:i+3)='.NON'
         endif
         open (15,file=filename)
        endif
      endif
C $IF DEFINED (mswindows)
C         if (iwindow.ge.1) then
C             open (16,file='USER',title='Non-ideal Solutions')
C             iunit=16
C             winfo.type=qwin$set
C             winfo.H=7
C             winfo.W=iwidth-6
C             if (iwidth.gt.90) winfo.W=82
C             winfo.X=8
C             winfo.Y=iheight-9
C             if (iheight-9.lt.6) winfo.Y=6
C             ii=setwsizeqq(iunit,winfo)
C         endif
C $ENDIF
C     write program name, control & data file names, comments, etc.
      write(8,'(A10,A50)') 'Program : ',version
      if (isummary.ge.1) write(7,'(A10,A50)') 'Program : ',version
      read(9,'(A80)')a80
      comcntl=a80
      write(8,'(A15,A12)') 'Control File : ',contfile
      if(isummary.ge.1) write(7,'(A15,A12)') 'Control File : ',contfile
      print *, 'Control File : ',contfile
      write(8,'(A79)') A80
      if(isummary.ge.1) write(7,'(A79)') A80
      write(*,'(1X,A79)') A80
      iord=3
      ibar=0
      read(9,*) datafile,itype,iord,ibar
      write(8,'(A26,A12)') 'Thermodynamic Data File : ',datafile
      if(isummary.ge.1) write(7,'(A26,A12)')
     >'Thermodynamic Data File : ',datafile
      print *, 'Thermodynamic Data File : ',datafile
      if (ibar.eq.1) then
         rgas2=rgas1/1d2
      else
         rgas2=rgas1/1.01325d2
      endif
      filename=datafile
      if (itype.eq.1) call readcoeff(filename)
      if (itype.eq.2) call readdg(filename)
      read(9,*) imode
      if ((imode.lt.0).or.(imode.gt.3)) then
         print *,'Invalid calculation mode !'
         stop
      endif
      nssorg=nss
      do i=1,nsslimit
         nssghost(i)=0
      enddo
C     initial total pressure if only monatomic gases are present
      read(9,*) pstart,pstop,dpres,dpmin
      ptot=pstart
      p0=ptot
C     start temperature, stop temperature, temperature step, resolution
      read(9,*) tstart,tstop,dtemp,dtmin
C     expansion parameter gam=ç/(ç-1), ç=Cp/Cv=<5/3 ("=" for ideal gas)
C     if gam=0 no expansion; and temperature at which expansion starts
C     initial guess for TG required for mode=2 or 3
      read(9,*) gam,ti,tg
c      if (imode.ne.0) gam=0d0
      if (tg.le.0d0) tg=1.8d0
      t0=tstart
      ti=ti/1d2
C     read # of elements
      read(9,*)mtemp
      do k=1,m
         r(k)=0d0
         mflag(k)=0
      enddo
C     read names of elements, abundances, initial guesses of pressures
      do i=1,mtemp
         read(9,*) we2,rtemp,ptemp,awtemp
C         print *, mtemp
C        pause
         do k=1,m
            if (we(k).eq.we2) then
               if (mflag(k).ne.0) then
                  print *,' --- More than two abundances are given ',
     >                    'to one element !'
                  stop
               endif
               mflag(k)=1
               r(k)=rtemp
               p(k)=ptemp
                print *, r(k)
               aw(k)=awtemp
            endif
         enddo
      enddo
C                print *, r(21)
C                print *, r(22)
C                pause

C       OPEN (99,FILE='super.abu')
      
C      do i=1, m
C      read (99,*) pstore(i)
C      print *, pstore(i)
C      enddo
C      close(99)
C      pause
C        do k=1,m
C        r(k)=pstore(k)
C        print *, r(k)
C        enddo
C        pause
      iflag=0
      do k=1,m
         if (mflag(k).eq.0) iflag=1
      enddo
      if (iflag.eq.1) then
         print *,' Warning: Abundances of some elements are not ',
     >           'given; 0.0 assumed'
      endif
C     abundances are normalized to their sum=1
      R0=0d0
      DO 12 I=1,M
 12   R0=R0+R(I)
      DO 14 I=1,M
 14   R(I)=R(I)/R0
      do k=1,m
         if (p(k).le.0d0) then
C           if initial guess is not given (=<0)
            p(k)=r(k)*1d-5
         else
            p(k)=p(k)/ptot
         endif
      enddo
      IC=0
      ID=0
      ICONT=0
      if (itable.ne.0) then
         ntable=0
         ist=0
         do j=1,ntlimit
            do i=1,istlimit
               pstable(i,j)=0d0
            enddo
         enddo
         do j=1,nslimit
            do i=1,m
               itbflag(i,j)=0
            enddo
         enddo
      endif
C     reduce # of elements if no abundances are given for some elements
C     also delete species which contain the elements with no abundances
      mtemp=m
      do i=mtemp,1,-1
      if (r(i).eq.0d0)  then
         nstemp=ns+msExtra
         do j=nstemp,1,-1
            if(nu(j,i).gt.0) then
               if (j.gt.ngp) then
                  print *, ' No abundances are given for elements ',
     >            'used in solid solutions.'
c                  stop
               endif
               if (j.le.n) n=n-1
               if (j.le.ngp) ngp=ngp-1
               ns=ns-1
               do k=1,nss
                  nssp(k)=nssp(k)-1
               enddo
               if (j-1.lt.ns+msExtra) then
                  do k=j,ns+msExtra
                     nspecies(k)=nspecies(k+1)
                     wc(k)=wc(k+1)
                     wref(k)=wref(k+1)
                     do l=1,m
                       nu(k,l)=nu(k+1,l)
                       dnu(k,l)=dnu(k+1,l)
                     enddo
                     do l=1,iord
                       xc(l,k)=xc(l,k+1)
                     enddo
                  enddo
               endif
            endif
         enddo
         m=m-1
c         do j=1,nss
c            if (mss(j).ne.0) mss(j)=mss(j)-1
c         enddo
         if (i-1.lt.m) then
            do j=i,m
               r(j)=r(j+1)
               p(j)=p(j+1)
               aw(j)=aw(j+1)
               we(j)=we(j+1)
               wee(j)=wee(j+1)
               do k=1,ns+msExtra
                  nu(k,j)=nu(k,j+1)
                  dnu(k,j)=dnu(k,j+1)
               enddo
            enddo
         endif
      endif
      enddo
      call checkposition()
      do I=1,M
      DLP(I)=DLOG10(P(I))
      enddo
      if(nss.ne.0) then
         do k=1,nsslimit
            mss(k)=0
         enddo
      endif
      N1=N+1
C     read species which already exist as condensates at the starting temp.
      read(9,*)is
      NIS=N+IS
      IF(IS.EQ.0) GO TO 7
      if(is.gt.islimit) then
         print *,' Exceed IS limit !'
         stop
      endif
      do I=1,IS
        J=i+m+n
        nq(i+n)=0
        READ(9,*) NP(I),P(J)
        do k=1,ns+msExtra
           if(np(i).eq.nspecies(k)) nq(i+n)=k
        enddo
        if(nss.ne.0) then
           do k=1,nss
              if (nssp(k).eq.nq(i+n)) mss(k)=i+m
           enddo
        endif
        if (nssspinel.ne.0) then
        if ((i+m.eq.mss(nssspinel)).or.(i+m.eq.mss(nssspinel)+1)) then
           if (p(j).le.0d0) then
              p(j)=1d-5
           endif
           goto 400
        endif
        endif
        if (p(j).le.0d0) then
              FAC=1d0
              DO 65 k=1,M
              IF(NU(nq(i+n),k)) 63,65,63
 63           FAC=DMIN1(R(k)/dnu(nq(i+n),k),FAC)
 65           CONTINUE
              p(j)=fac*1d-5
        else
              P(j)=P(j)*t0*rgas2/ptot
        endif
 400    continue
      enddo
      do i=1,is
         if (nq(i+n).eq.0) then
            print *,' Cannot find the initial condensed species !'
            stop
         endif
      enddo
      if (itable.ne.0) then
         ist=is
         do i=1,ist
            nqt(i)=nq(i+n)
         enddo
      endif
 7    continue
      call readseed()
      close(9)
C     end of reading control file
      return
      end

      subroutine Psave(indx)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
C     save values of variables temporarily
      issave(indx)=is
      tgsave(indx)=tg
      nis=n+is
      do 50 ii=1,m+nis
 50   psaves(indx,ii)=p(ii)
      do 51 ii=1,nis
 51   nqsave(indx,ii)=nq(ii)
      do 52 ii=1,nsslimit
 52   msssave(indx,ii)=mss(ii)
      do 53 ii=1,nss
      do 53 jj=1,ms(ii)
 53   sxsave(indx,ii,jj)=sx(ii,jj)
      xcrsave(indx)=xcr
      xfesave(indx)=xfe
      do 54 ii=1,ns+msExtra-n
         xk0save(indx,ii)=xk0(ii)
         xk1save(indx,ii)=xk1(ii)
         if (indx.eq.1) then
            xk0save(indx,ii)=xk1(ii)
         endif
 54   continue
      return
      end

      subroutine Prestore(indx)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
C     restore values of variables
      is=issave(indx)
      tg=tgsave(indx)
      nis=n+is
      do 50 ii=1,m+nis
 50   p(ii)=psaves(indx,ii)
      DO 42 K=1,M
 42   DLP(K)=DLOG10(P(K))
      do 51 ii=1,nis
 51   nq(ii)=nqsave(indx,ii)
      do 52 ii=1,nsslimit
 52   mss(ii)=msssave(indx,ii)
      do 53 ii=1,nss
      do 53 jj=1,ms(ii)
 53   sx(ii,jj)=sxsave(indx,ii,jj)
      xcr=xcrsave(indx)
      xal=1d0-xcr
      xfe=xfesave(indx)
      xmg=1d0-xcr
      do 54 ii=1,ns+msExtra-n
      xk0(ii)=xk0save(indx,ii)
 54   xk1(ii)=xk1save(indx,ii)
C     total pressure
      Ptotal=0d0
      do i=1,m+n
         Ptotal=Ptotal+P(i)
      enddo
      Ptotal=Ptotal*PTOT
C     total volume for 1 mol of elements
      Vtotal=rgas2*t0/tg/ptot
C     mole fractions
      call calcsx()
C     free energy of the system
      call freeenergy(G)
      return
      end

      subroutine packghost()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      if (nss.gt.nssorg) then
         do i=nss,nssorg+1,-1
            if (mss(i).gt.0) then
               do j=1,i-1
                  if (((nssghost(i).eq.nssghost(j)).or.
     >                (nssghost(i).eq.j)).and.(mss(j).eq.0)) then
                     mss(j)=mss(i)
                     mss(i)=0
                     do jj=1,ms(i)
                        nq(n+mss(j)-m+jj-1)=nssp(j)+jj-1
                        do k=1,nseed(i)
                           sxold(j,jj,k)=sxold(i,jj,k)
                        enddo
                     enddo
                  endif
               enddo
            endif
         enddo
      endif
      return
      end

      subroutine polymorphs()
      IMPLICIT REAL*8(A-H,O-Z)
      dimension ispoly(10),gpoly(10)
      INCLUDE 'con0621i.f'
C     search for polymorphs of condensed pure solids
C     if one of polymorphs is more stable than the condensate, replace it
      icpst=0
      icpdt=0
      if (is.eq.0) return
      do i=1,is
        if (nq(i+n).gt.ngp) goto 20
          icpoly=0
          f1=0d0
          do k=1,m
            if (dnu(nq(n+i),k).eq.0d0) goto 10
            if (f1.ne.0d0) goto 10
            f1=dnu(nq(n+i),k)
 10         continue
          enddo
        do j=n1,ngp
          f2=0d0
          do k=1,m
            if (dnu(j,k).eq.0d0) goto 11
            if (f2.ne.0d0) goto 11
            f2=dnu(j,k)
 11         continue
          enddo
          iflag=0
          do k=1,m
            if (dabs(dnu(nq(n+i),k)/f1-dnu(j,k)/f2).gt.eps) iflag=1
          enddo
          if (iflag.eq.0) then
            icpoly=icpoly+1
            ispoly(icpoly)=j
            gpoly(icpoly)=xc(1,j)
            do ii=2,iord
              gpoly(icpoly)=gpoly(icpoly)+xc(ii,j)/((t0/1d2)**(ii-1))
            enddo
            gpoly(icpoly)=-gpoly(icpoly)*dlog(10d0)*t0*rgas1/1d3/f2
          endif
        enddo
          if (icpoly.eq.0) then
            print *,'Error in Polymorphs subroutine !!!'
            stop
          elseif (icpoly.gt.1) then
            gmin=gpoly(1)
            iimin=1
            do ii=2,icpoly
              if (gpoly(ii).lt.gmin) then
                gmin=gpoly(ii)
                iimin=ii
              endif
            enddo
            if (ispoly(iimin).ne.nq(n+i)) then
              print 100,wc(nq(n+i)),wc(ispoly(iimin))
 100          format('   ',A8,' was replaced by ',A8)
              icpst=icpst+1
              ispst(icpst)=ispoly(iimin)
              icpdt=icpdt+1
              ispdt(icpdt)=nq(n+i)
              nq(n+i)=ispoly(iimin)
            endif
          endif
 20     continue
      enddo
      return
      end

      subroutine sspolymorphs(j,ipoly)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
C     search for polymorphs of solid solution
C     return solid solution # as ipoly
      ipoly=0
        if (is.eq.0) return
        if (nss.eq.0) return
      do iss=1,nss
        if (mss(iss).eq.0) goto 1000
        if (ms(j).ne.ms(iss)) goto 1000
          iflag2=0
        do jj=1,ms(j)
          f1=0d0
          do k=1,m
            if (dnu(nssp(j)+jj-1,k).eq.0d0) goto 10
            if (f1.ne.0d0) goto 10
            f1=dnu(nssp(j)+jj-1,k)
 10         continue
          enddo
          f2=0d0
          do k=1,m
            if (dnu(nssp(iss)+jj-1,k).eq.0d0) goto 11
            if (f2.ne.0d0) goto 11
            f2=dnu(nssp(iss)+jj-1,k)
 11         continue
          enddo
          iflag=0
          do k=1,m
            if (dabs(dnu(nssp(j)+jj-1,k)/f1-dnu(nssp(iss)+jj-1,k)/f2)
     >          .gt.eps) iflag=1
          enddo
          if (iflag.eq.1) iflag2=1
        enddo
        if (iflag2.eq.0) ipoly=iss
 1000 continue
      enddo
      return
      end

      subroutine calcsx()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
C     mole fractions for gaseous species
      pt=0d0
      do i=1,m+n
         pt=pt+p(i)
      enddo
      do i=1,m+n
         sxgas(i)=p(i)/pt
      enddo
C     mole fractions for solutions
      do j=1,nss
        if (mss(j).ne.0) then
         if (wss(j)(1:8).eq.'Spinel  ') then
            xcr=p(mss(j)+n)
            xal=1d0-xcr
            xfe=p(mss(j)+1+n)
            xmg=1d0-xfe
            sx(j,4)=xmg*xal
            sx(j,3)=xfe*xal
            sx(j,2)=xmg*xcr
            sx(j,1)=xfe*xcr
         else
            pt=0d0
            do i=1,ms(j)
               pt=pt+p(mss(j)+n+i-1)
            enddo
            do i=1,ms(j)
               sx(j,i)=p(mss(j)+n+i-1)/pt
            enddo
c            print *,sx(2,1)
         endif
 30     endif
      enddo
      return
      end

      subroutine calcG0()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
C     calculate the standard state free energies of formation
C     from monatomic gases at current T: ag0(i) kJ/mol
      gasconst=rgas1/1d3
      do i=1,m
            ag0(i)=0d0
      enddo
      do i=1,ns+msExtra
            gg=xc(1,i)
            do ii=2,iord
              gg=gg+xc(ii,i)/((t0/1d2)**(ii-1))
            enddo
            ag0(i+m)=-gg*dlog(10d0)*gasconst*t0
c           if(t0.eq.1000) then
c      species no., logK, dG      
c            print *,i+m
c            print *,gg
c            print *,ag0(i+m)
c           endif
      enddo
      do i=n1,ngp
            amu(i+m)=ag0(i+m)
      enddo
      if (nss.ne.0) then
        do i=1,nss
         if (wss(i)(1:8).eq.'Spinel  ') then
            do j=1,4
               ag0sp(j)=ag0(m+nssp(i)+ispinel(j)-1)
            enddo
           if (nonid(i).eq.0) then
            if (ispinel(1).eq.4) then
               ag0sp(1)=ag0sp(2)+ag0sp(3)-ag0sp(4)
            elseif (ispinel(2).eq.4) then
               ag0sp(2)=ag0sp(1)+ag0sp(4)-ag0sp(3)
            elseif (ispinel(3).eq.4) then
               ag0sp(3)=ag0sp(1)+ag0sp(4)-ag0sp(2)
            elseif (ispinel(4).eq.4) then
               ag0sp(4)=ag0sp(2)+ag0sp(3)-ag0sp(1)
            endif
           endif
         endif
        enddo
      endif
      return
      end

      subroutine freeenergy(G)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
C     calculate the free energy of the system (kJ) for 1 mole of elements
C      if (t0.eq.tstart) return
      G=0
C     Gas phase
C        total moles of gaseous species = Ptotal/Ptot/TG
C        mole of each species = partial press./Ptot/TG =p(j)/TG
      i=0
         call ssfreeenergy(i,gg)
         aNtotal=ptotal/ptot/tg
C         print *,aNtotal
C         print *,ptotal
C         print *,ptot
C         print *,tg
C         pause
         G=G+gg*aNtotal
      if (is.eq.0) goto 10
C     pure solids
      do i=n1,nis
         if (nq(i).le.ngp) then
C            print *,'pure solid',i,wc(i)
            G=G+amu(m+nq(i))*p(m+i)/tg
         endif
      enddo
      if (nss.eq.0) goto 10
C     solid or liquid solutions
      do i=1,nss
         if (mss(i).ne.0) then
C            print *,'solid solution',i,wss(i)
            call ssfreeenergy(i,gg)
            aNtotal=0d0
            if (wss(i)(1:8).eq.'Spinel  ') then
               aNtotal=p(n+mss(i)+2)/tg
            else
               do j=1,ms(i)
                  aNtotal=aNtotal+p(n+mss(i)+j-1)/tg
               enddo
            endif
            G=G+gg*aNtotal
         endif
      enddo
 10   continue
      if (idebug.ge.1) print 3030,G
      if (idebug.ge.2) write(8,3030) G
 3030 format (1X,'Total free energy of the system =',1pe12.5,' kJ')
      return
      end

      subroutine ssfreeenergy(i,G)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
C     calculate free energy of solutions (kJ/mol)
      gasconst=rgas1/1d3
      Ginfinite=1d20
      G=0d0
c      print *,ag0(30)
      if (i.eq.0) then
C        gas phase; use sxgas(m+n) = partial pressure/Ptotal = mole fraction
C                   cf. p(m+n)     = partial pressure/Ptot
C        Ptot is a given pressure at the start
C        Ptotal is a real gas pressure
C        If mode=2 or 3, Ptotal=Ptot (isobaric)
         do j=1,m+n
           if (sxgas(j).gt.1d-200) then
            amu(j)=ag0(j)+gasconst*t0*dlog(sxgas(j)*ptotal)
           endif
         enddo
         do j=1,m+n
           if (sxgas(j).gt.1d-200) then
            G=G+amu(j)*sxgas(j)
           endif
         enddo
         amugas=G
c          print *,wss(4)(1:8)
      elseif (wss(i)(1:8).eq.'Spinel  ') then
C        Spinel solid solution; use xfe and xcr
            xal=1d0-xcr
            xmg=1d0-xfe
        if ((xal.gt.0).and.(xcr.gt.0)
     >         .and.(xmg.gt.0).and.(xfe.gt.0)) then
            sx(i,4)=xmg*xal
            sx(i,3)=xfe*xal
            sx(i,2)=xmg*xcr
            sx(i,1)=xfe*xcr
            amuss(i,1)=ag0sp(1)+gasconst*t0*dlog(sx(i,1)*xcr)
            amuss(i,2)=ag0sp(2)+gasconst*t0*dlog(sx(i,2)*xcr)
            amuss(i,3)=ag0sp(3)+gasconst*t0*dlog(sx(i,3)*xal)
            amuss(i,4)=ag0sp(4)+gasconst*t0*dlog(sx(i,4)*xal)
         if (nonid(i).ne.0) then
            idf=0
            call calcgamma(i,idf)
            do j=1,4
               amuss(i,j)=amuss(i,j)+gasconst*t0*dlog(10d0)*gamma(j)
            enddo
         endif
         do j=1,4
            amu(m+nssp(i)+j-1)=amuss(i,j)
            G=G+amuss(i,j)*sx(i,j)
         enddo
        else
         G=Ginfinite
        endif
         amusstotal(i)=G
      else
C        other solid solutions or liquid; use sx(i,j)
         iflag=0
         do j=1,ms(i)
            if (sx(i,j).le.0) iflag=1
         enddo
        if (iflag.eq.0) then
         if (nonid(i).ne.0) then
            idf=0
            call calcgamma(i,idf)
         else
            do j=1,ms(i)
               gamma(j)=0d0
            enddo
         endif
         do j=1,ms(i)
            amuss(i,j)=ag0(m+nssp(i)+j-1)+dble(mdb(i))*gasconst*t0*
     >                 (dlog(sx(i,j))+dlog(10d0)*gamma(j))
            amuss(2,1)=ag0(m+nssp(2)+1-1)+dble(mdb(2))*gasconst*t0*
     >                 (dlog(sx(2,1))+dlog(10d0)*0.5*gamma(1))
            amuss(2,2)=ag0(m+nssp(2)+2-1)+dble(mdb(2))*gasconst*t0*
     >                 (dlog(sx(2,2))+dlog(10d0)*0.5*gamma(2))
            amuss(3,1)=ag0(m+nssp(3)+1-1)+dble(mdb(3))*gasconst*t0*
     >                 (dlog(sx(3,1))+dlog(10d0)*0.5*gamma(1))
            amuss(3,2)=ag0(m+nssp(3)+2-1)+dble(mdb(3))*gasconst*t0*
     >                 (dlog(sx(3,2))+dlog(10d0)*0.5*gamma(2))
c            amuss(i,j)=ag0(m+nssp(i)+j-1)+gasconst*t0*
c     >                 (2d0*dlog(sx(i,j))+dlog(10d0)*gamma(j))
c         if (wss(i)(1:8).eq.'Olivine ') then
c            amuss(2,j)=ag0(m+nssp(2)+j-1)+gasconst*t0*
c     >                 dlog(sx(i,j))+gasconst*t0*dlog(10d0)*gamma(j)
c         endif
         enddo
c         if(t0.eq.1700d0) then
c         print *,ag0(2)
c         stop
c         endif
         do j=1,ms(i)
            amu(m+nssp(i)+j-1)=amuss(i,j)
            G=G+amuss(i,j)*sx(i,j)
         enddo
        else
         G=Ginfinite
        endif
         amusstotal(i)=G
c         print *,amusstotal(4)
       endif
      if (i.ne.0) then
         if (idebug.ge.2) print 3030,wss(i),G
         if (idebug.ge.2) write(8,3030) wss(i),G
      endif
 3030 format (1X,'Free energy of ',A8,' =',1pe12.5,' kJ/mol')
      return
      end

      subroutine printout()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      dimension pw(islimit)
      CHARACTER*80 wcc(islimit)
ccc   print out info for this temperature step
C     print out a new condensate and/or a phase that disappeared
         if ((itable.ne.0).and.(ntable.lt.ntlimit)) then
            wcond(ntable+1)='        '
            wdisap(ntable+1)='        '
         endif
      if (icps.eq.1) then
         print 4105,t0,wc(isps(1))
         write(8,4005)t0,wc(isps(1))
         if (isummary.ge.1) write(7,4005)t0,wc(isps(1))
         if (isummary.ge.2) write(15,4005)t0,wc(isps(1))
         if (iwindow.ge.1) write(16,4105)t0,wc(isps(1))
         if (iwindow.ge.1) write(17,4105)t0,wc(isps(1))
         if ((itable.ne.0).and.(ist.lt.istlimit)
     >       .and.(ntable.lt.ntlimit)) then
            iflag=0
            if (ist.ne.0) then
               do i=1,ist
                  if (nqt(i).eq.isps(1)) iflag=1
               enddo
            endif
            if (iflag.eq.0) then
               ist=ist+1
               nqt(ist)=isps(1)
            endif
            wcond(ntable+1)=wc(isps(1))
         endif
      endif
      if (icss.eq.1) then
         print 4105,t0,wss(isss(1))
         write(8,4005)t0,wss(isss(1))
         if (isummary.ge.1) write(7,4005)t0,wss(isss(1))
         if (isummary.ge.2) write(15,4005)t0,wss(isss(1))
         if (iwindow.ge.1) write(16,4105)t0,wss(isss(1))
         if (iwindow.ge.1) write(17,4105)t0,wss(isss(1))
         if ((itable.ne.0).and.(ist.lt.istlimit-ms(isss(1))+1)
     >       .and.(ntable.lt.ntlimit)) then
            iflag=0
            if (ist.ne.0) then
               do i=1,ist
                  if (nqt(i).eq.nssp(isss(1))) iflag=1
               enddo
            endif
            if (iflag.eq.0) then
               do i=1,ms(isss(1))
                  ist=ist+1
                  nqt(ist)=nssp(isss(1))+i-1
               enddo
            endif
            wcond(ntable+1)=wss(isss(1))
         endif
      endif
      if (icpd.eq.1) then
         if (ispd(1).gt.ngp) then
            do iss=1,nss
               if (ispd(1).eq.nssp(iss)+ms(iss)-1) then
                  print 4108,t0,wss(iss)
                  write(8,4008)t0,wss(iss)
                  if (isummary.ge.1) write(7,4008)t0,wss(iss)
                  if (isummary.ge.2) write(15,4008)t0,wss(iss)
                  if (iwindow.ge.1) write(16,4108)t0,wss(iss)
                  if (iwindow.ge.1) write(17,4108)t0,wss(iss)
                  if ((itable.ne.0).and.(ntable.lt.ntlimit)) then
                     wdisap(ntable+1)=wss(iss)
                  endif
               endif
            enddo
         else
            print 4108,t0,wc(ispd(1))
            write(8,4008)t0,wc(ispd(1))
            if (isummary.ge.1) write(7,4008)t0,wc(ispd(1))
            if (isummary.ge.2) write(15,4008)t0,wc(ispd(1))
            if (iwindow.ge.1) write(16,4108)t0,wc(ispd(1))
            if (iwindow.ge.1) write(17,4108)t0,wc(ispd(1))
                  if ((itable.ne.0).and.(ntable.lt.ntlimit)) then
                     wdisap(ntable+1)=wc(ispd(1))
                  endif
         endif
      endif
      if (icps+icss+icpd.gt.0) then
         write(8,4002)
      endif
 4005 FORMAT('New Condensate (',F7.2,'K)  : ',A50)
 4008 FORMAT('PHASE DISAPPEARS (',F7.2,'K): ',A50)
 4105 FORMAT(' New Condensate (',F7.2,'K)  : ',A50)
 4108 FORMAT(' PHASE DISAPPEARS (',F7.2,'K): ',A50)
C     WRITE(10,3) (P(M+I),I=1,N)
C3    FORMAT(1X,1P12E10.3)
C     total pressure
 61   PT=0d0
      DO 25 I=1,m+n
 25   PT=PT+P(I)
      PT=PT*PTOT
      ptotal=pt
C     log10(P(O2))
      PO2=log10(P(m+nO2)*PTOT)
      if (ibar.eq.0) then
      print 299,T0,pt,NITER,po2,IS
      else
      print 1299,T0,pt,NITER,po2,IS
      endif
 299  FORMAT(1X,'T= ',F7.2,'K  P=',1PE12.5,' atm',
     1'  Iterations=',I3,'  LogPO2=',0pf7.2,'  # Solids=',I3)
1299  FORMAT(1X,'T= ',F7.2,'K  P=',1PE12.5,' bar',
     1'  Iterations=',I3,'  LogPO2=',0pf7.2,'  # Solids=',I3)
C      print 299,T0,NITER,IS
C 299  FORMAT(1X,'T = ',F7.2,'K  # ITERATIONS =',I3,
C     >' # SOLIDS =',I3)
      if (ibar.eq.0) then
      WRITE(8,2)T0,PT,rgas2*t0,po2,is
      WRITE(8,30)Vtotal,pt/ptot/TG
      if (isummary.ge.1) WRITE(7,2)T0,PT,rgas2*t0,po2,sx(2,1),is
      if (isummary.ge.2) WRITE(15,2)T0,PT,rgas2*t0,po2,sx(2,1),is
      if (iwindow.ge.1) WRITE(16,12)T0,PT,rgas2*t0,po2,sx(2,1),is
      if (iwindow.ge.1) WRITE(17,12)T0,PT,rgas2*t0,po2,sx(2,1),is
      else
      WRITE(8,1002)T0,PT,rgas2*t0,po2,sx(2,1),is
      WRITE(8,30)Vtotal,pt/ptot/TG
      if (isummary.ge.1) WRITE(7,1002)T0,PT,rgas2*t0,po2,sx(2,1),is
      if (isummary.ge.2) WRITE(15,1002)T0,PT,rgas2*t0,po2,sx(2,1),is
      if (iwindow.ge.1) WRITE(16,1012)T0,PT,rgas2*t0,po2,sx(2,1),is
      if (iwindow.ge.1) WRITE(17,1012)T0,PT,rgas2*t0,po2,sx(2,1),is
      endif
 2    FORMAT('T= ',F7.2,'K  P=',1PE12.5,' atm',
     >'  RT=',1pe12.5,'  LogPO2=',0pf7.2,'  Xfa=',f7.5,'  # Solids=',I3)
 12   FORMAT(1X,'T= ',F7.2,'K  P=',1PE12.5,' atm',
     >'  RT=',1pe12.5,'  LogPO2=',0pf7.2,'  Xfa=',f7.5,'  # Solids=',I3)
1002  FORMAT('T= ',F7.2,'K  P=',1PE12.5,' bar',
     >'  RT=',1pe12.5,'  LogPO2=',0pf7.2,'  Xfa=',f7.5,'  # Solids=',I3)
1012  FORMAT(1X,'T= ',F7.2,'K  P=',1PE12.5,' bar',
     >'  RT=',1pe12.5,'  LogPO2=',0pf7.2,'  Xfa=',f7.5,'  # Solids=',I3)
 30   FORMAT('V=',1PE12.5,' l (for 1 mol of elements)  N gas=',
     >1PE12.5,' mol')
c**   print out pure solid phases
      if (is.eq.0) goto 4004
      j=0
      do 4003 i=n1,nis
      if(nq(i).le.ngp) then
         j=j+1
         wcc(j)=wc(nq(i))
         pw(j)=P(I+M)*ptot/(rgas2*t0)
      endif
 4003 continue
      if (j.eq.0) goto 4004
      write(8,3)(wcc(i),pw(i),i=1,j)
c      print *,(pw(i),i=1,j)
      if (iwindow.ge.1) write(17,13)(wcc(i),pw(i),i=1,j)
 3    format('Condensed phases (mol/l): ',2(a8,1pE12.5,'      ',:),
     >8(/3(a8,E12.5,'      ',:)))
 13   format(' Condensed phases (mol/l): ',2(a8,1pE12.5,'      ',:),
     >8(/1X,3(a8,E12.5,'      ',:)))
 4004 continue
c**   print out solid solutions
      DO 19 I=1,NSS
      K=MS(I)+mse(i)
c      kik=ms(i)+nssp(i)-1
      if(i.eq.nssspinel) then
C      print out spinel solid solution
       if(mss(i).ne.0) then
         pal=p(mss(i)+ms(i)-1+n)*ptot/(t0*rgas2)
         xcr=p(mss(i)+n)
         xal=1d0-xcr
         xfe=p(mss(i)+1+n)
         xmg=1d0-xfe
         xksp(2)=dlp(mmg)+2*dlp(mcr)+4*dlp(mO)-dlog10(xmg)-2*dlog10(xcr)
         xksp(1)=dlp(mfe)+2*dlp(mcr)+4*dlp(mO)-dlog10(xfe)-2*dlog10(xcr)
         xksp(3)=dlp(mfe)+2*dlp(mal)+4*dlp(mO)-dlog10(xfe)-2*dlog10(xal)
         xksp(4)=dlp(mmg)+2*dlp(mal)+4*dlp(mO)-dlog10(xmg)-2*dlog10(xal)
         if (nonid(i).ne.0) then
            call calcG0()
            idf=0
            call calcgamma(i,idf)
            do j=1,4
               xksp(j)=xksp(j)-gamma(j)
            enddo
         endif
         write(8,225) wss(i),pal,xal,xmg,xfe,xcr
         if (iwindow.ge.1) write(17,226) wss(i),pal,xal,xmg,xfe,xcr
 225     format(A8,' :',1pe12.5,' mol/l',/,
     >   ' Xal=',1pe12.5,'  Xmg=',e12.5,'  Xfe=',e12.5,'  Xcr=',e12.5)
 226     format(1X,A8,' :',1pe12.5,' mol/l',/,
     >   '  Xal=',1pe12.5,'  Xmg=',e12.5,'  Xfe=',e12.5,'  Xcr=',e12.5)
         do 3006 isp=4,1,-1
         if (ispinel(isp).ne.4) goto 3006
         t1t=1d2/t0
         xkt=xc1(ispinel(isp)+nssp(i)-1)
         do kt=2,iord
            xkt=xkt-xc(kt,ispinel(isp)+nssp(i)-1)*(t1t**(kt-1))
         enddo
         difference=xksp(isp)-xkt
         write(8,3007)wsp(isp),xksp(isp),xkt,difference
 3007    format(' checks: logK(',A7,') out/in:',1p2e13.5, ' diff.:',
     >   e13.5)
 3006    continue
         write(8,321) (wsp(j),sx(i,j),j=1,4)
         if (iwindow.ge.1)
     >    write(17,323) (wsp(j),sx(i,j),j=1,4)
 321     format(3(3(' X(',A8,')=',1pe12.5,:' '),/))
 323     format(3(1X,3(' X(',A8,')=',1pe12.5,:' '),/))
C          if non-ideal solid solution, print activity coeff.
           if(nonid(i).ne.0) then
            write(8,322) (wsp(j),10**gamma(j),j=1,4)
            if (iwindow.ge.1)
     >      write(17,324) (wsp(j),10**gamma(j),j=1,4)
 322        format(3(3(' g(',A8,')=',1pe12.5,:' '),/))
 324        format(3(1X,3(' g(',A8,')=',1pe12.5,:' '),/))
           endif
       endif
      else
C      print out other solid solutions
       if(mss(i).gt.0) then
c         print *,p(2)
         pal=p(mss(i)+ms(i)-1+n)*ptot/(sx(i,ms(i))*t0*rgas2)
         write(8,221) wss(i),pal,(wc(nssp(i)+j-1),sx(i,j),j=1,k)
         if (iwindow.ge.1)
     >    write(17,223) wss(i),pal,(wc(nssp(i)+j-1),sx(i,j),j=1,k)
c         print*,wc(nssp(3)+2-1)
 221     format(A8,' :',1pe12.5,' mol/l',/,
     >      3(3(' X(',A8,')=',1pe12.5,:' '),/))
 223     format(1X,A8,' :',1pe12.5,' mol/l',/,
     >      3(1X,3(' X(',A8,')=',1pe12.5,:' '),/))
C        if non-ideal solid solution, print activity coeff.
         if(nonid(i).ne.0) then
            idf=0
            call calcgamma(i,idf)
            write(8,222) (wc(nssp(i)+j-1),10**gamma(j),j=1,k)
            if (iwindow.ge.1)
     >      write(17,224) (wc(nssp(i)+j-1),10**gamma(j),j=1,k)
 222        format(3(3(' g(',A8,')=',1pe12.5,:' '),/))
 224        format(3(1X,3(' g(',A8,')=',1pe12.5,:' '),/))
         endif
       endif
      endif
 19   continue
      write(8,4002)
      if (iwindow.ge.1) write(17,4002)
 4002 format()
c**   print out each element
      DO 20 I=1,M
      RP=R(I)*TG
      NN=1
      pw(1)=P(I)/RP
      difference=1d0-P(I)/RP
      wcc(1)=wee(i)
      IF(pw(1).LT.1.D-4) NN=0
      DO 21 J=1,NIS
      IF(NU(NQ(J),I)) 22,21,22
 22   RNT=P(J+M)*dnu(NQ(J),I)/RP
      if(nssspinel.ne.0) then
      if(mss(nssspinel).ne.0) then
        if((j+m-n.ge.mss(nssspinel)).and.
     >   (j+m-n.le.mss(nssspinel)+ms(nssspinel)-1)) go to 21
      endif
      endif
      difference=difference-rnt
      IF(RNT-1.d-4) 21,23,23
C     print only species which consume >=.0001 of the abundance of the element
 23   NN=NN+1
      pw(NN)=RNT
      wcc(nn)=wc(nq(j))
      if (itable.ne.0) then
         itbflag(i,nq(j))=1
      endif
 21   CONTINUE
ccc   Add spinel here, if any.
      if(nssspinel.eq.0) goto 2002
      if(mss(nssspinel).eq.0) goto 2002
         rnt=0d0
         pt=p(mss(nssspinel)+ms(nssspinel)-1+n)
         if(i.eq.mal) rnt=2d0*xal*pt/rp
         if(i.eq.mcr) rnt=2d0*xcr*pt/rp
         if(i.eq.mmg) rnt=xmg*pt/rp
         if(i.eq.mfe) rnt=xfe*pt/rp
         if(i.eq.mO)  rnt=4d0*pt/rp
         difference=difference-rnt
         if(rnt .ge. 1.d-4) then
            nn=nn+1
            pw(nn)=rnt
            wcc(nn)='Spinel'
            if (itable.ne.0) then
               do j=1,ms(nssspinel)
                  itbflag(i,nssp(nssspinel)+j-1)=1
               enddo
            endif
         endif
 2002 continue
      if (ibar.eq.0) then
      WRITE(8,4)WE(I),r(i)*tg*ptot/(rgas2*t0),we(i),p(i)*ptot,
     >difference*r(i)*tg*ptot/(rgas2*t0),(pw(J),wcc(j),J=1,NN)
      else
      WRITE(8,1004)WE(I),r(i)*tg*ptot/(rgas2*t0),we(i),p(i)*ptot,
     >difference*r(i)*tg*ptot/(rgas2*t0),(pw(J),wcc(j),J=1,NN)
      endif
 20   continue
 4    FORMAT(A2,' N total=',1PE12.5,' mol/l ',
     >' P(',A2,'(g))=',e12.5,' atm ',
     >'  N-Sum=',e8.1,' mol/l',/,' Fractions :        ',
     >3(e10.3,1x,a8,:,1X),20(/4(e10.3,1x,a8,:,1X)))
1004  FORMAT(A2,' N total=',1PE12.5,' mol/l ',
     >' P(',A2,'(g))=',e12.5,' bar ',
     >'  N-Sum=',e8.1,' mol/l',/,' Fractions :        ',
     >3(e10.3,1x,a8,:,1X),20(/4(e10.3,1x,a8,:,1X)))
      write(8,4002)
      if (itable.ne.0) then
         ntable=ntable+1
         if (ntable.le.ntlimit) then
            ttable(ntable)=t0
            tgtb(ntable)=tg
            pttb(ntable)=ptotal
            pO2tb(ntable)=pO2
            istb(ntable)=is
            ptottb(ntable)=ptot
            do i=1,m
               pmtable(i,ntable)=p(i)
            enddo
            if (is.ne.0) then
               do i=1,is
                  do j=1,ist
                     if (nqt(j).eq.nq(n+i)) then
                         pstable(j,ntable)=p(i+m+n)
                     endif
                  enddo
               enddo
            endif
         endif
         if (itable.ge.4) call wttables()
      endif
      return
      END

      subroutine readcoeff(filename)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      CHARACTER*12 filename
      CHARACTER*80 a
      INCLUDE 'con0621i.f'
C     read thermodynamic data file which consists of coefficients of logK
      OPEN(1, FILE = filename)
      read(1,'(A80)') a
      comdata=a
      print '(1X,A79)',a
      write(8,'(A79)') a
      if (isummary.ge.1) write(7,'(A79)') a
      read(1,*) m
      if (m.gt.mlimit) then
         print *,' Exceed M limit !!!'
         stop
      endif
      read(1,'(40A2)') (we(k),k=1,m)
      do k=1,m
         wee(k)=we(k)
         wee(k)(3:5)='(g)'
      enddo
      read(1,*) n,npure,nss
      if (n.gt.nlimit) then
         print *,' Exceed N limit !!!'
         stop
      endif
      if (n+npure.gt.nslimit) then
         print *,' Exceed NS limit !!!'
         stop
      endif
      if (nss.gt.nsslimit) then
         print *,' Exceed NSS limit !!!'
         stop
      endif
      ngp=npure+n
      nsst=0
      msExtra=0
      if (nss.eq.0) goto 100
      read(1,'(10A8)') (wss(k),k=1,nss)
      read(1,*) (ms(k),mdb(k),mse(k),nonid(k),k=1,nss)
      do i=1,nss
      mss(i)=0
      nssp(i)=nsst+msExtra+ngp+1
      nsst=nsst+ms(i)
      msExtra=msExtra+mse(i)
      if (ms(i)+mse(i).gt.isslimit) then
         print *,' Exceed ISS limit !!!'
         stop
      endif
      enddo
 100  ns=ngp+nsst
      if (ns+msExtra.gt.nslimit) then
         print *,' Exceed NS limit !!!'
         stop
      endif
      read(1,*)iord
      if (iord.gt.iordlimit) then
         print *,' Exceed iord limit !!!'
         stop
      endif
      DO 11 J=1,NS+msExtra
      READ(1,103) (XC(I,J),I=1,IORD),(NU(J,I),I=1,M),wc(j)
 103  FORMAT(1P3E12.5,15I2,1x,a8)
      nspecies(j)=j
      wref(j)='  '
      DO 45 K=1,M
 45   dnu(j,k)=dble(NU(J,K))
 11   continue
      close(1)
      return
      end

      subroutine readdg(filename)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      INCLUDE 'con0621i.f'
      parameter(ntemplimit=15)
      DOUBLE PRECISION F(mlimit,ntemplimit),G(nslimit,ntemplimit),
     >   DG(nslimit,ntemplimit),ECF(nslimit,ntemplimit)
C                      F(m,ntemp),G(ns,ntemp),DG(ns,ntemp),ECF(ns,ntemp)
      DOUBLE PRECISION DACT(ntemplimit),DPRIME(iordlimit),
     >                 DEV(ntemplimit)
C                      DACT(ntemp),DPRIME(iord),DEV(ntemp)
      DOUBLE PRECISION AA(ntemplimit,iordlimit),
     >            AT(iordlimit,ntemplimit),ATA(iordlimit,iordlimit+1)
C                      AA(ntemp,iord),AT(iord,ntemp),ATA(iord,iord+1)
      INTEGER IT(ntemplimit)
C             IT(ntemp)
      character*15 a10(ntemplimit)
C                  a10(ntemp)
      CHARACTER*12 filename
      CHARACTER*80 A,wee2(mlimit)
      m2=m2wt
C     read thermodynamic data file which consists of delta Gs of formation,
C     and calculate the coefficients of logKf from monatomic gases.
      ifit=iord
      if (iord.gt.iordlimit) then
         print *,' Exceed iord limit !!!'
         stop
      endif
      OPEN(1, FILE = filename)
      i=index(filename,'.')
      j=index(filename,' ')
      if (idebug.ge.1) then
         if (i.eq.0) then
            if ((j.eq.0).or.(j.gt.9)) then
               filename(9:12)='.CEF'
            else
               filename(j:j+3)='.CEF'
            endif
         else
            filename(i:i+3)='.CEF'
         endif
         open (2,file=filename)
      endif
      nelement=0
      ngas=0
      npure=0
      nsolution=0
      icount=0
   10 format(40I2)
   11 format(D10.3)
   12 FORMAT(5(A10,1X))
   13 format(I6)
      call geta(a)
      comdata=a
      print '(1X,A79)',a
      write(8,'(A79)') a
      if (isummary.ge.1) write(7,'(A79)') a
      call geta(a)
        print '(A80)', a
      read(a,'(I2,1X,I2)') m,ntemp
      ntstep=(ntemp+4)/5
      print *, m,ntemp
      if (m.gt.mlimit) then
         print *,' Exceed M limit !!!'
         stop
      endif
      if (ntemp.gt.ntemplimit) then
         print *,' Exceed ntemp limit !!!'
         stop
      endif
      call geta(a)
      read(a,'(40A2)') (we(k),k=1,m)
      do k=1,m
         if(we(k).eq.'Fe') mfe=k
C        print *,we(k)
      enddo
        print *, m,mfe
C          pause
C     read data for monatomic gases
  100 call geta(a)
      print *,A(1:5)
      IF (A(1:5) .ne. 'Start') goto 100
  101 call geta(a)
      if (A(1:3) .eq. 'End') goto 150
      nelement=nelement+1
      read(a,13) nspecies(nelement)
      print *,nspecies(nelement)
C      print *,'Hi'
C      pause

      if (nelement.eq.1) goto 102
      do k=1,nelement-1
         if (nspecies(nelement).eq.nspecies(k)) then
            print *,' Duplicate species number ! : ',nspecies(nelement)
            stop
         endif
      enddo
  102 call geta(a)
      wee2(nelement)=a(3:)
  103 call geta(a)
      wref(nelement)=a(3:)
  104 call geta(a)
      read(a,10)(nu(nelement,k),k=1,m)
      do nt=1,ntstep
        call geta(a)
        if (nt.eq.ntstep) then
           read(A,12) (a10(k),k=1+(nt-1)*5,ntemp)
        else
           read(A,12) (a10(k),k=1+(nt-1)*5,nt*5)
        endif
      enddo
      DO J=1,ntemp
        IF(a10(j) .EQ. '          ') then
           G(nelement,J)=999999.999D0
        else
           READ(A10(J),11) G(nelement,J)
        endif
      enddo
      goto 101
  150 continue
c      print *,G(1,15) 
      do i=1,nelement
      do k=1,nelement
         if (nu(i,k).ne.0) then
            do jj=1,ntemp
               f(k,jj)=g(i,jj)
c      print *,g(1,1)
            enddo
            wee(k)=wee2(i)
         endif
      enddo
      enddo
C     read data for gaseous compounds
  200 call geta(a)
      IF (A(1:5) .ne. 'Start') goto 200
  201 call geta(a)
      if (A(1:3) .eq. 'End') goto 250
      icount=icount+1
      ngas=ngas+1
      if (ngas.gt.nlimit) then
         print *,' Exceed N limit !!!'
         stop
      endif
      read(a,13) nspecies(icount)
c      print *,nspecies(icount)
c      print *,icount
      if (icount.eq.1) goto 202
      do k=1,icount-1
         if (nspecies(icount).eq.nspecies(k)) then
            print *,' Duplicate species number ! : ',nspecies(icount)
            stop
         endif
      enddo
  202 call geta(a)
      wc(icount)=a(3:)
  203 call geta(a)
      wref(icount)=a(3:)
  204 call geta(a)
      read(a,10)(nu(icount,k),k=1,m)
      do nt=1,ntstep
        call geta(a)
        if (nt.eq.ntstep) then
           read(A,12) (a10(k),k=1+(nt-1)*5,ntemp)
        else
           read(A,12) (a10(k),k=1+(nt-1)*5,nt*5)
        endif
      enddo
      DO J=1,ntemp
        IF(a10(j) .EQ. '          ') then
           G(icount,J)=999999.999D0
        else
           READ(A10(J),11) G(icount,J)
C           print *, G(icount,J)
C           print *, icount
        endif
      enddo

      goto 201
  250 continue
C      print *,G(29,15) 
C      pause
C     read data for pure solids or liquids
  300 call geta(a)
      IF (A(1:5) .ne. 'Start') goto 300
  301 call geta(a)
      if (A(1:3) .eq. 'End') goto 350
      icount=icount+1
      npure=npure+1
      if (icount.gt.nslimit) then
         print *,' Exceed NS limit !!!'
         stop
      endif
      read(a,13) nspecies(icount)
c      print *,nspecies(icount)
      if (icount.eq.1) goto 302
      do k=1,icount-1
         if (nspecies(icount).eq.nspecies(k)) then
            print *,' Duplicate species number ! : ',nspecies(icount)
            stop
         endif
      enddo
  302 call geta(a)
      wc(icount)=a(3:)
c        print *,wc(icount)(1:7)
c      if(wc(icount)(1:7).eq.'FeS sol') print *,icount
      if(wc(icount)(1:7).eq.'Fe.947O') nwustite=icount
      if(wc(icount)(1:7).eq.'Fe.98S') nfe98s=icount
      if(wc(icount)(1:7).eq.'Fe.875S') nfe875s=icount
  303 call geta(a)
      wref(icount)=a(3:)
  304 call geta(a)
      read(a,10)(nu(icount,k),k=1,m)
      call geta(a)
      read(a,10) kind(icount),(nu2(icount,k),k=1,m2)
      do nt=1,ntstep
        call geta(a)
        if (nt.eq.ntstep) then
           read(A,12) (a10(k),k=1+(nt-1)*5,ntemp)
        else
           read(A,12) (a10(k),k=1+(nt-1)*5,nt*5)
        endif
      enddo
      DO J=1,ntemp
        IF(a10(j) .EQ. '          ') then
           G(icount,J)=999999.999D0
        else
           READ(A10(J),11) G(icount,J)
        endif
      enddo
      goto 301
  350 continue
c      print *,G(135,15) 
C     read data for solid solutions
  400 call geta(a)
      IF (A(1:5) .ne. 'Start') goto 400
  401 call geta(a)
      if (A(1:3) .eq. 'End') goto 450
      nsolution=nsolution+1
      if (nsolution.gt.nsslimit) then
         print *,' Exceed NSS limit !!!'
         stop
      endif
      wss(nsolution)=a(3:)
  410 call geta(a)
      read(a,'(4I2)') ms(nsolution),mdb(nsolution),mse(nsolution),
     >                nonid(nsolution)
      if (ms(nsolution)+mse(nsolution).gt.isslimit) then
         print *,' Exceed ISS limit !!!'
         stop
      endif
c     read data for each endmenber
      do j=1,ms(nsolution)+mse(nsolution)
  411 call geta(a)
      icount=icount+1
      if (icount.gt.nslimit) then
         print *,' Exceed NS limit !!!'
         stop
      endif
      read(a,13) nspecies(icount)
c      print *,nspecies(icount)
      if (icount.eq.1) goto 402
      do k=1,icount-1
         if (nspecies(icount).eq.nspecies(k)) then
            print *,' Duplicate species number ! : ',nspecies(icount)
            stop
         endif
      enddo
      
  402 call geta(a)
      wc(icount)=a(3:)
  403 call geta(a)
      wref(icount)=a(3:)
  404 call geta(a)
      read(a,10)(nu(icount,k),k=1,m)
      call geta(a)
      read(a,10) kind(icount),(nu2(icount,k),k=1,m2)
      do nt=1,ntstep
        call geta(a)
        if (nt.eq.ntstep) then
           read(A,12) (a10(k),k=1+(nt-1)*5,ntemp)
        else
           read(A,12) (a10(k),k=1+(nt-1)*5,nt*5)
        endif
      enddo
      
      DO k=1,ntemp
        IF(a10(k) .EQ. '          ') then
           G(icount,k)=999999.999D0
        else
           READ(A10(k),11) G(icount,k)
        endif
      enddo
      enddo
      goto 401
  450 continue
      close(1)
      m=nelement
      n=ngas
      ngp=ngas+npure
      nss=nsolution
      nsst=0
      msExtra=0
      if (nss.eq.0) goto 500
      do i=1,nss
      mss(i)=0
      nssp(i)=nsst+msExtra+ngp+1
      nsst=nsst+ms(i)
      msExtra=msExtra+mse(i)
      enddo
      
 500  ns=ngp+nsst
c      print *,m,n,ns,nsst,icount
      
      do i=1,icount
        do k=1,m
          dnu(i,k)=dble(nu(i,k))
        enddo
      enddo
      
c       dnu(122,4)=0.947 : Fe in Wustite
c       print *,nfe875s,mfe
        dnu(nwustite,mfe)=0.947d0
        
c        dnu(nfe98s,mfe)=0.98d0
        
        dnu(nfe875s,mfe)=0.875d0
c       print *,dnu(138,4)
        
C     calculate free energies of reaction from monatomic gases
  700 DO 7250 J=1,ntemp
        do i=1,icount
          IF (G(i,J) .EQ. 999999.999D0) then
             DG(i,J)=999999.999D0
          else
             DG(i,J)=G(i,J)
             do k=1,m
                dg(i,j)=dg(i,j)-dnu(i,k)*f(k,j)
             enddo
c          print *,dg(903,10)
          endif
        enddo
7195  ITEMP=2000-(J*100)
      if (idebug.ge.1) then
      WRITE(2,'(17X,A)')'Free Energies of Reaction from Monatomic Gases'
      WRITE(2,'(32X,''(Kilojoules/mole)'')')
7200  WRITE(2,'(/,33X,''Temperature = '',I4,'' K'')')ITEMP
7205  WRITE(2,'(/)')
7210  WRITE(2,'(4(1X,I6,1X,F10.3,1X))')(nspecies(i),DG(I,J),I=1,icount)
      endif
C     calculate logarithms of equilibrium constants of formation reactions
7215  DO 7227 I=1,icount
c      print *,DG(I,10)
7217  IF (DG(I,J) .EQ. 999999.999D0) GOTO 7224
c7219  ECD(I,J)=DG(I,J)/(dLOG(10d0)*rgas1/1d3*dble(ITEMP))
7220  ECF(I,J)=(-DG(I,J))/(dLOG(10d0)*rgas1/1d3*dble(ITEMP))
7222  GOTO 7227
7224  continue
c7224  ECD(I,J)=999999.999D0
7226  ECF(I,J)=999999.999D0
7227  CONTINUE
c      print *,ECF(203,10)
      if (idebug.ge.1) then
7230  WRITE(2,'(/,/,11X,''Logarithms of Equilibrium Constants of Formati
     *on Reactions'',/)')
7232  WRITE(2,'(4(1X,I6,1X,F10.4,1X))')(nspecies(i),ECF(I,J),I=1,icount)
7234  WRITE(2,'(/,/)')
C7236  WRITE(2,'(9X,''Logarithms of Equilibrium Constants of Decompositio
C     *n Reactions'',/)')
C7238  WRITE(2,'(4(1X,I6,1X,F10.4,1X))')(nspecies(i),ECD(I,J),I=1,icount)
C7240  WRITE(2,'(/,/)')
      endif
7250  CONTINUE
C     calculate coefficients of logKf by least squares fitting
      if (idebug.ge.1) then
      WRITE(2,'(15X,''Coefficients and Deviations from Least Squares Fit
     *'')')
      endif
c      print *,icount
      DO 7700 I=1,icount
      Ndata=0
      DO 50 J=1,ntemp
      IF (ECF(I,J) .EQ. 999999.999D0) GO TO 50
      Ndata=Ndata+1
40    IT(Ndata)=2000-(100*J)
      DACT(Ndata)=ECF(I,J)
50    CONTINUE
      ifittemp=0
      if (ndata.lt.ifit) then
         ifittemp=ifit
         ifit=ndata
      endif
7500  CALL FIT(ntemplimit,iordlimit,
     >         DPRIME,Ndata,DACT,IT,ifit,aa,at,ata)
      if (ifittemp.ne.0) then
         ifit=ifittemp
         do k=ndata+1,ifit
            DPRIME(K)=0d0
         enddo
      endif
      if (idebug.ge.1) then
      write(2,'(/,1X,A75)') wc(i)
      WRITE(2,'(1X,I6,5(3X,D11.5))') nspecies(i),(DPRIME(K),K=1,IFIT)
      endif
      SUM=0.0
c      print *,Ndata
      DO 7600 L=1,Ndata
      DEV(L)=DACT(L)
      DO 7550 k=1,IFIT
      DEV(L)=DEV(L)-(DPRIME(k)/(DFLOAT(IT(L))**(k-1)))
7550  CONTINUE
      SUM=SUM+DABS(DEV(L))
7600  CONTINUE
      SUM=SUM/Ndata
      if (idebug.ge.1) then
      WRITE(2,'(4(1X,I4,2X,F10.5,1X))')(IT(L),DEV(L),L=1,Ndata)
      WRITE(2,'(1X,''Mean deviation = '',F10.5)') SUM
      endif
      if (sum.gt.0.1) then
         print *,' Bad fit for Log K !!!'
         stop
      endif
      do k=1,ifit
        xc(k,i)=dprime(k)/(1d2**(k-1))
      enddo
7700  CONTINUE
      return
      END
      
      SUBROUTINE FIT(Ndim,Idim,DPRIME,N,DACT,IT,ifit,aa,at,ata)
      DOUBLE PRECISION AA(Ndim,Idim),AT(Idim,Ndim)
      DOUBLE PRECISION DACT(Ndim),DPRIME(Idim),ATA(Idim,Idim+1)
      INTEGER IT(N)
      DO 200 K=1,N
      DO 200 M=1,IFIT
      AA(K,M)=(1./DFLOAT(IT(K)))**(M-1)
      AT(M,K)=(1./DFLOAT(IT(K)))**(M-1)
200   CONTINUE
      DO 300 K=1,IFIT
      DO 300 L=1,IFIT
      ATA(K,L)=0.0
      DO 300 M=1,N
      ATA(K,L)=ATA(K,L)+(AT(K,M)*Aa(M,L))
300   CONTINUE
      DO 400 K=1,IFIT
      DPRIME(K)=0.0
      DO 400 M=1,N
      DPRIME(K)=DPRIME(K)+(AT(K,M)*DACT(M))
400   CONTINUE
      do k=1,ifit
         ata(k,ifit+1)=dprime(k)
      enddo
      CALL MATINV(Idim,Idim+1,ATA,IFIT,IFIT+1)
      do k=1,ifit
         dprime(k)=ata(k,ifit+1)
      enddo
      RETURN
      end

      subroutine geta(a)
      character*80 a
      common/flag/idebug,isummary
  100 READ(1,'(A80)') A
      print *,a(1:33)
      if (idebug.ge.3) print '(1X,A79)',a(1:79)
      IF (A(1:1) .EQ. 'C') GOTO 100
      return
      end

      SUBROUTINE MATINV(Mdim,Ndim,A,M,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(Mdim,Ndim)
      common/flag/idebug,isummary
C     scaling
      do i=1,M
         r=0d0
         do j=1,M
            c=dabs(A(i,j))
            if (c.gt.r) r=c
         enddo
         do j=1,N
            A(i,j)=A(i,j)/r
         enddo
      enddo
C     Gaussian elimination with partial pivoting
      do k=1,M-1
         l=k
         do i=k+1,M
            if (dabs(A(i,k)).gt.dabs(A(l,k))) l=i
         enddo
         if (l.ne.k) then
            do j=k,N
               temp=A(k,j)
               A(k,j)=A(l,j)
               A(l,j)=temp
            enddo
         endif
         if (dabs(A(k,k)).lt.1d-150) then
            if (idebug.ge.1) then
               iflag=0
               do i=1,M
                  print 3031,i,(A(i,j),j=1,N)
                  if (idebug.ge.2) then
                     write(8,3031) i,(A(i,j),j=1,N)
                  endif
               enddo
 3031          format (1X,I2,9(1pe8.1),10(/,3X,9(e8.1,:)))
               write(*,'(1X,A,I3,A)')
     >           'Pivot (',k,') is too small. Stop here ? (Yes=1/No=0)'
               read(*,*) iflag
               if (iflag.ne.0) then
                  stop
               endif
            endif
            print *,' Warning --- Pivot for Matrix inv. is too small !'
            if (idebug.ge.2)  write(8,*)
     >           ' Warning --- Pivot for Matrix inv. is too small !'
            A(k,k)=dsign(1d-150,A(k,k))
         endif
         do i=k+1,M
            f=A(i,k)/A(k,k)
            do j=k+1,N
               A(i,j)=A(i,j)-f*A(k,j)
            enddo
         enddo
      enddo
         if (dabs(A(M,M)).lt.1d-150) then
            if (idebug.ge.1) then
               iflag=0
               do i=1,M
                  print 3031,i,(A(i,j),j=1,N)
                  if (idebug.ge.2) then
                     write(8,3031) i,(A(i,j),j=1,N)
                  endif
               enddo
c 3031          format (1X,I2,9(1pe8.1),10(/,3X,9(e8.1,:)))
               write(*,'(1X,A,I3,A)')
     >           'Pivot (',M,') is too small. Stop here ? (Yes=1/No=0)'
               read(*,*) iflag
               if (iflag.ne.0) then
                  stop
               endif
            endif
            print *,' Warning --- Pivot for Matrix inv. is too small !'
            if (idebug.ge.2)  write(8,*)
     >           ' Warning --- Pivot for Matrix inv. is too small !'
            A(M,M)=dsign(1d-150,A(M,M))
         endif
C     backsubstitution
      do k=M+1,N
         A(M,k)=A(M,k)/A(M,M)
         do i=M-1,1,-1
            S=0d0
            do j=i+1,M
               S=S+A(i,j)*A(j,k)
            enddo
            A(i,k)=(A(i,k)-S)/A(i,i)
         enddo
      enddo
      return
      end
