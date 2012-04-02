C     Condensation Calculation Program "cwin1115" for Windows
C     "con1115A.for" source file; part one of six (A,B,C,D,E and F)
C     Designed by          L. Grossman  GCA 36,597-619         (1972)
C     Original program by  J. Lattimer  Ap. J. 219,230-249     (1978)
C     Modified by          S. Yoneda    LPS XXV 1533-1534      (1994)
C                                       Meteoritics 29,554-555 (1994)
C
$DEFINE mswindows
C $UNDEFINE mswindows
$IF DEFINED (mswindows)
      INCLUDE 'FLIB.FI'
$ENDIF
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
C
      CHARACTER*80 A80
      CHARACTER*12 filename,files(nfilelimit)
C
      data version/'cwin1115.exe 11/16/94'/
C
      data disappear/1d-11/
      data sscheck,ss2check,factor/1.0d0,0.99d0,1d-4/
      data eps/1.0D-12/
      data m2wt,nfeas/9,4/
      data rgas1/8.314510d0/
C      gas constant (J/mol-K) the 1986 recommended value
C      Cohen & Taylor J. Res. Natl Bureau Standards,92,85 (1987)
C      rgas2 (l-atm/mol-K) or (l-bar/mol-K); rgas3 (cal/mol-K)
C      ibar=1
C      if (ibar.eq.1) then
C         rgas2=rgas1/1d2
C      else
C         rgas2=rgas1/1.01325d2
C      endif
      rgas3=rgas1/4.184d0
$IF DEFINED (mswindows)
      iwindow=1
$ELSE
      iwindow=0
$ENDIF
C
C     Main
C
      print *, 'Program : ',version
      WRITE(*,'(A)') ' ENTER Control File Name : '
C     Control file name is usually "????????.CON"
      READ (*,'(A)') contfile
      filename=contfile
      i=index(filename,'.')
      j=index(filename,' ')
      if ((i.gt.9).or.
     >    ((i.eq.0).and.(j.ge.12)).or.
     >    (j-i.gt.4)) then
         write(*,*) filename,' --- name is too long !'
         stop
      endif
      OPEN (9,FILE=filename)
      read(9,'(A80)') a80
      read(a80,'(I2)') idebug
      if (idebug.lt.10) then
         nfile=1
         files(1)=contfile
         iauto=0
      else
C        read many control file names if idebug>=10
         nfile=0
         iauto=1
         do ifile=1,nfilelimit
            read(9,'(A12)',END=5000) filename
            i=index(filename,'.')
            j=index(filename,' ')
            if ((i.gt.9).or.
     >          ((i.eq.0).and.(j.ge.12)).or.
     >          (j-i.gt.4)) then
               write(*,*) filename,' --- name is too long !'
               stop
            endif
            nfile=nfile+1
            files(nfile)=filename
         enddo
 5000    continue
      endif
      close(9)
C
C    Loop for different control files
C
      do 5001 ifile=1,nfile
      istopflag=0
      contfile=files(ifile)
      filename=contfile
C
      call readcontfile(filename)
C
      DLOGP=DLOG10(PTOT)
      DO 11 J=1,NS+msExtra
      dJJ=0d0
      DO 45 K=1,M
 45   dJJ=dJJ+dNU(J,K)
C     XC1(): coeff. of logK for the use of variables (partial pressure/ptot)
 11   XC1(J)=-XC(1,J)-DLOGP*dJJ
      DO 24 I=1,N
      XC1(I)=XC1(I)+DLOGP
 24   NQ(I)=I
C     write calculation parameters
      if ((imode.eq.0).or.(imode.eq.2)) then
         write(8,300) dtmin
         if (isummary.ge.1) write(7,300) dtmin
         print 301,dtmin
 300  format('Temperature resolution was set to ',F5.2,' K')
 301  format(' Temperature resolution was set to ',F5.2,' K')
      elseif ((imode.eq.1).or.(imode.eq.3)) then
        if (ibar.eq.0) then
         write(8,302) dpmin
         if (isummary.ge.1) write(7,302) dpmin
         print 303,dpmin
        else
         write(8,1302) dpmin
         if (isummary.ge.1) write(7,1302) dpmin
         print 1303,dpmin
        endif
 302  format('Pressure resolution was set to ',F5.2,' log10(atm)')
 303  format(' Pressure resolution was set to ',F5.2,' log10(atm)')
1302  format('Pressure resolution was set to ',F5.2,' log10(bar)')
1303  format(' Pressure resolution was set to ',F5.2,' log10(bar)')
      endif
      if (ibar.eq.0) then
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(8,101)PTOT
        if (isummary.ge.1) write(7,101) ptot
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(8,99)PTOT
        if (isummary.ge.1) write(7,99) ptot
      endif
      else
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(8,1101)PTOT
        if (isummary.ge.1) write(7,1101) ptot
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(8,199)PTOT
        if (isummary.ge.1) write(7,199) ptot
      endif
      endif
 101  FORMAT(/,'P initial =',1PE10.3,' atm ',
     >'(if H2 is dominant, the actual pressure is half)',
     >/,'Atomic Abundances relative to their sum = 1 :')
 99   FORMAT(/,'P initial =',1PE10.3,' atm ',
     >/,'Atomic Abundances relative to their sum = 1 :')
1101  FORMAT(/,'P initial =',1PE10.3,' bar ',
     >'(if H2 is dominant, the actual pressure is half)',
     >/,'Atomic Abundances relative to their sum = 1 :')
199   FORMAT(/,'P initial =',1PE10.3,' bar ',
     >/,'Atomic Abundances relative to their sum = 1 :')
      WRITE(8,49)(WE(I),R(I),I=1,M)
      if (isummary.ge.1) write(7,49)(WE(I),R(I),I=1,M)
 49   format(5(1x,a2,1pe12.5,:,1x))
      if (ibar.eq.0) then
         write(8,48)
      else
         write(8,148)
      endif
      if (isummary.ge.1) write(7,'()')
 48   format(/,'To calculate mol/l of any species, multiply ',
     >'(N total) times (Fraction) ',/,'divided by ',
     >'(stoichiometric coeff.). To calculate the partial ',
     >'pressure (atm), ',/,'multiply mol/l times RT.',/)
148   format(/,'To calculate mol/l of any species, multiply ',
     >'(N total) times (Fraction) ',/,'divided by ',
     >'(stoichiometric coeff.). To calculate the partial ',
     >'pressure (bar), ',/,'multiply mol/l times RT.',/)
      if (itable.ge.4) call openwttables()
      if (idebug.ge.3) then
C        write initial values of variables and stop
         write(8,3014) m,n,ngp,ns,msExtra,nss
         do i=1,m
         write(8,3011) i,wee(i),r(i),p(i)
         enddo
         if (nss.ne.0) then
            do i=1,nss
            write(8,3013) i,wss(i),ms(i),mdb(i),mse(i),nssp(i),mss(i)
            enddo
         endif
         if (is.ne.0) then
            do i=n1,nis
            write(8,*)nq(i),p(i+m),wc(nq(i))
            enddo
         endif
         do i=1,ns+msExtra
         print 3011,i,wc(i),(xc(j,i),j=1,iord),p(m+i)
         print 3012,nspecies(i),(nu(i,j),j=1,m)
         write(8,3011) i,wc(i),(xc(j,i),j=1,iord),p(m+i)
         write(8,3012) nspecies(i),(nu(i,j),j=1,m)
         enddo
         stop
 3011 format(1X,I3,1X,A8,6(1pe11.4,:))
 3012 format(4X,I3,25(I3,:))
 3013 format(1X,I3,1X,A8,6(I4,:))
 3014 format(1X,'m=',I2,'  n=',I3,'  ngp=',I3,'  ns=',I3,
     >'  msExtra=',I2,'  nss=',I2)
      endif
C
C   Main Loop
C
C     Calculate equilibrium composition from Tstart to Tstop by DTemp step.
C     T0 is a current temperature to be checked.
C     If new condensates appear or old ones disappear, seek a temperature
C     at which only one such event is found such that
C     the difference from the temperature where no new condensates are found
C     becomes less than the temperature resolution specified.
C
C     Temperature variable
C
      if ((imode.eq.0).or.(imode.eq.2)) then
          t0=tstart
          dt=dtemp
          Tprevious=t0+dt
          Tdone=t0+dt
          Ttobe=t0
          Tnew=Ttobe
          imustsee=0
          istart=0
          call psave(1)
 201  dt=Tdone-t0
      if (dt.eq.0d0) dt=dtemp
      ifailconv=0
      if (t0.eq.Ttobe) then
         temp=0d0
         print 100,T0,Tdone
         if (idebug.ge.2) write(8,100) t0,Tdone
         if (isummary.ge.2) write(15,100) t0,Tdone
         if (iwindow.ge.1) write(16,100) t0,Tdone
         if (iwindow.ge.1) write(18,120) t0,Ttobe,Tprevious,temp,Tdone
 100     format(' Checking',F9.3,' K  ( Done',F9.3,' K )')
      else
         print 110,T0,Tdone,Tnew
         if (idebug.ge.2) write(8,110) t0,Tdone,Tnew
         if (isummary.ge.2) write(15,110) t0,Tdone,Tnew
         if (iwindow.ge.1) write(16,110) t0,Tdone,Tnew
         if (iwindow.ge.1) write(18,120) t0,Ttobe,Tprevious,Tnew,Tdone
 110     format(' Checking',F9.3,' K  ( Done',F9.3,' K;',
     >          ' New condns./disap. found',F9.3,' K )')
 120     format(' Current Temp.  =',F9.3,' K',/,
     >          ' T to be/T prev.=',F9.3,' /',F9.3,/,
     >          ' Cond. Bracket  =',F9.3,' /',F9.3)
      endif
C     check # of new condensates that appeared and old ones that disappeared
C     at this temperature t0
C
      call checktemp()
C
      if (ifailconv.ne.0) then
          if ((istopflag.eq.1).or.(t0.eq.tstart)) then
             goto 1000
          elseif (dt.lt.dtmin) then
             print *,' Failed to converge !!!'
             goto 1000
          else
             t0=(t0+Tdone)/2d0
             goto 201
          endif
      endif
          if ((icps.eq.0).and.(icpd.eq.0).and.(icss.eq.0)) then
                if (t0.eq.Ttobe) then
C                  Nothing happened, move to the next temperature
                        Tdone=Ttobe
                        call prestore(2)
                        call psave(1)
                        call printout()
                        Tprevious=Ttobe
                        Ttobe=Ttobe-dtemp
                        t0=Ttobe
                        Tnew=Ttobe
                elseif (t0-Tnew.lt.dtmin) then
                        Tdone=t0
                        call prestore(2)
                        call psave(1)
                        t0=Tnew
                   print *,'No new condensates/disappearance found, ',
     >                     'seeking the temperature'
                   if (isummary.ge.2) write(15,*)
     >                  'No new condensates - seeking the temperature'
                   if (iwindow.ge.1) write(16,*)
     >                  'No new condensates - seeking the temperature'
                   if (iwindow.ge.1) write(18,*)
     >                  'No new condensates - seeking the temperature'
                else
                        Tdone=t0
                        call prestore(2)
                        call psave(1)
                        t0=(Tdone+Tnew)/2d0
                   print *,'No new condensates/disappearance found, ',
     >                     'seeking the temperature'
                   if (isummary.ge.2) write(15,*)
     >                  'No new condensates - seeking the temperature'
                   if (iwindow.ge.1) write(16,*)
     >                  'No new condensates - seeking the temperature'
                   if (iwindow.ge.1) write(18,*)
     >                  'No new condensates - seeking the temperature'
                endif
          elseif (t0.eq.tstart) then
                print *,' New condensate(s) at starting temperature !'
                   if (iwindow.ge.1) write(18,*)
     >                  ' New condensate(s) at starting temperature !'
                istart=istart+1
                call prestore(2)
                call psave(1)
                if (istart.le.20) goto 201
                goto 1000
          elseif ((icpd.le.1).and.(icps+icss.le.1)) then
                if (Tdone-t0.lt.dtmin) then
                   print *,'1 new condensate/disappearance found ',
     >                     '--- done'
                   if (isummary.ge.2) write(15,*)
     >                     '1 new condensate --- done'
                   if (iwindow.ge.1) write(16,*)
     >                     '1 new condensate --- done'
                   if (iwindow.ge.1) write(18,*)
     >                     '1 new condensate --- done'
                   if ((Tprevious.ne.Tdone).and.(imustsee.eq.0)) then
C                        Output results of the previous temperature
                            t0t=t0
                            icpst=icps
                            icpdt=icpd
                            icsst=icss
                            t0=Tdone
                            icps=0
                            icpd=0
                            icss=0
                         call prestore(1)
                         call printout()
                            t0=t0t
                            icps=icpst
                            icpd=icpdt
                            icss=icsst
                   endif
                   if ((t0.eq.Tdone).and.(imustsee.ne.0)) then
                      print 501,wss(imustsee)
                      write (8,501)wss(imustsee)
                      if (isummary.ge.1) write(7,501)wss(imustsee)
                      if (isummary.ge.2) write(15,501)wss(imustsee)
                      if (iwindow.ge.1) write(16,501)wss(imustsee)
                      if (iwindow.ge.1) write(17,501)wss(imustsee)
                      if (iwindow.ge.1) write(18,501)wss(imustsee)
 501                  format(' ***WARNING*** Condensation of ',A8,
     >                    ' must be earlier than this step!')
                      call unsetmustsee(imustsee)
                      imustsee=0
                   elseif (imustsee.ne.0) then
                      call unsetmustsee(imustsee)
                      imustsee=0
                   elseif (imustsee.eq.0) then
                    if (icss.eq.1) then
                      if ((nonid(isss(1)).ne.0).or.
     >                    (isss(1).eq.nssspinel)) then
                         imustsee=isss(1)
                         print 500,wss(imustsee)
                      if (isummary.ge.2) write(15,500)wss(imustsee)
                      if (iwindow.ge.1) write(16,500)wss(imustsee)
                      if (iwindow.ge.1) write(18,500)wss(imustsee)
 500                     format(' Make sure that ',A8,
     >                     ' is not stable in the previous step')
                         call prestore(2)
                         call setmustsee(imustsee)
                         Tnew=t0
                         t0=Tdone
                         goto 201
                      endif
                    endif
                   endif
                   if (t0.eq.Ttobe) then
                         Tdone=Ttobe
                         call prestore(2)
                         call psave(1)
                         call printout()
                         Tprevious=Ttobe
                         Ttobe=Ttobe-dtemp
                         t0=Ttobe
                         Tnew=Ttobe
                   else
                         Tdone=t0
                         Tprevious=Tdone
                         call prestore(2)
                         call psave(1)
                         call printout()
                         t0=Ttobe
                         Tnew=Ttobe
                   endif
                else
                   Tnew=t0
                   t0=(Tdone+Tnew)/2d0
                   print *,'1 new condensate/disappearance found, ',
     >                     'seeking the temperature'
                   if (isummary.ge.2) write(15,*)
     >                     '1 new condensate - seeking the temperature'
                   if (iwindow.ge.1) write(16,*)
     >                     '1 new condensate - seeking the temperature'
                   if (iwindow.ge.1) write(18,*)
     >                     '1 new condensate - seeking the temperature'
                endif
          else
                Tnew=t0
                t0=(Tdone+Tnew)/2d0
                print *,'2 or more new condensates/disappearance found'
                   if (isummary.ge.2) write(15,*)
     >                  '2 or more new condensates/disappearance found'
                   if (iwindow.ge.1) write(16,*)
     >                  '2 or more new condensates/disappearance found'
                   if (iwindow.ge.1) write(18,*)
     >                  '2 or more new condensates'
          endif
      IF(T0.LT.(TSTOP-1d-10)) GO TO 1000
      GO TO 201
C
C     Pressure variable
C
      elseif ((imode.eq.1).or.(imode.eq.3)) then
          p0=ptot
          dp=dpres
          Pprevious=10**(dlog10(p0)-dp)
          Pdone=Pprevious
          Ptobe=p0
          Pnew=Ptobe
          call psave(1)
          istart=0
          imustsee=0
 202  continue
          ptot=p0
          DLOGP=DLOG10(PTOT)
          do J=1,NS+msExtra
             dJJ=0
             do K=1,M
                dJJ=dJJ+dNU(J,K)
             enddo
             XC1(J)=-XC(1,J)-DLOGP*dJJ
          enddo
          do I=1,N
             XC1(I)=XC1(I)+DLOGP
          enddo
      dp=-Pdone+p0
      dp2=-dlog10(Pdone)+dlog10(p0)
      if (dp.eq.0d0) dp=-10**(dlog10(p0)-dpres)+p0
      if (dp2.eq.0d0) dp2=dpres
      ifailconv=0
      if (p0.eq.Ptobe) then
         temp=0d0
        if (ibar.eq.0) then
         print 102,p0,Pdone
         if (idebug.ge.2) write(8,102) p0,Pdone
         if (isummary.ge.2) write(15,102) p0,Pdone
         if (iwindow.ge.1) write(16,102) p0,Pdone
         if (iwindow.ge.1) write(18,121) p0,Pprevious,Ptobe,temp,Pdone
        else
         print 102,p0,Pdone
         if (idebug.ge.2) write(8,1102) p0,Pdone
         if (isummary.ge.2) write(15,1102) p0,Pdone
         if (iwindow.ge.1) write(16,1102) p0,Pdone
         if (iwindow.ge.1) write(18,1121) p0,Pprevious,Ptobe,temp,Pdone
        endif
 102     format(' Checking',1pE10.3,' atm  ( Done',E10.3,' atm )')
1102     format(' Checking',1pE10.3,' bar  ( Done',E10.3,' bar )')
      else
        if (ibar.eq.0) then
         print 111,p0,Pdone,Pnew
         if (idebug.ge.2) write(8,111) p0,Pdone,Pnew
         if (isummary.ge.2) write(15,111) p0,Pdone,Pnew
         if (iwindow.ge.1) write(16,111) p0,Pdone,Pnew
         if (iwindow.ge.1) write(18,121) p0,Pprevious,Ptobe,Pdone,Pnew
        else
         print 1111,p0,Pdone,Pnew
         if (idebug.ge.2) write(8,1111) p0,Pdone,Pnew
         if (isummary.ge.2) write(15,1111) p0,Pdone,Pnew
         if (iwindow.ge.1) write(16,1111) p0,Pdone,Pnew
         if (iwindow.ge.1) write(18,1121) p0,Pprevious,Ptobe,Pdone,Pnew
        endif
 111     format(' Checking',1pE10.3,' atm  ( Done',E10.3,' atm;',
     >          ' New condns./disap. found',E10.3,' atm )')
 121     format(' Current Pres.  =',1pE10.3,' atm',/,
     >          ' Previous/To be =',E10.3,' /',E10.3,/,
     >          ' Cond. Bracket  =',E10.3,' /',E10.3)
1111     format(' Checking',1pE10.3,' bar  ( Done',E10.3,' bar;',
     >          ' New condns./disap. found',E10.3,' bar )')
1121     format(' Current Pres.  =',1pE10.3,' bar',/,
     >          ' Previous/To be =',E10.3,' /',E10.3,/,
     >          ' Cond. Bracket  =',E10.3,' /',E10.3)
      endif
C     check # of new condensates that appeared and old ones that disappeared
C     at this pressure p0
C
      call checktemp()
C
      if (ifailconv.ne.0) then
          if ((istopflag.eq.1).or.(p0.eq.tstart)) then
             goto 1000
          elseif (dp2.lt.dpmin) then
             print *,' Failed to converge !!!'
             goto 1000
          else
             p0=10**((dlog10(p0)+dlog10(Pdone))/2d0)
             goto 202
          endif
      endif
          if ((icps.eq.0).and.(icpd.eq.0).and.(icss.eq.0)) then
                if (p0.eq.Ptobe) then
C                  Nothing happened, move to the next pressure
                        Pdone=Ptobe
                        call prestore(2)
                        call psave(1)
                        call printout()
                        Pprevious=Ptobe
                        Ptobe=10**(dlog10(Ptobe)+dpres)
                        p0=Ptobe
                        Pnew=Ptobe
                elseif (-dlog10(p0)+dlog10(Pnew).lt.dpmin) then
                        Pdone=p0
                        call prestore(2)
                        call psave(1)
                        p0=Pnew
                   print *,'No new condensates/disappearance found, ',
     >                     'seeking the pressure'
                   if (isummary.ge.2) write(15,*)
     >                     'No new condensates - seeking the pressure'
                   if (iwindow.ge.1) write(16,*)
     >                     'No new condensates - seeking the pressure'
                   if (iwindow.ge.1) write(18,*)
     >                     'No new condensates - seeking the pressure'
                else
                        Pdone=p0
                        call prestore(2)
                        call psave(1)
                        p0=10**((dlog10(Pnew)+dlog10(Pdone))/2d0)
                   print *,'No new condensates/disappearance found, ',
     >                     'seeking the pressure'
                   if (isummary.ge.2) write(15,*)
     >                     'No new condensates - seeking the pressure'
                   if (iwindow.ge.1) write(16,*)
     >                     'No new condensates - seeking the pressure'
                   if (iwindow.ge.1) write(18,*)
     >                     'No new condensates - seeking the pressure'
                endif
          elseif (p0.eq.pstart) then
                print *,' New condensate(s) at starting pressure !'
                   if (isummary.ge.2) write(15,*)
     >                    ' New condensate(s) at starting pressure !'
                   if (iwindow.ge.1) write(16,*)
     >                    ' New condensate(s) at starting pressure !'
                   if (iwindow.ge.1) write(18,*)
     >                    ' New condensate(s) at starting pressure !'
                istart=istart+1
                call prestore(2)
                call psave(1)
                if (istart.le.20) goto 202
                goto 1000
          elseif ((icpd.le.1).and.(icps+icss.le.1)) then
                if (-dlog10(Pdone)+dlog10(p0).lt.dpmin) then
                   print *,'1 new condensate/disappearance found ',
     >                     '--- done'
                   if (isummary.ge.2) write(15,*)
     >                     '1 new condensate --- done'
                   if (iwindow.ge.1) write(16,*)
     >                     '1 new condensate --- done'
                   if (iwindow.ge.1) write(18,*)
     >                     '1 new condensate --- done'
                   if ((Pprevious.ne.Pdone).and.(imustsee.eq.0)) then
C                        Output results of the previous pressure
                            p0t=p0
                            icpst=icps
                            icpdt=icpd
                            icsst=icss
                            p0=Pdone
                            ptot=Pdone
                            icps=0
                            icpd=0
                            icss=0
                         call prestore(1)
                         call printout()
                            p0=p0t
                            ptot=p0t
                            icps=icpst
                            icpd=icpdt
                            icss=icsst
                   endif
                   if ((p0.eq.Pdone).and.(imustsee.ne.0)) then
                      print 501,wss(imustsee)
                      write (8,501)wss(imustsee)
                      if (isummary.ge.1) write(7,501)wss(imustsee)
                      if (isummary.ge.2) write(15,501)wss(imustsee)
                      if (iwindow.ge.1) write(16,501)wss(imustsee)
                      if (iwindow.ge.1) write(17,501)wss(imustsee)
                      if (iwindow.ge.1) write(18,501)wss(imustsee)
c 501                  format(' ***WARNING*** Condensation of ',A8,
c     >                    ' must be earlier than this step!')
                      call unsetmustsee(imustsee)
                      imustsee=0
                   elseif (imustsee.ne.0) then
                      call unsetmustsee(imustsee)
                      imustsee=0
                   elseif (imustsee.eq.0) then
                    if (icss.eq.1) then
                      if ((nonid(isss(1)).ne.0).or.
     >                    (isss(1).eq.nssspinel)) then
                         imustsee=isss(1)
                         print 500,wss(imustsee)
                      if (isummary.ge.2) write(15,500)wss(imustsee)
                      if (iwindow.ge.1) write(16,500)wss(imustsee)
                      if (iwindow.ge.1) write(18,500)wss(imustsee)
c 500                     format(' Make sure that ',A8,
c     >                     ' is not stable in the previous step')
                         call prestore(2)
                         call setmustsee(imustsee)
                         Pnew=p0
                         p0=Pdone
                         goto 202
                      endif
                    endif
                   endif
                   if (p0.eq.Ptobe) then
                         Pdone=Ptobe
                         call prestore(2)
                         call psave(1)
                         call printout()
                         Pprevious=Ptobe
                         Ptobe=10**(dlog10(Ptobe)+dpres)
                         p0=Ptobe
                         Pnew=Ptobe
                   else
                         Pdone=p0
                         Pprevious=Pdone
                         call prestore(2)
                         call psave(1)
                         call printout()
                         p0=Ptobe
                         Pnew=Ptobe
                   endif
                else
                   Pnew=p0
                   p0=10**((dlog10(Pnew)+dlog10(Pdone))/2d0)
                   print *,'1 new condensate/disappearance found, ',
     >                     'seeking the pressure'
                   if (isummary.ge.2) write(15,*)
     >                    '1 new condensate - seeking the pressure'
                   if (iwindow.ge.1) write(16,*)
     >                    '1 new condensate - seeking the pressure'
                   if (iwindow.ge.1) write(18,*)
     >                    '1 new condensate - seeking the pressure'
                endif
          else
                Pnew=p0
                p0=10**((dlog10(Pnew)+dlog10(Pdone))/2d0)
                print *,'2 or more new condensates/disappearance found'
                   if (isummary.ge.2) write(15,*)
     >                 '2 or more new condensates/disappearance found'
                   if (iwindow.ge.1) write(16,*)
     >                 '2 or more new condensates/disappearance found'
                   if (iwindow.ge.1) write(18,*)
     >                 '2 or more new condensates'
          endif
      IF(p0.GT.(Pstop+1d-10)) GO TO 1000
      GO TO 202
      endif
C     End of main loop
 1000 CONTINUE
      close(8)
      if(isummary.ge.1) close(7)
      if(isummary.ge.2) close(15)
      if(itable.ge.4) close(3)
      if(iwindow.ge.1) close(16)
      if(iwindow.ge.1) close(17)
      if(iwindow.ge.1) close(18)
      if(itable.ge.1) call tables()
C     End of the calculation for each control file
 5001 continue
      end

      subroutine checktemp()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.for'
      dimension sxtemp(isslimit)
C     Check # of new condensates that appeared and old ones that disappeared
C     at this temperature t0.
C     This routine sets up the condensates and calls subroutine temperature();
C     first with the same set of condensates as the previous temperature,
C     then with each solid solution added to the set if that solid solution
C     is about to condense, since the temperature() subroutine searches only
C     for new pure solids and we are not sure whether the solid solution
C     actually condenses or not until solving the exact solution with it.
C     Combining both results, return:
C            icps: # of new pure solids
C            icss: # of new solid solutions
C            icpd: # of phases that disappeared
      eps2=dsqrt(eps)
      call prestore(1)
      icps=0
      icpd=0
      icss=0
c      print 100,T0
c      if (idebug.ge.2) write(8,100) t0
c 100  format(' Checking temperature :',F10.4,' K')
      IC=0
      ID=0
      ICONT=0
c**   replace pure solid with one of polymorphs if it's more stable
      call polymorphs()
      icps=icpst
      icpd=icpdt
      if(icpst.eq.0) goto 130
      do 131 ii=1,icpst
 131     isps(ii)=ispst(ii)
 130  if(icpdt.eq.0) goto 132
      do 133 ii=1,icpdt
 133     ispd(ii)=ispdt(ii)
 132  continue
      isave2flag=0
C
      call temperature()
C
C      print *,icpst,icpdt,icsst
      if (ifailconv.ne.0) return
      if (idebug.ge.1) call freeenergy(G)
      if(icpst.eq.0) goto 30
      do 31 ii=1,icpst
 31      isps(ii+icps)=ispst(ii)
 30   if(icpdt.eq.0) goto 32
      do 33 ii=1,icpdt
 33      ispd(ii+icpd)=ispdt(ii)
 32   if(icsst.eq.0) goto 34
      do 35 ii=1,icsst
 35      isss(ii+icss)=issst(ii)
 34   continue
      icps=icps+icpst
      icpd=icpd+icpdt
      icss=icss+icsst
      call packghost()
      call psave(2)
C
c**   check stability of solid solutions
C
      if(nssorg.eq.0) go to 11
      do 10 j2=1,nssorg
      call prestore(3)
      j=j2
      ighflag=0
      do jj=1,ms(j)
         sxtemp(jj)=0d0
      enddo
      if (MSS(J2).GT.0) then
         if ((nonid(j2).eq.0).or.(j2.eq.nssspinel)) then
            GO TO 10
         else
            iflag=0
            if (nss.gt.nssorg) then
               do jj=nssorg+1,nss
                  if ((nssghost(jj).eq.j2).and.(mss(jj).eq.0)) then
                     iflag=jj
                     goto 51
                  endif
               enddo
  51           continue
            endif
           if (iflag.ne.0) then
C           use an old ghost
            j=iflag
            do jj=1,ms(j)
               do k=1,nseed(j)
                  sxold(j,jj,k)=sxold(j2,jj,k)
               enddo
            enddo
           else
C           set up a new ghost solid solution
            ighflag=1
            j=nss+1
            if (j.gt.nsslimit) then
               print *,' Exceed NSS limit !!!'
               stop
            endif
            nss=j
            nssp(j)=ns+msExtra+1
            nssghost(j)=j2
            wss(j)=wss(j2)
            ms(j)=ms(j2)
            mss(j)=0
            mdb(j)=mdb(j2)
            mse(j)=mse(j2)
            nonid(j)=nonid(j2)
            nseed(j)=nseed(j2)
            amusstotal(j)=amusstotal(j2)
            if (nssp(j)+ms(j2)+mse(j2).gt.nslimit) then
               print *,' Exceed NS limit !!!'
               stop
            endif
            ns=ns+ms(j2)
            msExtra=msExtra+mse(j2)
            do jj=1,ms(j)
               wc(nssp(j)+jj-1)=wc(nssp(j2)+jj-1)
               wref(nssp(j)+jj-1)=wref(nssp(j2)+jj-1)
               nspecies(nssp(j)+jj-1)=nspecies(nssp(j2)+jj-1)
               kind(nssp(j)+jj-1)=kind(nssp(j2)+jj-1)
               ag0(m+nssp(j)+jj-1)=ag0(m+nssp(j2)+jj-1)
               amu(m+nssp(j)+jj-1)=amu(m+nssp(j2)+jj-1)
               amuss(j,jj)=amuss(j2,jj)
               do k=1,m
                  nu(nssp(j)+jj-1,k)=nu(nssp(j2)+jj-1,k)
                  dnu(nssp(j)+jj-1,k)=dnu(nssp(j2)+jj-1,k)
               enddo
               do k=1,m2wt
                  nu2(nssp(j)+jj-1,k)=nu2(nssp(j2)+jj-1,k)
               enddo
               do k=1,iord
                  xc(k,nssp(j)+jj-1)=xc(k,nssp(j2)+jj-1)
               enddo
                  xc1(nssp(j)+jj-1)=xc1(nssp(j2)+jj-1)
                  xk(nssp(j)+jj-1)=xk(nssp(j2)+jj-1)
               do k=1,nseed(j)
                  sxold(j,jj,k)=sxold(j2,jj,k)
                  seed(j,jj,k)=seed(j2,jj,k)
C                  mustsee(j,k)=mustsee(j2,k)
               enddo
               do k=1,nfeas
                  feas(j,jj,k)=feas(j2,jj,k)
               enddo
            enddo
           endif
         endif
      endif
      do 12 iseed=1,nseed(j)
c      print *,j2,j,iseed,nssp(j2),nssp(j)
      iflag=0
      ipoly=0
      I=nssp(j)
      ibase=nssp(j)+ms(j)-1
      if((icps+icss.gt.1).or.(icpd.gt.1)) goto 11
      call prestore(3)
c**   search candidates of solid solutions
      if (icss.eq.0) goto 41
      do 40 ii=1,icss
      if ((nonid(j).eq.0).and.(j.eq.isss(ii))) goto 10
 40   continue
 41   if (icpd.eq.0) goto 43
      do 42 ii=1,icpd
      if ((nonid(j).eq.0).and.(ispd(ii).eq.nssp(j)+ms(j)-1)) goto 10
 42   continue
 43   KCH=MS(j)
C      if ((wss(j)(1:8).eq.'Spinel  ').and.(nonid(j).eq.0)) then
Cc**   Spinel --- check only MgAl2O4 (old way)
CC     If Sum((stoichiom. coeff.)x(logP(element)))-logK > 0 (actually > -0.01),
CC     then check the stability of spinel.
C        ibase=nssp(nssspinel)+ispinel(4)-1
C        P0=-XK(IBASE)
C        DO 89 JJ=1,M
C 89     P0=P0+DLP(JJ)*dnu(IBASE,JJ)
C        if (idebug.ge.1) print 510,wss(j),p0
C        if (idebug.ge.2) write(8,510)wss(j),p0
C        IF(P0.GT.-0.01) iflag=1
C      else
C
C     Solid solutions
C     If sum of mole fractions of endmembers > 1 (actually sscheck),
C     then check the stability.
       if (iseed.eq.1) then
        if (wss(j)(1:8).eq.'Spinel  ') then
C         Spinel solid solution
C         Ideal estimation of X(Cr) and X(Fe)
          do ii=1,4
            sx(j,ii)=-xk(I+ii-1)
            DO JJ=1,M
              sx(j,ii)=sx(j,ii)+DLP(JJ)*dnu(I+ii-1,JJ)
            end do
          end do
          xcrt=sx(j,ispinel(4))-sx(j,ispinel(2))
          xfet=sx(j,ispinel(4))-sx(j,ispinel(3))
          xcrt=1/(1+10d0**(xcrt/2d0))
          xfet=1/(1+10d0**(xfet))
          sx(j,1)=xcrt
          sx(j,2)=xfet
          sx(j,3)=0d0
        else
c**       ideal solid solutions
          do ii=1,kch
            sx(j,ii)=-xk(I+ii-1)
            DO JJ=1,M
              sx(j,ii)=sx(j,ii)+DLP(JJ)*dnu(I+ii-1,JJ)
            end do
            sx(j,ii)=10d0**sx(j,ii)
          end do
c          print *,(sx(j,ii),ii=1,kch)
        endif
       endif
c**     non-ideal solid solutions
C       Subroutine "sxnonideal" searches the miniumum free energy
C       for the solution and compare it with the current phase assemblage.
C       If the solution is more stable, it returns calculated mole fractions.
C       If not, all mole fractions are set to be zero.
        iseed2=iseed
        if ((nonid(j).ne.0).or.(wss(j)(1:8).eq.'Spinel  ')) then
          call sxnonideal(j,iseed2)
        endif
        sxt=0d0
        do ii=1,kch
          sxt=sxt+sx(j,ii)
        end do
       if (nonid(j).eq.0) then
        if (idebug.ge.1) print 510,wss(j),sxt,(sx(j,ii),ii=1,kch)
        if (idebug.ge.2) write(8,510)wss(j),sxt,(sx(j,ii),ii=1,kch)
       endif
        if (j.gt.nssorg) then
           do j3=1,j-1
              if (((nssghost(j).eq.nssghost(j3)).or.
     >            (nssghost(j).eq.j3)).and.(mss(j3).ne.0)) then
                 iflag2=0
                 do ii=1,kch
                    if (dabs(sx(j,ii)-sx(j3,ii)).gt.eps2) iflag2=1
                 enddo
                 if (iflag2.eq.0) then
                    if (idebug.ge.1) print 511,wss(j)
                    if (idebug.ge.1) write(8,511)wss(j)
                    if (isummary.ge.2) write(15,511)wss(j)
                    if (iwindow.ge.1) write(16,511)wss(j)
                    sxt=0d0
                 endif
              endif
           enddo
        endif
        sxtcheck=sscheck
        if (nonid(j).ne.0) sxtcheck=ss2check
        if (sxt.gt.sxtcheck-eps2) iflag=1
C      endif
 510  format(1x,A8,' total=',1pe9.2,5(1x,0pf7.4))
 511  format(1x,A8,' Identical composition exists')
      call sspolymorphs(j,ipoly)
      if ((nonid(j).eq.0).and.(ipoly.ne.0).and.(iflag.eq.1)) then
c**   replace solid solution with one of polymophs if it's more stable
        call ssfreeenergy(ipoly,G)
        Gcurrent=G
        FAC=1
        DO 165 ii=1,M
        IF(NU(ibase,ii)) 163,165,163
 163     FAC=DMIN1(R(ii)/dnu(ibase,ii),FAC)
 165     CONTINUE
        MSS(j)=mss(ipoly)
        mss(ipoly)=0
        DO 182 ii=1,KCH
c        if (j.eq.nssspinel) then
c           if ((ii.eq.1).or.(ii.eq.2)) then
c              p(mss(j)+ii-1+n)=factor
c           else
c              p(mss(j)+ii-1+n)=factor*fac
c           endif
c        else
c          P(mss(j)+ii-1+n)=factor*FAC*sx(j,ii)
c        endif
 182    NQ(mss(j)+ii-1-m+n)=i+ii-1
        xk0(i-n)=0d0
        xk1(i-n)=0d0
        print 1101,wss(j),wss(ipoly)
        if (idebug.ge.2) write(8,1101)wss(j),wss(ipoly)
 1101   format('   --- checking stability of ',A8,' replacing ',A8)
        IC=0
        ID=0
        ICONT=0
C
        call temperature()
C
        if (ifailconv.ne.0) return
        if (mss(j).gt.0) then
          call ssfreeenergy(j,G)
          if (G.lt.Gcurrent) then
            print 1102,wss(j),wss(ipoly)
            if (idebug.ge.2) write(8,1102)wss(j),wss(ipoly)
 1102       format('   --- ',A8,' is more stable than ',A8,' !')
            icss=icss+1
            isss(icss)=j
            icpdt=icpdt+1
            ispdt(icpdt)=nssp(ipoly)+ms(ipoly)-1
            if (icpdt.eq.0) goto 120
              iii=icpd
              do 122 ii=1,icpdt
              if (iii.eq.0) goto 124
              do 123 jj=1,iii
                 if (ispd(jj).eq.ispdt(ii)) goto 122
  123         continue
  124         icpd=icpd+1
              ispd(icpd)=ispdt(ii)
  122         continue
  120       if (icpst.eq.0) goto 121
              iii=icps
              do 125 ii=1,icpst
              if (iii.eq.0) goto 126
              do 127 jj=1,iii
                 if (isps(jj).eq.ispst(ii)) goto 125
  127         continue
  126         icps=icps+1
              isps(icps)=ispst(ii)
  125         continue
  121       continue
            call packghost()
            call psave(2)
          else
            print 1103,wss(j),wss(ipoly)
            if (idebug.ge.2) write(8,1103)wss(j),wss(ipoly)
 1103       format('   --- ',A8,' is not more stable than ',A8)
          endif
        else
          print 103,wss(j)
          if (idebug.ge.2) write(8,103)wss(j)
        endif
      elseif ((iflag.eq.1).and.((nonid(j).ne.0).or.(ipoly.eq.0))) then
c**   set a solid solution as a condensate if it may condense
        FAC=1
        DO 65 ii=1,M
        IF(NU(ibase,ii)) 63,65,63
 63     FAC=DMIN1(R(ii)/dnu(ibase,ii),FAC)
 65     CONTINUE
C       set the position of solid solution
        MSS(j)=IS+m+1
 81     DO 82 ii=1,KCH
        NIS=NIS+1
        if (nis-n.gt.islimit) then
         print *,' Exceed IS limit !!!'
         stop
        endif
C       set initial values for P
        if (j.eq.nssspinel) then
           if (ii.eq.1) then
              p(m+nis)=sx(j,1)
           elseif (ii.eq.2) then
              p(m+nis)=sx(j,2)
           else
              p(m+nis)=factor*fac
           endif
        else
          P(M+NIS)=factor*FAC*sx(j,ii)
        endif
C       set the species # (internal, not as originally fed in)
 82     NQ(NIS)=I+ii-1
        IS=KCH+IS
        if (is.gt.islimit) then
         print *,' Exceed IS limit !!!'
         stop
        endif
        xk0(i-n)=0d0
        xk1(i-n)=0d0
        print 101,wss(j),iseed
        if (idebug.ge.2) write(8,101)wss(j),iseed
        if (isummary.ge.2) write(15,101)wss(j),iseed
        if (iwindow.ge.1) write(16,101)wss(j),iseed
 101    format('   --- checking stability of ',A8,'  Seed#=',I3)
        IC=0
        ID=0
        ICONT=0
C
        call temperature()
C
        if (ifailconv.ne.0) return
C       check whether the solid solution is still among the stable condensates
        if (mss(j).gt.0) then
          print 102,wss(j)
          if (idebug.ge.2) write(8,102)wss(j)
          if (isummary.ge.2) write(15,102)wss(j)
          if (iwindow.ge.1) write(16,102)wss(j)
 102      format('   --- ',A8,' is stable !')
         iflag2=1
         if (iseed.gt.1) then
          iflag2=0
          do ii=1,kch
             if (dabs(sx(j,ii)-sxtemp(ii)).gt.eps2) iflag2=1
          enddo
         endif
         if (iflag2.eq.0) then
          if (idebug.ge.0) print 512,wss(j)
          if (idebug.ge.1) write(8,512)wss(j)
          if (isummary.ge.2) write(15,512)wss(j)
          if (iwindow.ge.1) write(16,512)wss(j)
 512      format(1x,'But this composition is already obtained ',A8)
         else
          do ii=1,kch
             sxtemp(ii)=sx(j,ii)
          enddo
          icss=icss+1
          isss(icss)=j
          if (icpdt.eq.0) goto 20
            iii=icpd
            do 22 ii=1,icpdt
            if (iii.eq.0) goto 24
            do 23 jj=1,iii
               if (ispd(jj).eq.ispdt(ii)) goto 22
  23        continue
  24        icpd=icpd+1
            ispd(icpd)=ispdt(ii)
  22        continue
  20      if (icpst.eq.0) goto 21
            iii=icps
            do 25 ii=1,icpst
            if (iii.eq.0) goto 26
            do 27 jj=1,iii
               if (isps(jj).eq.ispst(ii)) goto 25
  27        continue
  26        icps=icps+1
            isps(icps)=ispst(ii)
  25        continue
  21      continue
          call packghost()
          call psave(2)
         endif
        else
          print 103,wss(j)
          if (idebug.ge.2) write(8,103)wss(j)
          if (isummary.ge.2) write(15,103)wss(j)
          if (iwindow.ge.1) write(16,103)wss(j)
 103      format('   --- ',A8,' is not stable')
        endif
      endif
 12   continue
 10   continue
 11   return
      end

      subroutine temperature()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.for'
      DIMENSION DF(mlimit+islimit+1,mlimit+islimit+2)
      DIMENSION H(mlimit+islimit+1)
      dimension np(mlimit+islimit+1),ip(mlimit+islimit+1)
C     Calculate the equilibrium compositions at the temperature t0
C     with a given set of condensates.
C     This routine calculates logKs and calls subroutine grater() where
C     actual Newton-Raphson iterations take place, then searches for new pure
C     solids.  If new solids are found, add the one which has the highest
C     condensation temperature to the set of condensates, and go back to
C     the line where grater() is called. Repeat until no more new solids are
C     found. Note that this routine does not search for solid solutions.
C     Search for solid solutions is performed in the latter half of the
C     checktemp() subroutine.
      icpst=0
      icpdt=0
      icsst=0
      IC=0
      ID=0
ccc   Step in temperature
 201  T=T0*1.d-2
      T1=1d0/T
C     calc. logKs
      do I=1,NS+msExtra
        xk(i)=xc1(i)
        do k=2,iord
           xk(i)=xk(i)-xc(k,i)*(t1**(k-1))
        enddo
      enddo
      if ((imode.eq.0).or.(imode.eq.1)) then
         TG=(T/TI)**GAM
      endif
      call calcg0()
C
 51   MM=M+IS
      if((imode.eq.2).or.(imode.eq.3)) mm=m+is+1
      CALL GRATER(MM,mlimit+islimit+1,DF,H)
C
ccc  See whether to redo a step, continue, or advance T w/o checking.
      if (ifailconv.ne.0) return
      mmtemp=m+is
      if((imode.eq.2).or.(imode.eq.3)) mmtemp=m+is+1
      IF(mmtemp-MM+ICONT) 51,52,61
 52   L=1
      K=0
      if (isave2flag.eq.0) then
         call psave(3)
C        total pressure
         Ptotal=0d0
         do i=1,m+n
            Ptotal=Ptotal+P(i)
         enddo
         Ptotal=Ptotal*PTOT
C        total volume for 1 mol of elements
         Vtotal=rgas2*t0/tg/ptot
C        mole fractions
         call calcsx()
C        free energy of the system
         call freeenergy(G)
         GwoNew=G
      endif
      isave2flag=isave2flag+1
ccc  Check to see what phases have >eps abundance fraction, and label
C      print *,'n=',n,'is=',is,'nis=',nis
      DO 16 I=N1,ngp
C      print *,i
      if (is.eq.0) goto 100
      DO 70 J=N1,NIS
      IF(I-NQ(J)) 70,16,70
 70   CONTINUE
C     exclude polymorphs
          f2=0d0
          do k2=1,m
            if (dnu(i,k2).eq.0d0) goto 11
            if (f2.ne.0d0) goto 11
            f2=dnu(i,k2)
 11         continue
          enddo
      do ii=1,is
          f1=0d0
          do k2=1,m
            if (dnu(nq(n+ii),k2).eq.0d0) goto 10
            if (f1.ne.0d0) goto 10
            f1=dnu(nq(n+ii),k2)
 10         continue
          enddo
          iflag=0
          do k2=1,m
            if(dabs(dnu(nq(n+ii),k2)/f1-dnu(i,k2)/f2).gt.eps) iflag=1
          enddo
          if (iflag.eq.0) goto 16
      enddo
 100  continue
C     If Sum((stoichiom. coeff.)x(logP(element)))-logK > 0 (actually eps),
C     then label it as a candidate for a new condensate.
      P0=-XK(I)
      DO 17 J=1,M
 17   P0=P0+dnu(I,J)*DLP(J)
      XK1(I-N)=P0
      IF(P0-EPS) 16,18,18
 18   K=K+1
C      print *,'K=',k
      NP(K)=I
      IP(K)=0
 16   CONTINUE
ccc   for solid solutions
c**   disabled search for s.s.
c     np(k): species # (ns)
c     ip(k): solid solution # (nss)
 85   IF(K) 61,61,60
 61   continue
C     No new solids; calculate some variables and return.
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
      if (G-GwoNew.ne.0d0) then
         print 3031,GwoNew,G,G-GwoNew
         if (idebug.ge.2) write(8,3031) GwoNew,G,G-GwoNew
      endif
 3031 format(1X,'Total G without new condensates =',1pe18.11,
     >' kJ',/,1X,'Total G of the system (current) =',
     >1pe18.11,' kJ','  DG =',1pe11.4,' kJ')
      return
C
ccc   Check to see if any new phases should condense
 60   continue
      if ((imode.eq.0).or.(imode.eq.2)) then
         DELT=-T0
      elseif ((imode.eq.1).or.(imode.eq.3)) then
         DELT=-ptot
      endif
      KCH=1
      DO 62 I=1,K
      J=NP(I)-N
C     estimate the condensation temperatures
      if ((imode.eq.0).or.(imode.eq.2)) then
         DELT0=XK1(J)*DT/(XK1(J)-XK0(J))
      elseif ((imode.eq.1).or.(imode.eq.3)) then
         DELT0=XK1(J)*dp/(XK1(J)-XK0(J))
      endif
      if(idebug.ge.1) then
        if((imode.eq.0).or.(imode.eq.2)) then
         if(ip(i).eq.0) then
            print 108,np(i),wc(np(i)),I,DELT0+t0,XK1(J),XK0(J)
         else
            print 108,np(i),wss(ip(i)),I,DELT0+t0,XK1(J),XK0(J)
         endif
        elseif((imode.eq.1).or.(imode.eq.3)) then
         if(ip(i).eq.0) then
            print 109,np(i),wc(np(i)),I,-DELT0+ptot,XK1(J),XK0(J)
         else
            print 109,np(i),wss(ip(i)),I,-DELT0+ptot,XK1(J),XK0(J)
         endif
        endif
      endif
 108  FORMAT(1X,'New Condensate   :',I3,1x,A8,1x,I3,' T=',f10.4,
     >1x,1P2E11.3)
 109  FORMAT(1X,'New Condensate   :',I3,1x,A8,1x,I3,' P=',1pE10.3,
     >1x,1P2E11.3)
C     pick up the one which has the highest condensation temperature
      IF(DELT0-DELT) 62,62,64
 64   DELT=DELT0
      L=NP(I)
      LK=IP(I)
 62   CONTINUE
      IC=IC+1
      KPC(IC)=L
      if (lk.eq.0) then
        icpst=icpst+1
        ispst(icpst)=L
      else
        icsst=icsst+1
        issst(icsst)=lk
      endif
      XK1(L-N)=0d0
      if((imode.eq.0).or.(imode.eq.2)) then
         if(lk.eq.0) then
              WRITE(*,4006)wc(L),DELT0+t0
         else
              WRITE(*,4006)wss(lk),DELT0+t0
         endif
         if (idebug.ge.2) then
           if(lk.eq.0) then
              WRITE(8,4005)wc(L),DELT0+t0
           else
              WRITE(8,4005)wss(lk),DELT0+t0
           endif
         endif
      elseif((imode.eq.1).or.(imode.eq.3)) then
         if(lk.eq.0) then
              WRITE(*,4008)wc(L),-DELT0+ptot
         else
              WRITE(*,4008)wss(lk),-DELT0+ptot
         endif
         if (idebug.ge.2) then
           if(lk.eq.0) then
              WRITE(8,4007)wc(L),-DELT0+ptot
           else
              WRITE(8,4007)wss(lk),-DELT0+ptot
           endif
         endif
      endif
 4005 FORMAT('New Condensate   : ',A8,'    T= ',f10.4)
 4006 FORMAT(' New Condensate   : ',A8,'    T= ',f10.4)
 4007 FORMAT('New Condensate   : ',A8,'    P= ',1pE10.3)
 4008 FORMAT(' New Condensate   : ',A8,'    P= ',1pE10.3)
      FAC=1d0
      DO 65 I=1,M
      IF(NU(L,I)) 63,65,63
 63   FAC=DMIN1(R(I)/dnu(l,i),FAC)
 65   CONTINUE
      IF(LK.EQ.0) GO TO 81
      KCH=MS(LK)
      MSS(LK)=IS+m+1
 81   DO 82 I=1,KCH
      NIS=NIS+1
      if (nis-n.gt.islimit) then
         print *,' Exceed IS limit !!!'
         stop
      endif
C     set initial values for P
      if (lk.eq.0) then
        P(M+NIS)=factor*FAC
      elseif(nssspinel.ne.0) then
        if((m+nis-n.eq.mss(nssspinel)).or.
     >       (m+nis-n.eq.mss(nssspinel)+1)) then
           p(m+nis)=factor
        elseif(m+nis-n.eq.mss(nssspinel)+ms(nssspinel)-1) then
           p(m+nis)=factor*fac
        else
           P(M+NIS)=factor*FAC*sx(lk,i)
        endif
      else
           P(M+NIS)=factor*FAC*sx(lk,i)
      endif
C     set the species # (internal)
 82   NQ(NIS)=L+I-1
      IS=KCH+IS
      if (is.gt.islimit) then
         print *,' Exceed IS limit !!!'
         stop
      endif
C     repeat grater()
      GO TO 51
      end

      SUBROUTINE GRATER(MM,MMdim,DF,H)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.for'
      DIMENSION DF(MMdim,MMdim+1),H(MMdim)
      dimension np(nsslimit+1),H2(mlimit+islimit+1),
     >          Hslope(mlimit+islimit+1),
     >          Fslope(mlimit+islimit+1),Pold(mlimit+islimit+1)
C     Newton-Raphson iterations; also searches for disappearing phases.
C     This routine evaluates equations and derivatives and calls MatInv()
C     subroutine, then calculates new Ps.  This is repeated until the Ps
C     converge and until the residues of the equations become smaller
C     than some criterion.
C     After each iteration, the routine searches for condensed phases which
C     have very small Ps and are going to disappear. If it finds that any
C     phases should disappear, it removes those phases from the set of
C     condensates and returns to the temperature() routine.  The temperature()
C     routine then finds that some phases disappeared and calls this routine
C     again.
C
C     Newton-Raphson method:
C        to solve simultaneous equations: f(x)=0
C        x(new)=x(old)-inv(J)df
C          x : variables that must converge --- ln(P) in this routine
C          J : Jacobian matrix              --- DF=-J
C          df: residues of equations        --- H
C                                       then P(new)=P(old)*exp(inv(DF)H)
C
      inumerical=0
      ibacktrack=0
      hlimit=2d0
      eps2=dsqrt(eps)
      eps3=eps2*10
ccc   Begin Newton-Raphson iteration on simultaneous equations
      NITER=0
 50   NITER=NITER+1
C     evaluate equations
      call eqMassBalance(MM,MMdim,H,DF,inumerical)
      if (inumerical.ne.0) then
C     numerical calculation of Jacobian matrix for mass action law eqs.
      if (is.ne.0) then
      do i=m+1,mm
         H2(i)=H(i)
      enddo
      do j=1,mm
         k=j
         if (j.gt.m) k=j+n
         temp1=p(k)
         temp2=dlog(p(k))
         dx=eps2*dabs(temp2)
         if (dx.eq.0d0) dx=eps2
         x=temp2+dx
         p(k)=dexp(temp2+dx)
         dx=dlog(p(k))-temp2
         call eqMassAction(mm,MMdim,H,DF,inumerical)
         p(k)=temp1
         do i=m+1,mm
            df(i,j)=-(H(i)-H2(i))/dx
         enddo
      enddo
C     make sure that all values are the same as before calc. of Jacobian
      call eqMassAction(mm,MMdim,H,DF,inumerical)
      endif
      endif
C      if (idebug.ge.1) call freeenergy(G)
      if (idebug.ge.2) then
         do i=1,mm
            print 3031,i,(df(i,j),j=1,mm),h(i)
            write(8,3031) i,(df(i,j),j=1,mm),h(i)
         enddo
 3031    format (1X,I2,9(1pe8.1),10(/,3X,9(e8.1,:)))
      endif
C     check whether residues are small enough (< eps)
      do i=1,mm
         DF(i,mm+1)=H(i)
      enddo
      idflag=0
      idelta=0
      do i=1,m
         delta=dabs(h(i)/(r(i)*tg))
         if (delta.gt.eps) then
            idflag=1
            idelta=i
            difference=delta
         endif
      enddo
      if (is.ne.0) then
         do i=1,is
            delta=dabs(h(i+m)/xk(nq(i+n)))
            if (delta.gt.eps) then
               idflag=1
               idelta=nq(i+n)
               difference=delta
            endif
         enddo
      endif
      if((imode.eq.2).or.(imode.eq.3)) then
         delta=dabs(h(mm))
            if (delta.gt.eps) then
               idflag=1
               idelta=0
               difference=delta
            endif
      endif
      Hsum=0d0
      do i=1,mm
         Hsum=Hsum+H(i)*H(i)
         Fslope(i)=0d0
         do j=1,mm
            Fslope(i)=Fslope(i)-DF(j,i)*H(j)
         enddo
      enddo
      Hsum=0.5d0*Hsum
      CALL MATINV(MMdim,MMdim+1,DF,MM,MM+1)
      do i=1,mm
         Hslope(i)=DF(i,mm+1)
      enddo
C     Hslope is -dx=inv(DF)H(old)
      if (idebug.ge.2) then
         print 3032,(hslope(i),i=1,mm)
         write(8,3032) (hslope(i),i=1,mm)
 3032    format (1X,'d=',9(1pe8.1),10(/,3X,9(e8.1,:)))
      endif
C     the change of ln(P) (=dx) is restricted to less than hlimit (=2) since
C     abs(dx) may > (say) 10^5 which would make it impossible to calculate
C     new P (not ln(P)).
      HMAX=0d0
      do i=1,mm
        hmax=hmax+hslope(i)*hslope(i)
      enddo
      hmax=dsqrt(hmax)
      fnorm=1d0
      if (hmax.gt.hlimit) then
        fnorm=hlimit/hmax
      endif
      slope=0d0
      do i=1,mm
         slope=slope+Fslope(i)*Hslope(i)*fnorm
      enddo
      FCON=1d0
      do i=1,MM
        j=i
        if (i.gt.m) j=i+n
        if((i.eq.mm).and.((imode.eq.2).or.(imode.eq.3))) then
           TGold=TG
        else
           Pold(i)=P(j)
        endif
      enddo
 71   continue
      DO 40 I=1,MM
        j=i
        if (i.gt.m) j=i+n
C     new Ps
        if((i.eq.mm).and.((imode.eq.2).or.(imode.eq.3))) then
           TG=TGold*DEXP(Hslope(i)*fnorm*FCON)
        else
           P(J)=Pold(i)*DEXP(Hslope(i)*fnorm*FCON)
        endif
 40   continue
      if (ibacktrack.ne.0) then
C     Check whether Hsum decreases or not.
C     If not, backtrack along the Newton direction.
      call eqMassBalance(MM,MMdim,H,DF,inumerical)
      HsumNew=0d0
      do i=1,mm
         HsumNew=HsumNew+h(i)*h(i)
      enddo
      HsumNew=0.5d0*HsumNew
      if (idebug.ge.1) then
         print 3033,(h(i),i=1,mm)
         write(8,3033) (h(i),i=1,mm)
 3033    format (1X,'H=',9(1pe8.1),10(/,3X,9(e8.1,:)))
      endif
      ibflag=0
      if ((idflag.eq.0).or.(Hsum.lt.eps*eps)) then
C        Residues are small; need not line search
         if (idebug.ge.1) then
            print 75,fcon,HsumNew,Hsum,slope
            write(8,75) fcon,HsumNew,Hsum,slope
            print *,'Residues are small'
            write(8,*) 'Residues are small'
         endif
      elseif (fcon.lt.eps2) then
C        looks like all variables converged, but check to see it only at
C        local minimum of the function Hsum
         ibflag=1
c         if (idflag.ne.0) then
            if (idebug.ge.1) then
              print 3034,(Fslope(i),i=1,kch)
              write(8,3034) (Fslope(i),i=1,kch)
 3034         format (1X,'F=',9(1pe8.1,:),10(/,3X,9(e8.1,:)))
            endif
            ifalg=0
            do i=1,mm
               if (dabs(Fslope(i)).ge.eps2) iflag=1
            enddo
            if (iflag.eq.0) then
               print *,'Converged to a local minimum'
               if (idebug.ge.1) then
                  write(8,*) 'Converged to a local minimum'
               endif
               goto 300
            endif
C           something wrong --- consider roundoff error
            print *,'Problem in the backtrack routine'
            if (idebug.ge.1) then
               write(8,*) 'Problem in the backtrack routine'
            endif
            goto 300
c         endif
      elseif (HsumNew.le.Hsum+eps3*fcon*slope) then
C        sufficient decrease
         if (idebug.ge.1) then
            print 75,fcon,HsumNew,Hsum,slope
            write(8,75) fcon,HsumNew,Hsum,slope
            print *,'Sufficient decrease'
            write(8,*) 'Sufficient decrease'
         endif
      elseif ((hmax.gt.hlimit).and.(HsumNew.lt.Hsum)) then
C        Newton step is limited
         if (idebug.ge.1) then
            print 75,fcon,HsumNew,Hsum,slope
            write(8,75) fcon,HsumNew,Hsum,slope
            print *,'Newton step is limited'
            write(8,*) 'Newton step is limited'
         endif
      else
C        backtracking
         if (idebug.ge.1) then
            print 75,fcon,HsumNew,Hsum,slope
            write(8,75) fcon,HsumNew,Hsum,slope
            print *,'Backtracking'
            write(8,*) 'Backtracking'
 75         format(1X,'FCON=',1pE9.2,' HNew=',E9.2,
     >          ' HOld=',E9.2,' Slope=',E9.2)
         endif
         if (Fcon.eq.1d0) then
            fcontmp=-slope/(2d0*(HsumNew-Hsum-slope))
         else
            r1=HsumNew-Hsum-fcon*slope
            r2=Hsum2-Hsum-fcon2*slope
            a1=(r1/fcon/fcon-r2/fcon2/fcon2)/(fcon-fcon2)
            a2=(-fcon2*r1/fcon/fcon+fcon*r2/fcon2/fcon2)/(fcon-fcon2)
            if (a1.eq.0d0) then
               fcontmp=-slope/(2d0*a2)
            else
               fcontmp=(-a2+sqrt(a2*a2-3d0*a1*slope))/(3d0*a1)
            endif
         endif
         fcon2=fcon
         Hsum2=HsumNew
         fcontmp=min(fcontmp,fcon*0.5d0)
         fcon=max(fcontmp,fcon*0.1d0)
         goto 71
      endif
      endif
      if (is.eq.0) goto 500
ccc   Search for disappearing phases
      iis=is+1
      issflag=0
 600    iis=iis-1
        i=m+iis
        j=i+n
        K=n+iis
        FAC=1
        iss=0
        issp=0
C       see if this condensate is an endmember of a solid solution
        if (nss.ne.0) then
          do ii=1,nss
            if(mss(ii).ne.0) then
            if(i.eq.mss(ii)) then
                iss=ii
                issp=1
            elseif(i.eq.mss(ii)+ms(ii)-1) then
                iss=ii
                issp=3
            elseif((i.gt.mss(ii)).and.(i.lt.mss(ii)+ms(ii)-1)) then
                iss=ii
                issp=2
            endif
            endif
          enddo
        endif
C       FAC is an elemental abundance which determines the limit of
C       the abundance of species
        DO 165 ii=1,M
          IF(NU(nq(k),ii)) 163,165,163
 163      FAC=DMIN1(R(ii)/dnu(nq(k),ii),FAC)
 165    CONTINUE
C       criteria of phase disappearance:
C         (1) P/FAC < 1E-13
C         or
C         (2) P/FAC < disappear(1E-11), and dx (now Hslope) < -1E6
C       Usually only (2) occurs since Ps of the other species converge
C       before P/FAC of the disappearing phase becomes < 1E-13, causing
C       its dx < -1E6.
C       Although our selection of limits for these criteria were selected
C       so that we can be sure that a phase is really disappearing, this
C       routine might delete phases which should not disappear due to the
C       use of temporary data for the Ps (prior to convergence --- since
C       phases should disappear which haven't yet, convergence will not
C       occur and Ps are not final).
C       Actually, two phases sometimes disappear simultaneously,
C       but the temperature seeking routine handles this situation
C       by raising the temperature until only one phase disappears.
C
        iflag=0
        if(iss.eq.0) then
C          pure solid
           if((P(J)/FAC.LT.1.D-13).or.
     >     (Hslope(I).lt.-1.D6.and.P(J)/fac.lt.disappear)) iflag=1
        else
C          check all endmembers of the solid solution except for spinel
C          (check total P for spinel)
           if((P(J)/FAC.LT.1.D-13).or.
     >     (Hslope(I).lt.-1.D6.and.P(J)/fac.lt.disappear)) then
              if (issp.eq.1) then
                 if (issflag.eq.1) then
                    iflag=1
                 endif
              elseif (issp.eq.3) then
                 issflag=1
              endif
           else
              if (nssspinel.ne.0) then
                 if (iss.eq.nssspinel) then
                    if ((issp.eq.1).and.(issflag.eq.1)) then
                       iflag=1
                    endif
                 else
                    issflag=0
                 endif
              else
                 issflag=0
              endif
           endif
        endif
        imrgflag=0
        if ((iss.gt.nssorg).and.(issp.eq.1)) then
           pw=0d0
           do l=1,ms(iss)
              pw=pw+p(n+mss(iss)-1+l)
           enddo
           do l=1,ms(iss)
              sx(iss,l)=p(n+mss(iss)-1+l)/pw
           enddo
           do iss2=1,iss-1
              if (((nssghost(iss).eq.nssghost(iss2)).or.
     >          (nssghost(iss).eq.iss2)).and.(mss(iss2).gt.0)) then
                 pw=0d0
                 do l=1,ms(iss2)
                    pw=pw+p(n+mss(iss2)-1+l)
                 enddo
                 do l=1,ms(iss2)
                    sx(iss2,l)=p(n+mss(iss2)-1+l)/pw
                 enddo
                 imrgflag=1
                 do l=1,ms(iss)
                 if(dabs(sx(iss,l)-sx(iss2,l)).gt.eps2/10d0)imrgflag=0
                 enddo
                 if (imrgflag.eq.1) then
                    iss2t=iss2
                    goto 603
                 endif
              endif
           enddo
        endif
        if (iflag.eq.0) goto 602
 603    continue
 145    issflag=0
C       print the disappearing phase
        if (iss.ne.0) then
           pw=0d0
           if (iss.eq.nssspinel) then
              pw=p(n+mss(iss)+ms(iss)-1)
           else
              do l=1,ms(iss)
                 pw=pw+p(n+mss(iss)-1+l)
              enddo
           endif
          if ((imrgflag.eq.1).and.(iflag.eq.0)) then
           do l=1,ms(iss)
              p(n+mss(iss2t)+l-1)=p(n+mss(iss2t)+l-1)
     >                         +p(n+mss(iss)+l-1)
           enddo
           if (idebug.eq.0) print 4009,wss(iss),pw
           if (idebug.ge.1) print 3,nq(k),wss(iss),I-m,pw,Hslope(I)
           if (idebug.ge.2) write(8,3)nq(k),wss(iss),I-m,pw,Hslope(I)
          else
           if (idebug.eq.0) print 4008,wss(iss),pw
           if (idebug.ge.1) print 2,nq(k),wss(iss),I-m,pw,Hslope(I)
           if (idebug.ge.2) write(8,2)nq(k),wss(iss),I-m,pw,Hslope(I)
          endif
        else
          if (idebug.eq.0) print 4008,wc(nq(k)),P(J)
          if (idebug.ge.1) print 2,nq(k),wc(nq(k)),I-m,P(J),Hslope(I)
          if (idebug.ge.2) write(8,2)nq(k),wc(nq(k)),I-m,P(J),Hslope(I)
        endif
 2      FORMAT(' PHASE DISAPPEARS :',I3,' ',a8,' FUNC.=',I3,
     1  ' PRESS.=',1pe10.3,' CHANGE=',1PE10.3)
 4008   FORMAT(' PHASE DISAPPEARS : ',a8,' P/Ptotal=',1pe10.3)
 3      FORMAT(' PHASE MERGED :',I3,' ',a8,' FUNC.=',I3,
     1  ' PRESS.=',1pe10.3,' CHANGE=',1PE10.3)
 4009   FORMAT(' PHASE MERGED : ',a8,' P/Ptotal=',1pe10.3)
        ID=ID+1
        KPD(ID)=NQ(K)
        icpdt=icpdt+1
        if (iss.eq.0) then
           ispdt(icpdt)=nq(k)
        else
           ispdt(icpdt)=nssp(iss)+ms(iss)-1
        endif
C       delete the disappearing phase
        KCH=1
        if(nss.ne.0) then
           DO 61 L=1,NSS
              NP(L)=0
              IF(MSS(L).EQ.0) GO TO 61
              IF(I.GT.MSS(L)+MS(L)-1) GO TO 61
              IF(I.LT.MSS(L)) GO TO 62
              K=MSS(L)+n-m
              KCH=MS(L)
              MSS(L)=KCH
 62           NP(L)=1
 61        CONTINUE
c*         print *,(nq(ii),ii=n1,nis)
           DO 64 L=1,NSS
 64        MSS(L)=MSS(L)-KCH*NP(L)
c          print *,(mss(l),l=1,nss)
        endif
        IS=IS-KCH
        NIS=N+IS
        IF(K.GT.NIS) goto 602
        DO 43 II=K,NIS
            II1=II+KCH
            P(II+M)=P(II1+M)
            NQ(II)=NQ(II1)
 43     continue
 602  if (iis.gt.1) goto 600
C     if any phases disappeared, return
      mmtemp=m+is
      if ((imode.eq.2).or.(imode.eq.3)) mmtemp=m+is+1
      if (mmtemp.lt.mm) return
 500  continue
C     end of the search for disappearing phases
      DO 42 K=1,M
 42   DLP(K)=DLOG10(P(K))
      do J=1,N
         P0=-XK(J)
         do II=1,M
            P0=P0+DLP(II)*dnu(J,II)
         enddo
         P(M+J)=10d0**P0
      enddo
C     check for convergence
 44   DO 20 II=1,MM
      IF(DABS(Hslope(II))-EPS2) 20,20,22
 20   CONTINUE
      if((idflag.ne.0).and.(niter.gt.10)
     >.and.(idebug.ge.1)) then
        if (idelta.eq.0) then
           write(*,'(1x,A2,2x,''Diff.='',e10.4)')'TG',difference
           if (idebug.ge.2)
     >     write(8,'(1x,A2,2x,''Diff.='',e10.4)')'TG',difference
        elseif (idelta.le.m) then
           write(*,'(1x,A2,2x,''Diff.='',e10.4)')we(idelta),difference
           if (idebug.ge.2)
     >     write(8,'(1x,A2,2x,''Diff.='',e10.4)')we(idelta),difference
        else
           write(*,'(1x,A8,2x,''Diff.='',e10.4)')wc(idelta),difference
           if (idebug.ge.2)
     >     write(8,'(1x,A8,2x,''Diff.='',e10.4)')wc(idelta),difference
        endif
      endif
C     also residues must be small enough
      if(idflag.ne.0) goto 22
      RETURN
  22  continue
C     allow 100 iterations
      IF(NITER-100) 50,50,300
C
cccc  Need for operator intervention.
 300  if (idebug.ge.1) WRITE(8, 105) NITER
 105  FORMAT(26H FAILED TO CONVERGE AFTER ,I3,12H ITERATIONS.)
      PRINT 105,NITER
      if (idflag.ne.0) then
        if (idelta.eq.0) then
           write(*,'(1x,A2,2x,''Diff.='',e10.4)')'TG',difference
           if (idebug.ge.2)
     >     write(8,'(1x,A2,2x,''Diff.='',e10.4)')'TG',difference
        elseif (idelta.le.m) then
           write(*,'(1x,A2,2x,''Diff.='',e10.4)')we(idelta),difference
           if (idebug.ge.2)
     >     write(8,'(1x,A2,2x,''Diff.='',e10.4)')we(idelta),difference
        else
           write(*,'(1x,A8,2x,''Diff.='',e10.4)')wc(idelta),difference
           if (idebug.ge.2)
     >     write(8,'(1x,A8,2x,''Diff.='',e10.4)')wc(idelta),difference
        endif
      endif
 99   print 3020,(nq(k+n),wc(nq(k+n)),p(k+m+n),k=1,is)
      if ((imode.eq.2).or.(imode.eq.3))
     > print 3020,is+1,'TG',tg
 3020 format (3(1x,I3,1x,A8,1pe10.3,3x))
      PRINT 98,'Condensate: ',(KPC(I),I=1,IC)
      PRINT 98,'Disappear:  ',(KPD(I),I=1,ID)
      PRINT 1,T0,ptot
  1   format(' T = ',F9.3, 'K,   P =',1pe10.3)
 98   FORMAT(1X,a12,30I4)
      iflag=0
      if ((idebug.ge.1).or.(iauto.eq.0)) then
         write (*,'(A)') ' Try higher T (or lower P)? (yes=1/No=0) '
         read (*,*) iflag
      else
         iflag=1
      endif
      if (iflag.ne.0) then
         ifailconv=1
         WRITE(*,'(A)')' Try higher T automatically (Yes=1 / No=0) ? '
         READ (*,*) iauto
         return
      endif
      print *,' Tidy up phases (i,n,nss,P), then input (0,0,0,0) ',
     >'and proper IS.'
      print *,' To stop here, input -1 as IS.'
      print *,' To continue as if it converged, input (-1,-1,0,0).'
      print *,' To erase solid solution, input (0,0,-nss,0).'
      print *,'  i: function #    n: species #    nss: solid sol. #'
      print *,'  P: partial pressure;   IS: total no. of condensates'
 92   READ *,I,J,K,PX
ccc   To continue to next T, tidy up phases, then enter -1,0,proper IS,0
      IF(I) 93,91,94
 93   continue
      IF(J.LT.0) goto 96
      ICONT=1
      IS=K
      GO TO 96
 94   IF(K) 191,192,193
 191   MSS(-K)=0
      GO TO 92
 193  MSS(K)=M+I
 192  continue
      if (i.gt.is) then
         tg=px
      else
         NQ(I+N)=J
         P(m+n+I)=PX
         WRITE(8,95) I,J
      endif
 95   FORMAT(' FUNCTION ',I2,' REPLACED BY SPECIE ',I3)
      GO TO 92
 91   READ *,IS
      IF(IS.EQ.-1) then
         ifailconv=1
         istopflag=1
         return
      endif
      MM=100
 96   NIS=N+IS
      WRITE(*,'(A)') ' Debug mode (Level 1-2 / No=0) ? '
      READ (*,*) idebug
      RETURN
      END

      SUBROUTINE eqMassBalance(MM,MMdim,H,DF,inumerical)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(MMdim),DF(MMdim,MMdim+1)
      INCLUDE 'con0621i.for'
ccc   Evaluate equations and their derivatives
      do I=1,M
        DLP(I)=DLOG10(P(I))
      enddo
C     new Ps for gaseous species
      DO 15 J=1,N
        P0=-XK(J)
        DO 6 II=1,M
 6      P0=P0+DLP(II)*dnu(J,II)
 15   P(M+J)=10d0**P0
      if (idebug.ge.1) then
        print *,niter
        if (idebug.ge.2) write(8,*) niter
        if (idebug.ge.2) print 3011,(k,wee(k),p(k),k=1,m)
        if (idebug.ge.2) write(8,3011) (k,wee(k),p(k),k=1,m)
        print 3011,(nq(k+n),wc(nq(k+n)),p(k+m+n),k=1,is)
        if (idebug.ge.2)
     >     write(8,3011) (nq(k+n),wc(nq(k+n)),p(k+m+n),k=1,is)
 3011   format (3(1x,I3,1x,A8,1pe10.3,3x))
      endif
C     Mass balance equations for elements
C     (Abund.)-P(element)-Sum((stoichiom.coeff.)xP)(gas&condensed phases) = 0
      DO 25 I=1,M
        F0=TG*R(I)-P(I)
        DO 3008 J=1,NIS
          if(nssspinel.eq.0) goto 26
          if(mss(nssspinel).eq.0) goto 26
          if((j-n.ge.mss(nssspinel)-m).and.
     >      (j-n.le.mss(nssspinel)+ms(nssspinel)-1-m)) goto 3008
 26       F0=F0-P(J+M)*dnu(NQ(J),I)
 3008   continue
        DO 27 J=1,M
          P0=0d0
          DO 28 II=1,N
 28       P0=P0+P(M+II)*dnu(II,J)*dnu(II,I)
          IF(I.EQ.J) P0=P0+P(J)
 27     DF(I,J)=P0
 25   H(I)=F0
      mn=n-m
      IF(IS.EQ.0) GO TO 51
C     equilibrium equations (mass action law) for condensates
C     logK - Sum((stoichiom.coeff.) x logP(element)) = 0
      DO 32 I=N1,NIS
        F0=XK(NQ(I))
        DO 33 J=1,M
          DF(I-MN,J)=dnu(NQ(I),J)/dlog(10d0)
          DF(J,I-MN)=P(I+M)*dnu(NQ(I),J)
 33     F0=F0-DLP(J)*dnu(NQ(I),J)
        DO 37 J=m+1,MM
 37     DF(I-MN,J)=0d0
 32   H(I-MN)=F0
      IF(NSS.EQ.0) GO TO 51
C     for solid solutions
C     logK - Sum((stoichiom.coeff.) x logP(element)) + (logX + logç) = 0
      DO 52 ISS=1,nss
        IF(MSS(ISS).EQ.0) GO TO 52
        K=MSS(ISS)
        IMAX=MS(ISS)+K-1
        PT=0d0
      if (iss.eq.nssspinel) then
cccc    For spinel solid solution
C       spinel solid solution is treated separately because Ps are different
C       P(#1): X(Cr)=Cr/(Cr+Al) molar ratio
C       P(#2): X(Fe)=Fe/(Fe+Mg) molar ratio
C       P(#3): total P of spinel s.s.
        if (p(n+k+1).ge.1d0) then
         p(n+k+1)=1.-eps
         print *,'  Warning: X(Fe) in Spinel >= 1'
         if (idebug.ge.2) write(8,*) '  Warning: X(Fe) in Spinel >= 1'
        endif
        xfe=p(n+k+1)
        xmg=1d0-xfe
        if (p(n+k).ge.1d0) then
         p(n+k)=1.-eps
         print *,'  Warning: X(Cr) in Spinel >= 1'
         if (idebug.ge.2) write(8,*) '  Warning: X(Cr) in Spinel >= 1'
        endif
        xcr=p(n+k)
        xal=1d0-xcr
        pt=p(n+imax)
c        print *,xal,xcr,xmg,xfe,pt
        a=1d0/dlog(1.d1)
 107    format(1x,24i3)
        DXMG=DLOG10(XMG)
        DXFE=DLOG10(XFE)
        DXAL=DLOG10(XAL)
        dxcr=dlog10(xcr)
        do 3012 i=1,m
        do 3012 j=k,imax
        df(i,j)=0d0
 3012   continue
        h(mal)=h(mal)-2d0*xal*pt
        df(mal,imax)=2d0*xal*pt
        df(mal,k)=-2d0*xcr*pt
        h(mcr)=h(mcr)-2d0*xcr*pt
        df(mcr,imax)=2d0*xcr*pt
        df(mcr,k)=2d0*xcr*pt
        h(mmg)=h(mmg)-xmg*pt
        df(mmg,imax)=xmg*pt
        df(mmg,k+1)=-xfe*pt
        h(mfe)=h(mfe)-xfe*pt
        df(mfe,imax)=xfe*pt
        df(mfe,k+1)=xfe*pt
        h(mO)=h(mO)-4d0*pt
        df(mO,imax)=4d0*pt
        do 3009 i=k,imax
        imn=i-m+n
        h(i)=h(i)+dnu(nq(imn),mal)*dxal+dnu(nq(imn),mcr)*dxcr
     >  +dnu(nq(imn),mmg)*dxmg+dnu(nq(imn),mfe)*dxfe
        df(i,k)=-a*(dnu(nq(imn),mcr)-dnu(nq(imn),mal)*xcr/xal)
        df(i,k+1)=-a*(dnu(nq(imn),mfe)-dnu(nq(imn),mmg)*xfe/xmg)
        df(i,imax)=0d0
 3009   continue
      else
C       ideal solid solutions
C       add +logX (mole fraction) to mass action law equations
        DO 53 I=K,IMAX
 53     PT=PT+P(N+I)
        A=dble(MDB(ISS))/dlog(10d0)
        DO 54 I=K,IMAX
        DO 55 J=K,IMAX
 55     DF(I,J)=A*P(N+J)/PT
        sx(iss,i-k+1)=P(N+I)/PT
        H(I)=H(I)+DLOG10(sx(iss,i-k+1))*dble(MDB(ISS))
 54     DF(I,I)=DF(I,I)-A
      endif
C       non-ideality is treated in a separate routine
C       add +logç (activity coefficient) to mass action law equations
        idf=1
        if (inumerical.ne.0) idf=0
        if (nonid(iss).ne.0)
     >     call hnonideal(MMdim,H,DF,iss,inumerical,idf)
 52   continue
 51   continue
      if((imode.eq.2).or.(imode.eq.3)) then
C        Pressure regulation
         H(mm)=1d0
         do i=1,m+n
            H(mm)=H(mm)-P(i)
         enddo
         do i=1,m
            DF(i,mm)=-R(i)*TG
         enddo
         if (is.ne.0) then
            do i=m+1,m+is
               DF(i,mm)=0d0
               DF(mm,i)=0d0
            enddo
         endif
         DF(mm,mm)=0d0
         do i=1,m
            F0=P(i)
            do j=1,n
               F0=F0+P(J+M)*dnu(NQ(J),I)
            enddo
            DF(mm,i)=F0
         enddo
      endif
      return
      end

      SUBROUTINE eqMassAction(mm,MMdim,H,DF,inumerical)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(MMdim),DF(MMdim,MMdim+1)
      INCLUDE 'con0621i.for'
      do I=1,M
        DLP(I)=DLOG10(P(I))
      enddo
      mn=n-m
      IF(IS.EQ.0) GO TO 51
C     equilibrium equations (mass action law) for condensates
C     logK - Sum((stoichiom.coeff.) x logP(element)) = 0
      DO 32 I=N1,NIS
        F0=XK(NQ(I))
        DO 33 J=1,M
C          DF(I-MN,J)=dnu(NQ(I),J)/dlog(10d0)
C          DF(J,I-MN)=P(I+M)*dnu(NQ(I),J)
 33     F0=F0-DLP(J)*dnu(NQ(I),J)
C        DO 37 J=m+1,MM
C 37     DF(I-MN,J)=0d0
 32   H(I-MN)=F0
      IF(NSS.EQ.0) GO TO 51
C     for solid solutions
C     logK - Sum((stoichiom.coeff.) x logP(element)) + (logX + logç) = 0
      DO 52 ISS=1,nss
        IF(MSS(ISS).EQ.0) GO TO 52
        K=MSS(ISS)
        IMAX=MS(ISS)+K-1
        PT=0d0
      if (iss.eq.nssspinel) then
cccc    For spinel solid solution
C       spinel solid solution is treated separately because Ps are different
C       P(#1): X(Cr)=Cr/(Cr+Al) molar ratio
C       P(#2): X(Fe)=Fe/(Fe+Mg) molar ratio
C       P(#3): total P of spinel s.s.
        if (p(n+k+1).ge.1d0) then
         p(n+k+1)=1.-eps
         print *,'  Warning: X(Fe) in Spinel >= 1'
         if (idebug.ge.2) write(8,*) '  Warning: X(Fe) in Spinel >= 1'
        endif
        xfe=p(n+k+1)
        xmg=1d0-xfe
        if (p(n+k).ge.1d0) then
         p(n+k)=1.-eps
         print *,'  Warning: X(Cr) in Spinel >= 1'
         if (idebug.ge.2) write(8,*) '  Warning: X(Cr) in Spinel >= 1'
        endif
        xcr=p(n+k)
        xal=1d0-xcr
        pt=p(n+imax)
c        print *,xal,xcr,xmg,xfe,pt
        a=1d0/dlog(1.d1)
 107    format(1x,24i3)
        DXMG=DLOG10(XMG)
        DXFE=DLOG10(XFE)
        DXAL=DLOG10(XAL)
        dxcr=dlog10(xcr)
C        do 3012 i=1,m
C        do 3012 j=k,imax
C        df(i,j)=0d0
C 3012   continue
C        h(mal)=h(mal)-2d0*xal*pt
C        df(mal,imax)=2d0*xal*pt
C        df(mal,k)=-2d0*xcr*pt
C        h(mcr)=h(mcr)-2d0*xcr*pt
C        df(mcr,imax)=2d0*xcr*pt
C        df(mcr,k)=2d0*xcr*pt
C        h(mmg)=h(mmg)-xmg*pt
C        df(mmg,imax)=xmg*pt
C        df(mmg,k+1)=-xfe*pt
C        h(mfe)=h(mfe)-xfe*pt
C        df(mfe,imax)=xfe*pt
C        df(mfe,k+1)=xfe*pt
C        h(mO)=h(mO)-4d0*pt
C        df(mO,imax)=4d0*pt
        do 3009 i=k,imax
        imn=i-m+n
        h(i)=h(i)+dnu(nq(imn),mal)*dxal+dnu(nq(imn),mcr)*dxcr
     >  +dnu(nq(imn),mmg)*dxmg+dnu(nq(imn),mfe)*dxfe
C        df(i,k)=-a*(dnu(nq(imn),mcr)-dnu(nq(imn),mal)*xcr/xal)
C        df(i,k+1)=-a*(dnu(nq(imn),mfe)-dnu(nq(imn),mmg)*xfe/xmg)
C        df(i,imax)=0d0
 3009   continue
      else
C       ideal solid solutions
C       add +logX (mole fraction) to mass action law equations
        DO 53 I=K,IMAX
 53     PT=PT+P(N+I)
C        A=dble(MDB(ISS))/dlog(10d0)
        DO 54 I=K,IMAX
C        DO 55 J=K,IMAX
C 55     DF(I,J)=A*P(N+J)/PT
        sx(iss,i-k+1)=P(N+I)/PT
        H(I)=H(I)+DLOG10(sx(iss,i-k+1))*dble(MDB(ISS))
C 54     DF(I,I)=DF(I,I)-A
  54    continue
      endif
C       non-ideality is treated in a separate routine
C       add +logç (activity coefficient) to mass action law equations
        idf=0
        if (nonid(iss).ne.0)
     >     call hnonideal(MMdim,H,DF,iss,inumerical,idf)
 52   continue
 51   continue
      if((imode.eq.2).or.(imode.eq.3)) then
C        Pressure regulation
         H(mm)=1d0
         do i=1,m+n
            H(mm)=H(mm)-P(i)
         enddo
      endif
      return
      end
