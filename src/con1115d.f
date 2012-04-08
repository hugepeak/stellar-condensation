C     "con1115d.for"   11/16/94
      subroutine tables()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      dimension pw(nslimit+26),itbwflag(nslimit),p2(nslimit)
      CHARACTER*8  ww(nslimit+26),wex*4
      CHARACTER*12 filename
ccc   make table-like output files (TAB-delimited text files),
c     which are readable by spreadsheet programs such as MS-Excel
c     itable=1 : master table only
c     itable=2 : master + specified elements
c     itable=3 : master + specified elements + solid solutions
c     itable=4 : master + specified elements + solid solutions + wt%
c     itable=5 : master + all elements + solid solutions + wt%
      if (ntable.eq.0) return
      rgas=rgas2
c      print *,itable
c     master table
      do j=1,ns+msExtra
         itbwflag(j)=0
         do i=1,m
            if (itbflag(i,j).eq.1) itbwflag(j)=1
         enddo
      enddo
      filename=contfile
      i=index(filename,'.')
      j=index(filename,' ')
      if (i.eq.0) then
         if ((j.eq.0).or.(j.gt.9)) then
            filename(9:12)='.TAB'
         else
            filename(j:j+3)='.TAB'
         endif
      else
         filename(i:i+3)='.TAB'
      endif
      OPEN (8,FILE=filename)
      print *,'Making : ',filename
      write(8,'(A10,A50)') 'Program : ',version
      write(8,'(A15,A12)') 'Control File : ',contfile
      write(8,'(A79)') comcntl
      write(8,'(A26,A12)') 'Thermodynamic Data File : ',datafile
      write(8,'(A79)') comdata
      if ((nssspinel.ne.0).and.(nonid(nssspinel).eq.0)) then
      write(8,3002) wc(nssp(nssspinel)+ms(nssspinel))
 3002 format(A8,'was ignored for the calculation of the spinel solid',
     >' solution ')
      endif
      if ((imode.eq.0).or.(imode.eq.2)) then
         write(8,300) dtmin
 300  format('Temperature resolution was set to ',F5.2,' K')
      elseif ((imode.eq.1).or.(imode.eq.3)) then
        if (ibar.eq.0) then
         write(8,302) dpmin
        else
         write(8,1302) dpmin
        endif
 302  format('Pressure resolution was set to ',F5.2,' log10(atm)')
1302  format('Pressure resolution was set to ',F5.2,' log10(bar)')
      endif
      if (ibar.eq.0) then
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(8,101)PTOT
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(8,99)PTOT
      endif
      else
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(8,1101)PTOT
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(8,199)PTOT
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
      WRITE(8,48)(WE(I),char(9),R(I),I=1,M)
 48   format(a2,a,1pe12.5)
      write(8,'(/,''All data are mol/l unless otherwise noted'',/)')
      ww(1)='Condns.'
      ww(2)='Disapp.'
      ww(3)='T (K)'
      if (ibar.eq.0) then
         ww(4)='P (atm)'
      else
         ww(4)='P (bar)'
      endif
      ww(5)='RT'
      ww(6)='LogPO2'
      ww(7)='# Solids'
      iw=7
      do i=1,m
         ww(i+iw)=wee(i)
       enddo
      iw=iw+m
      do i=1,ns+msExtra
         if (itbwflag(i).eq.1) then
            iw=iw+1
            ww(iw)=wc(i)
C            print *, ww(iw)
         endif
      enddo
C       pause
      write(8,10) (ww(i),char(9),i=1,iw)
 10   format(500(a8,a))
      do it=1,ntable
         t0=ttable(it)
         T=T0*1.d-2
         T1=1d0/T
          ptot=ptottb(it)
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
         do I=1,N
            xk(i)=xc1(i)
            do k=2,iord
               xk(i)=xk(i)-xc(k,i)*(t1**(k-1))
            enddo
         enddo
         do i=1,M
            DLP(i)=DLOG10(pmtable(i,it))
         enddo
         do i=1,N
            P0=-XK(i)
            do j=1,M
               P0=P0+DLP(j)*dnu(i,j)
            enddo
            p2(i)=10d0**P0
         enddo
         do i=1,n
            p2(i)=p2(i)*ptot/(rgas*t0)
         enddo
         do i=1,ist
           if (nssspinel.eq.0) then
                p2(nqt(i))=pstable(i,it)*ptot/(rgas*t0)
           else
             isp=nssp(nssspinel)
             if ((nqt(i).eq.isp).or.(nqt(i).eq.isp+1)) then
                p2(nqt(i))=pstable(i,it)
             else
                p2(nqt(i))=pstable(i,it)*ptot/(rgas*t0)
             endif
           endif
         enddo
         ww(1)=wcond(it)
         ww(2)=wdisap(it)
         pw(1)=t0
         pw(2)=ptot 
c        pttb(it)
         pw(3)=rgas*t0
c         rgas*t0
         pw(4)=pO2tb(it)
         pw(5)=dble(istb(it))
C         print *,time
         iw=5
         do i=1,m
            pw(i+iw)=pmtable(i,it)*ptot/(rgas*t0)
         enddo
         iw=iw+m
         do i=1,ns+msExtra
            if (itbwflag(i).eq.1) then
               iw=iw+1
               pw(iw)=p2(i)
            endif
         enddo
         write(8,11) (ww(i),char(9),i=1,2),(pw(i),char(9),i=1,iw)
 11      format(2(a8,a),500(1pE12.5,a))
      enddo
      close(8)
c     table for each element
      if (itable.lt.2) goto 1000
      if (itable.ge.5) then
         do i=1,m
            wtb(i)=we(i)
         enddo
      endif
      do ie=1,m
         iflag=0
         do i=1,m
            if (wtb(i).eq.we(ie)) iflag=1
         enddo
      if (iflag.eq.1) then
         filename=contfile
         i=index(filename,'.')
         j=index(filename,' ')
         wex='.   '
         wex(2:3)=we(ie)
         if (i.eq.0) then
            if ((j.eq.0).or.(j.gt.9)) then
               filename(9:12)=wex
            else
               filename(j:j+3)=wex
            endif
         else
            filename(i:i+3)=wex
         endif
         OPEN (8,FILE=filename)
         print *,'Making : ',filename
         write(8,'(A10,A50)') 'Program : ',version
         write(8,'(A15,A12)') 'Control File : ',contfile
         write(8,'(A79)') comcntl
         write(8,'(A26,A12)') 'Thermodynamic Data File : ',datafile
         write(8,'(A79)') comdata
         if ((nssspinel.ne.0).and.(nonid(nssspinel).eq.0)) then
            write(8,3002) wc(nssp(nssspinel)+ms(nssspinel))
         endif
      if ((imode.eq.0).or.(imode.eq.2)) then
         write(8,300) dtmin
c 300  format('Temperature resolution was set to ',F5.2,' K')
      elseif ((imode.eq.1).or.(imode.eq.3)) then
        if (ibar.eq.0) then
         write(8,302) dpmin
        else
         write(8,1302) dpmin
        endif
c 302  format('Pressure resolution was set to ',F5.2,' log10(atm)')
      endif
      if (ibar.eq.0) then
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(8,101)PTOT
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(8,99)PTOT
      endif
      else
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(8,1101)PTOT
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(8,199)PTOT
      endif
      endif
c 101  FORMAT(/,'P initial =',1PE10.3,' atm ',
c     >'(if H2 is dominant, the actual pressure is half)',
c     >/,'Atomic Abundances relative to their sum = 1 :')
c 99   FORMAT(/,'P initial =',1PE10.3,' atm ',
c     >/,'Atomic Abundances relative to their sum = 1 :')
         WRITE(8,48)(WE(I),char(9),R(I),I=1,M)
         write(8,'(/,''Element : '',A,A2,/,
     >   ''All data are fractions of the element '',
     >   ''unless otherwise noted'',/)') char(9),we(ie)
         ww(1)='Condns.'
         ww(2)='Disapp.'
         ww(3)='T (K)'
         if (ibar.eq.0) then
            ww(4)='P (atm)'
         else
            ww(4)='P (bar)'
         endif
         ww(5)='RT'
         ww(6)='N(mol/l)'
         if (ibar.eq.0) then
            ww(7)='P(  )atm'
         else
            ww(7)='P(  )bar'
         endif
         ww(7)(3:4)=we(ie)
         iw=8
         ww(iw)=wee(ie)
         do i=1,ns+msExtra
            if (itbflag(ie,i).eq.1) then
               iw=iw+1
               ww(iw)=wc(i)
            endif
         enddo
         write(8,10) (ww(i),char(9),i=1,iw)
         do it=1,ntable
            t0=ttable(it)
            T=T0*1.d-2
            T1=1d0/T
          ptot=ptottb(it)
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
            do I=1,N
               xk(i)=xc1(i)
               do k=2,iord
                  xk(i)=xk(i)-xc(k,i)*(t1**(k-1))
               enddo
            enddo
            do i=1,M
               DLP(i)=DLOG10(pmtable(i,it))
            enddo
            do i=1,N
               P0=-XK(i)
               do j=1,M
                  P0=P0+DLP(j)*dnu(i,j)
               enddo
               p2(i)=10d0**P0
            enddo
            TG=tgtb(it)
            RP=r(ie)*TG
            do i=1,n
               p2(i)=p2(i)*dnu(nq(i),ie)/rp
            enddo
            do i=1,ist
             if (nssspinel.eq.0) then
                  p2(nqt(i))=pstable(i,it)*dnu(nqt(i),ie)/rp
             else
               isp=nssp(nssspinel)
               if ((nqt(i).eq.isp).or.(nqt(i).eq.isp+1)) then
                  p2(nqt(i))=pstable(i,it)
               elseif(nqt(i).eq.isp+2) then
                  rnt=0d0
                  pt=pstable(i,it)
                  xcr=p2(isp)
                  xal=1d0-xcr
                  xfe=p2(isp+1)
                  xmg=1d0-xfe
                  if(ie.eq.mal) rnt=2d0*xal*pt/rp
                  if(ie.eq.mcr) rnt=2d0*xcr*pt/rp
                  if(ie.eq.mmg) rnt=xmg*pt/rp
                  if(ie.eq.mfe) rnt=xfe*pt/rp
                  if(ie.eq.mO)  rnt=4d0*pt/rp
                  p2(nqt(i))=rnt
               else
                  p2(nqt(i))=pstable(i,it)*dnu(nqt(i),ie)/rp
               endif
             endif
            enddo
            ww(1)=wcond(it)
            ww(2)=wdisap(it)
            pw(1)=t0
            pw(2)=pttb(it)
            pw(3)=rgas*t0
            pw(4)=rp*ptot/(rgas*t0)
            pw(5)=pmtable(ie,it)*ptot
            iw=6
            pw(iw)=pmtable(ie,it)/rp
             do i=1,ns+msExtra
               if (itbflag(ie,i).eq.1) then
                  iw=iw+1
                  pw(iw)=p2(i)
               endif
            enddo
            write(8,11) (ww(i),char(9),i=1,2),(pw(i),char(9),i=1,iw)
         enddo
         close(8)
      endif
      enddo
c     table for solid solutions
      if (itable.lt.3) goto 1000
      if (nss.eq.0) goto 1000
         filename=contfile
         i=index(filename,'.')
         j=index(filename,' ')
         wex='.SS '
         if (i.eq.0) then
            if ((j.eq.0).or.(j.gt.9)) then
               filename(9:12)=wex
            else
               filename(j:j+3)=wex
            endif
         else
            filename(i:i+3)=wex
         endif
         OPEN (8,FILE=filename)
         print *,'Making : ',filename
         write(8,'(A10,A50)') 'Program : ',version
         write(8,'(A15,A12)') 'Control File : ',contfile
         write(8,'(A79)') comcntl
         write(8,'(A26,A12)') 'Thermodynamic Data File : ',datafile
         write(8,'(A79)') comdata
         if ((nssspinel.ne.0).and.(nonid(nssspinel).eq.0)) then
            write(8,3002) wc(nssp(nssspinel)+ms(nssspinel))
         endif
      if ((imode.eq.0).or.(imode.eq.2)) then
         write(8,300) dtmin
c 300  format('Temperature resolution was set to ',F5.2,' K')
      elseif ((imode.eq.1).or.(imode.eq.3)) then
        if (ibar.eq.0) then
         write(8,302) dpmin
        else
         write(8,1302) dpmin
        endif
c 302  format('Pressure resolution was set to ',F5.2,' log10(atm)')
      endif
      if (ibar.eq.0) then
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(8,101)PTOT
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(8,99)PTOT
      endif
      else
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(8,1101)PTOT
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(8,199)PTOT
      endif
      endif
c 101  FORMAT(/,'P initial =',1PE10.3,' atm ',
c     >'(if H2 is dominant, the actual pressure is half)',
c     >/,'Atomic Abundances relative to their sum = 1 :')
c 99   FORMAT(/,'P initial =',1PE10.3,' atm ',
c     >/,'Atomic Abundances relative to their sum = 1 :')
         WRITE(8,48)(WE(I),char(9),R(I),I=1,M)
         write(8,'(/,
     >   ''All data are mol/l of solid solutions and '',
     >   ''mole fractions of endmembers '',
     >   ''unless otherwise noted'',/)')
         ww(1)='Condns.'
         ww(2)='Disapp.'
         ww(3)='T (K)'
         if (ibar.eq.0) then
            ww(4)='P (atm)'
         else
            ww(4)='P (bar)'
         endif
         ww(5)='LogPO2'
         iw=5
      do iss=1,nss
            nmem=ms(iss)
            iflag=0
         do i=1,nmem
            if (itbwflag(nssp(iss)+i-1).eq.1) iflag=1
         enddo
         if (iflag.eq.1) then
            iw=iw+1
            ww(iw)=wss(iss)
            if (nssspinel.ne.0) then
               if (iss.eq.nssspinel)  nmem=nmem-1
            endif
            do i=1,nmem
                  iw=iw+1
                  ww(iw)=wc(nssp(iss)+i-1)
            enddo
            do i=1,nmem
               if ((nonid(iss).ne.0).and.(iss.ne.nssspinel)) then
                  iw=iw+1
                  ww(iw)='g'//wc(nssp(iss)+i-1)
c                  ww(iw)(2:8)=wc(nssp(iss)+i-1)
               endif
            enddo
         endif
      enddo
         write(8,10) (ww(i),char(9),i=1,iw)
      do it=1,ntable
         t0=ttable(it)
          ptot=ptottb(it)
         do i=1,ist
           if (nssspinel.eq.0) then
                p2(nqt(i))=pstable(i,it)*ptot/(rgas*t0)
           else
             isp=nssp(nssspinel)
             if ((nqt(i).eq.isp).or.(nqt(i).eq.isp+1)) then
                p2(nqt(i))=pstable(i,it)
             else
                p2(nqt(i))=pstable(i,it)*ptot/(rgas*t0)
             endif
           endif
         enddo
            ww(1)=wcond(it)
            ww(2)=wdisap(it)
            pw(1)=t0
            pw(2)=pttb(it)
            pw(3)=pO2tb(it)
            iw=3
         do iss=1,nss
               nmem=ms(iss)
               iflag=0
            do i=1,nmem
               if (itbwflag(nssp(iss)+i-1).eq.1) iflag=1
            enddo
            if (iflag.eq.1) then
               pt=0d0
               do i=1,ms(iss)
                  pt=pt+p2(nssp(iss)+i-1)
               enddo
               iw=iw+1
               if (nssspinel.ne.0) then
                  isp=nssp(nssspinel)
                  if (iss.eq.nssspinel) then
                     nmem=nmem-1
                     pw(iw)=p2(isp+2)
                     pt=1d0
                  else
                     pw(iw)=pt
                  endif
               else
                  pw(iw)=pt
               endif
               do i=1,nmem
                  iw=iw+1
                  if (pt.eq.0d0) then
                     pw(iw)=0d0
                  else
                     pw(iw)=p2(nssp(iss)+i-1)/pt
                     sx(iss,i)=pw(iw)
                  endif
               enddo
               do i=1,nmem
                  if ((nonid(iss).ne.0).and.(iss.ne.nssspinel)) then
                     iw=iw+1
                     if (pt.eq.0d0) then
                        pw(iw)=0d0
                     else
                        idf=0
                        call calcgamma(iss,idf)
                        pw(iw)=10**gamma(i)
                     endif
                  endif
               enddo
            endif
         enddo
         write(8,11) (ww(i),char(9),i=1,2),(pw(i),char(9),i=1,iw)
      enddo
      close(8)
 1000 continue
      return
      end

      subroutine openwttables()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      CHARACTER*8  ww(nslimit+26)
      CHARACTER*12 filename
      common/mex/mex
         iex=7
c         mex=5
         iwt=10
      filename= contfile
      i=index(filename,'.')
      j=index(filename,' ')
      if (i.eq.0) then
         if ((j.eq.0).or.(j.gt.9)) then
            filename(9:12)='.WT'
         else
            filename(j:j+3)='.WT'
         endif
      else
         filename(i:i+3)='.WT'
      endif
      OPEN (3,FILE= filename)
      print *,'Making : ',filename
      write(3,'(A10,A50)') 'Program : ',version
      write(3,'(A15,A12)') 'Control File : ',contfile
      write(3,'(A79)') comcntl
      write(3,'(A26,A12)') 'Thermodynamic Data File : ',datafile
      write(3,'(A79)') comdata
      if ((nssspinel.ne.0).and.(nonid(nssspinel).eq.0)) then
      write(3,3002) wc(nssp(nssspinel)+ms(nssspinel))
 3002 format(A8,'was ignored for the calculation of the spinel solid',
     >' solution ')
      endif
      if ((imode.eq.0).or.(imode.eq.2)) then
         write(3,300) dtmin
 300  format('Temperature resolution was set to ',F5.2,' K')
      elseif ((imode.eq.1).or.(imode.eq.3)) then
        if (ibar.eq.0) then
         write(3,302) dpmin
        else
         write(3,1302) dpmin
        endif
 302  format('Pressure resolution was set to ',F5.2,' log10(atm)')
1302  format('Pressure resolution was set to ',F5.2,' log10(bar)')
      endif
      if (ibar.eq.0) then
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(3,101)PTOT
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(3,99)PTOT
      endif
      else
      if ((imode.eq.0).or.(imode.eq.1)) then
        WRITE(3,1101)PTOT
      elseif ((imode.eq.2).or.(imode.eq.3)) then
        WRITE(3,199)PTOT
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
      WRITE(3,48)(WE(I),char(9),R(I),char(9),aw(i),I=1,M)
 48   format(a2,a,1pe12.5,a,e12.5)
      write(3,'(/,''data are mol/l & fraction in condensates,'',
     >         '' and mol% & wt% in oxides'',/)')
      ww(1)='Condns.'
      ww(2)='Disapp.'
      ww(3)='T (K)'
      if (ibar.eq.0) then
         ww(4)='P (atm)'
      else
         ww(4)='P (atm)'
      endif
      ww(5)='RT'
      ww(6)='LogPO2'
      ww(7)='# Solids'
      ww(8)='mol/l'
      iw=8
      mex=0
      do i=1,m
         if (we(i).eq.'Fe') then
            ww(i+iw)='Fe(0)'
            ww(i+iw+1)='Fe(+2)'
            ww(i+iw+2)='Fe(+3)'
            awox(5)=aw(i)+aw(mO)
            awox(6)=aw(i)+aw(mO)*1.5d0
            iw=iw+2
            mex=mex+2
         elseif (we(i).eq.'Si') then
            ww(i+iw)='Si(0)'
            ww(i+iw+1)='Si(+4)'
            awox(4)=aw(i)+aw(mO)*2d0
            iw=iw+1
            mex=mex+1
         elseif (we(i).eq.'Ti') then
            ww(i+iw)='Ti(+3)'
            ww(i+iw+1)='Ti(+4)'
            awox(8)=aw(i)+aw(mO)*1.5d0
            awox(9)=aw(i)+aw(mO)*2d0
            iw=iw+1
            mex=mex+1
         elseif (we(i).eq.'Cr') then
            ww(i+iw)='Cr(0)'
            ww(i+iw+1)='Cr(+3)'
            awox(10)=aw(i)+aw(mO)*1.5d0
            iw=iw+1
            mex=mex+1
         elseif (we(i).eq.'Ca') then
            ww(i+iw)=we(i)
            awox(1)=aw(i)+aw(mO)
         elseif (we(i).eq.'Mg') then
            ww(i+iw)=we(i)
            awox(2)=aw(i)+aw(mO)
         elseif (we(i).eq.'Al') then
            ww(i+iw)=we(i)
            awox(3)=aw(i)+aw(mO)*1.5d0
         elseif (we(i).eq.'Na') then
            ww(i+iw)=we(i)
            awox(7)=aw(i)+aw(mO)*0.5d0
         else
            ww(i+iw)=we(i)
         endif
      enddo
      iw=iex+1+m+mex+1
      ww(iw)='fraction'
      do i=1,m+mex
         ww(iw+i)=ww(iw+i-m-mex-1)
      enddo
      iw=iex+(m+mex+1)*2
      do i=1,2
         if (i.eq.1) ww(iw+1)='mole %'
         if (i.eq.2) ww(iw+1)='wt %'
         iw=iw+1
         ww(iw+1)='CaO'
         ww(iw+2)='MgO'
         ww(iw+3)='1/2Al2O3'
         ww(iw+4)='SiO2'
         ww(iw+5)='FeO'
         ww(iw+6)='1/2Fe2O3'
         ww(iw+7)='1/2Na2O'
         ww(iw+8)='1/2Ti2O3'
         ww(iw+9)='TiO2'
         ww(iw+10)='1/2Cr2O3'
         iw=iw+iwt
      enddo
      write(3,10) (ww(i),char(9),i=1,iw)
 10   format(500(a8,a))
      return
      end

      subroutine wttables()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      dimension pw(nslimit+26)
      CHARACTER*8  ww(2)
      common/mex/mex
      rgas=rgas2
      T=T0*1.d-2
c      TG=(T/TI)**GAM
         iex=5
c         mex=5
         iwt=10
c         iw=iex+1+(m+mex+1)*2+iwt+1+iwt
         iw=nslimit+26
      do i=1,iw
         pw(i)=0d0
      enddo
      it=ntable
         tg=tgtb(it)
         ww(1)=wcond(it)
         ww(2)=wdisap(it)
         pw(1)=t0
         pw(2)=pttb(it)
         pw(3)=rgas*t0
         pw(4)=pO2tb(it)
         pw(5)=dble(istb(it))
         iw=iex+1
         ii=iex+1+(m+mex+1)*2
         i2=m+mex+1
         pw(iw)=0d0
         pw(iw+i2)=0d0
      if (is.gt.0) then
      do j=1,is
         jj=j
         iw=iex+1
         ispflag=0
         if(nssspinel.ne.0) then
            if(j.eq.mss(nssspinel)-m) then
               ispflag=1
               xcr=p(m+n+j)
               xal=1d0-xcr
               xfe=p(M+n+j+1)
               xmg=1d0-xfe
               jj=j+ms(nssspinel)-1
            elseif((j.eq.mss(nssspinel)-m+1).or.
     >               (j.eq.mss(nssspinel)-m+2)) then
               goto 100
            endif
         endif
         pwt=p(m+n+jj)
         k=nq(n+jj)
      do i=1,m
         if (we(i).eq.'Fe') then
            pw(i+iw)=pw(i+iw)+pwt*nu2(k,1)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*nu2(k,1)/r(i)/tg
            pw(i+iw+2)=pw(i+iw+2)+pwt*nu2(k,3)
            pw(i+iw+2+i2)=pw(i+iw+2+i2)+pwt*nu2(k,3)/r(i)/tg
            if (ispflag.eq.1) then
            pw(i+iw+1)=pw(i+iw+1)+pwt*xfe
            pw(i+iw+1+i2)=pw(i+iw+1+i2)+pwt/r(i)/tg*xfe
              if ((kind(k).eq.1).or.(kind(k).eq.2)) then
                pw(ii+5)=pw(ii+5)+pwt*xfe
                pw(ii+6)=pw(ii+6)+pwt*nu2(k,3)
              endif
            else
            pw(i+iw+1)=pw(i+iw+1)+pwt*nu2(k,2)
            pw(i+iw+1+i2)=pw(i+iw+1+i2)+pwt*nu2(k,2)/r(i)/tg
              if ((kind(k).eq.1).or.(kind(k).eq.2)) then
                pw(ii+5)=pw(ii+5)+pwt*nu2(k,2)
                pw(ii+6)=pw(ii+6)+pwt*nu2(k,3)
              endif
            endif
            iw=iw+2
         elseif (we(i).eq.'Si') then
            pw(i+iw)=pw(i+iw)+pwt*nu2(k,4)
            pw(i+iw+1)=pw(i+iw+1)+pwt*nu2(k,5)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*nu2(k,4)/r(i)/tg
            pw(i+iw+1+i2)=pw(i+iw+1+i2)+pwt*nu2(k,5)/r(i)/tg
            if ((kind(k).eq.1).or.(kind(k).eq.2))
     >          pw(ii+4)=pw(ii+4)+pwt*nu2(k,5)
            iw=iw+1
         elseif (we(i).eq.'Ti') then
            pw(i+iw)=pw(i+iw)+pwt*nu2(k,6)
            pw(i+iw+1)=pw(i+iw+1)+pwt*nu2(k,7)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*nu2(k,6)/r(i)/tg
            pw(i+iw+1+i2)=pw(i+iw+1+i2)+pwt*nu2(k,7)/r(i)/tg
            if ((kind(k).eq.1).or.(kind(k).eq.2)) then
                pw(ii+8)=pw(ii+8)+pwt*nu2(k,6)
                pw(ii+9)=pw(ii+9)+pwt*nu2(k,7)
            endif
            iw=iw+1
         elseif (we(i).eq.'Cr') then
            pw(i+iw)=pw(i+iw)+pwt*nu2(k,8)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*nu2(k,8)/r(i)/tg
            if (ispflag.eq.1) then
            pw(i+iw+1)=pw(i+iw+1)+pwt*xcr*2d0
            pw(i+iw+1+i2)=pw(i+iw+1+i2)+pwt/r(i)/tg*xcr*2d0
            if ((kind(k).eq.1).or.(kind(k).eq.2))
     >          pw(ii+10)=pw(ii+10)+pwt*xcr*2d0
            else
            pw(i+iw+1)=pw(i+iw+1)+pwt*nu2(k,9)
            pw(i+iw+1+i2)=pw(i+iw+1+i2)+pwt*nu2(k,9)/r(i)/tg
            if ((kind(k).eq.1).or.(kind(k).eq.2))
     >          pw(ii+10)=pw(ii+10)+pwt*nu2(k,9)
            endif
            iw=iw+1
         elseif (we(i).eq.'Ca') then
            pw(i+iw)=pw(i+iw)+pwt*dnu(k,i)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*dnu(k,i)/r(i)/tg
            if ((kind(k).eq.1).or.(kind(k).eq.2))
     >          pw(ii+1)=pw(ii+1)+pwt*dnu(k,i)
         elseif (we(i).eq.'Mg') then
            if (ispflag.eq.1) then
            pw(i+iw)=pw(i+iw)+pwt*xmg
            pw(i+iw+i2)=pw(i+iw+i2)+pwt/r(i)/tg*xmg
            if ((kind(k).eq.1).or.(kind(k).eq.2))
     >          pw(ii+2)=pw(ii+2)+pwt**xmg
            else
            pw(i+iw)=pw(i+iw)+pwt*dnu(k,i)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*dnu(k,i)/r(i)/tg
            if ((kind(k).eq.1).or.(kind(k).eq.2))
     >          pw(ii+2)=pw(ii+2)+pwt*dnu(k,i)
            endif
         elseif (we(i).eq.'Al') then
            if (ispflag.eq.1) then
            pw(i+iw)=pw(i+iw)+pwt*xal*2d0
            pw(i+iw+i2)=pw(i+iw+i2)+pwt/r(i)/tg*xal*2d0
            if ((kind(k).eq.1).or.(kind(k).eq.2))
     >          pw(ii+3)=pw(ii+3)+pwt*xal*2d0
            else
            pw(i+iw)=pw(i+iw)+pwt*dnu(k,i)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*dnu(k,i)/r(i)/tg
            if ((kind(k).eq.1).or.(kind(k).eq.2))
     >          pw(ii+3)=pw(ii+3)+pwt*dnu(k,i)
            endif
         elseif (we(i).eq.'Na') then
            pw(i+iw)=pw(i+iw)+pwt*dnu(k,i)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*dnu(k,i)/r(i)/tg
            if ((kind(k).eq.1).or.(kind(k).eq.2))
     >          pw(ii+7)=pw(ii+7)+pwt*dnu(k,i)
         elseif (we(i).eq.'O ') then
            if (ispflag.eq.1) then
            pw(i+iw)=pw(i+iw)+pwt*4d0
            pw(i+iw+i2)=pw(i+iw+i2)+pwt/r(i)/tg*4d0
            else
            pw(i+iw)=pw(i+iw)+pwt*dnu(k,i)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*dnu(k,i)/r(i)/tg
            endif
         else
            pw(i+iw)=pw(i+iw)+pwt*dnu(k,i)
            pw(i+iw+i2)=pw(i+iw+i2)+pwt*dnu(k,i)/r(i)/tg
         endif
      enddo
 100  continue
      enddo
         ii=iex+1
         do i=1+ii,m+mex+ii
            pw(i)=pw(i)*ptot/(rgas*t0)
         enddo
         ii=iex+1+(m+mex+1)*2
         pw(ii)=0d0
         pwt=0d0
         do i=1+ii,iwt+ii
            pwt=pwt+pw(i)
         enddo
         if (pwt.ne.0d0) then
         do i=1+ii,iwt+ii
            pw(i)=pw(i)/pwt*100d0
         enddo
         ii=iex+1+(m+mex+1)*2+iwt+1
         pw(ii)=0d0
         pwt=0d0
         do i=1+ii,iwt+ii
            pw(i)=pw(i-iwt-1)*awox(i-ii)
            pwt=pwt+pw(i)
         enddo
         do i=1+ii,iwt+ii
            pw(i)=pw(i)/pwt*100d0
         enddo
         endif
         iw=iex+1+(m+mex+1)*2+iwt+1+iwt
      endif
         write(3,11) (ww(i),char(9),i=1,2),(pw(i),char(9),i=1,iw)
 11      format(2(a8,a),500(1pE12.5,a))
      return
      end
