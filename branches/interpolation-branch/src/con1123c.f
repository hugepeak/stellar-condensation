C     'con1123c.f' 11/23/94
      subroutine checkposition()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
C     check and store positions of some elements, species, etc. according to
C     their names.
      do i=1,m
         if(we(i).eq.'Fe') mfe=i
         if(we(i).eq.'Mg') mmg=i
         if(we(i).eq.'O ') mO =i
         if(we(i).eq.'Al') mal=i
         if(we(i).eq.'Si') msi=i
         if(we(i).eq.'Cr') mcr=i
         if(we(i).eq.'Ni') mni=i
      enddo
      do j=1,ns+msExtra
         if(wc(j)(1:3).eq.'O2 ') nO2=j
         if(wc(j)(1:3).eq.'Si ') nsi=j
         if(wc(j)(1:3).eq.'Cr ') ncr=j
         if(wc(j)(1:3).eq.'Ni ') nni=j
         if(wc(j)(1:3).eq.'Ab ') nab=j
         if(wc(j)(1:3).eq.'An ') nan=j
         if(wc(j)(1:3).eq.'Fa ') nfa=j
         if(wc(j)(1:3).eq.'Fo ') nfo=j
         if(wc(j)(1:3).eq.'Fs ') nfs=j
         if(wc(j)(1:3).eq.'Ort') noen=j
      enddo
      nssolivine=0
      nssortho=0
      nssmetal=0
      nssspinel=0
      nssfeldspar=0
      nsscmas=0
      nssalsp=0
      if (nss.ne.0) then
         do i=1,nss
            if(wss(i)(1:8).eq.'Olivine ') then
               nssolivine=i
            elseif(wss(i)(1:8).eq.'Ortho En') then
               nssortho=i
            elseif(wss(i)(1:8).eq.'Metal   ') then
               nssmetal=i
            elseif(wss(i)(1:8).eq.'Spinel  ') then
               nssspinel=i
            elseif(wss(i)(1:8).eq.'Feldspar') then
               nssfeldspar=i
            elseif(wss(i)(1:8).eq.'CMAS liq') then
               nsscmas=i
            elseif(wss(i)(1:8).eq.'AlSpinel') then
               nssalsp=i
            endif
         enddo
      endif
      if (nssolivine.ne.0) then
         nfa=nfa-nssp(nssolivine)+1
         nfo=nfo-nssp(nssolivine)+1
      endif
      if (nssortho.ne.0) then
         nfs=nfs-nssp(nssortho)+1
         noen=noen-nssp(nssortho)+1
      endif
      if (nssmetal.ne.0) then
         nsi=nsi-nssp(nssmetal)+1
         ncr=ncr-nssp(nssmetal)+1
         nni=nni-nssp(nssmetal)+1
      endif
      if (nssfeldspar.ne.0) then
         nab=nab-nssp(nssfeldspar)+1
         nan=nan-nssp(nssfeldspar)+1
      endif
      if (nsscmas.ne.0) then
         call CMASopen()
      endif
      if (nssspinel.eq.0) goto 3010
c      data wsp/'FeCr2O4','MgCr2O4','FeAl2O4','MgAl2O4'/
c     ispinel represents the real order in the data set
c      ispinel(4) for MgAl2O4
c      ispinel(3) for FeAl2O4
c      ispinel(2) for MgCr2O4
c      ispinel(1) for FeCr2O4
c      the last species in the data set will be ignored in the calculation
      wsp(1)='FeCr2O4'
      wsp(2)='MgCr2O4'
      wsp(3)='FeAl2O4'
      wsp(4)='MgAl2O4'
      do isp=1,4
         ispinel(isp)=0
      enddo
      do j=nssp(nssspinel),nssp(nssspinel)
     >                       +ms(nssspinel)+mse(nssspinel)-1
        do isp=1,4
        if(wc(j)(1:7).eq.wsp(isp)) ispinel(isp)=j-nssp(nssspinel)+1
        enddo
      enddo
      do 3003 isp=4,1,-1
      print 3004, wsp(isp),ispinel(isp)
 3003 continue
 3004 format(3X,A8,3X,I2)
      if (nonid(nssspinel).eq.0) then
        print 3005,wc(nssp(nssspinel)+ms(nssspinel))
        write(8,3002) wc(nssp(nssspinel)+ms(nssspinel))
        if(isummary.ge.1)write(7,3002)wc(nssp(nssspinel)+ms(nssspinel))
 3002 format(A8,'was ignored for the calculation of the spinel solid',
     >' solution ')
 3005 format(1X,A8,'was ignored for the calculation of the spinel solid'
     >,' solution ')
      endif
      do isp=1,4
        if(ispinel(isp).eq.0) then
          print *,wsp(isp),' is missing in the thermodynamic data file'
          stop
        endif
      enddo
      wc(nssp(nssspinel))='X(Cr)sp'
      wc(nssp(nssspinel)+1)='X(Fe)sp'
      wc(nssp(nssspinel)+2)='Spinel'
 3010 continue
      return
      end

      subroutine readseed()
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*8  wtemp
      INCLUDE 'con0621i.f'
      eps2=dsqrt(eps)
      if (nss.ne.0) then
         do i=1,nsslimit
            nseed(i)=1
            do j=1,isslimit
               do k=1,nseedlimit
                  seed(i,j,k)=0d0
                  sxold(i,j,k)=0d0
                  mustsee(i,k)=0
               enddo
               feas(i,j,1)=0d0
               feas(i,j,2)=1d10
               feas(i,j,3)=0d0
               feas(i,j,4)=1d10
            enddo
         enddo
C        seeds for end-members
         do i=1,nss
            if (i.eq.nssspinel) then
               seed(i,1,2)=0d0
               seed(i,2,2)=0d0
               seed(i,1,3)=1d0
               seed(i,2,3)=0d0
               seed(i,1,4)=0d0
               seed(i,2,4)=1d0
               seed(i,1,5)=1d0
               seed(i,2,5)=1d0
               nseed(i)=nseed(i)+4
            elseif (nonid(i).ge.1) then
               do k=1,ms(i)
                  do j=1,ms(i)
                     if (k.eq.j) then
                        seed(i,j,k+1)=1d0
                     else
                        seed(i,j,k+1)=0d0
                     endif
                  enddo
               enddo
               nseed(i)=nseed(i)+ms(i)
            endif
         enddo
C        seeds for the mid-point
         do i=1,nss
            if (i.eq.nssspinel) then
               seed(i,1,nseed(i)+1)=0.5d0
               seed(i,2,nseed(i)+1)=0.5d0
               nseed(i)=nseed(i)+1
            elseif (nonid(i).ge.1) then
               do j=1,ms(i)
                  seed(i,j,nseed(i)+1)=1d0/dble(ms(i))
               enddo
               nseed(i)=nseed(i)+1
            endif
         enddo
         read(9,*) nn
         if (nn.ne.0) then
            do i=1,nn
               read(9,*) wtemp,nseedtemp,nmustsee
               ii=0
               do j=1,nss
                  if (wss(j)(1:8).eq.wtemp) ii=j
               enddo
               if (ii.eq.0) then
                  print *,wtemp,' is not in the current data set.'
                  stop
               else
                  if (nseedtemp+nseed(ii).gt.nseedlimit) then
                     print *,'Exceed nseedlimit!'
                     stop
                  endif
                  if (nmustsee.gt.nseedtemp) nmustsee=nseedtemp
                  do j=1,nseedtemp
                     sxt=0d0
                     read(9,*) (seed(ii,k,j+nseed(ii)),k=1,ms(ii))
                     do k=1,ms(ii)
                        sxt=sxt+seed(ii,k,j+nseed(ii))
                     enddo
                     do k=1,ms(ii)
                      seed(ii,k,j+nseed(ii))=seed(ii,k,j+nseed(ii))/sxt
                     enddo
                  enddo
                  if (nonid(ii).eq.0) then
                     print 10,wss(ii)
                     write(8,10)wss(ii)
  10                 format (1X,A8,' is not a non-ideal solution')
                  else
                     if (nmustsee.ge.1) then
                        do j=1,nmustsee
                           mustsee(ii,j+nseed(ii))=1
                        enddo
                     endif
                     nseed(ii)=nseedtemp+nseed(ii)
                  endif
               endif
            enddo
         endif
         do i=1,nss
            if (nseed(i).gt.1) then
               do j=2,nseed(i)
                  iflag=0
                  do k=1,ms(i)
                     if (seed(i,k,j).le.0d0) then
                        seed(i,k,j)=eps2
                        iflag=1
                     elseif (seed(i,k,j).ge.1d0) then
                        seed(i,k,j)=1d0-eps2
                        iflag=1
                     endif
                  enddo
                  if ((iflag.eq.1).and.(i.ne.nssspinel)) then
                     sxt=0d0
                     do k=1,ms(i)
                        sxt=sxt+seed(i,k,j)
                     enddo
                     do k=1,ms(i)
                        seed(i,k,j)=seed(i,k,j)/sxt
                     enddo
                  endif
               enddo
            endif
         enddo
      endif
      return
      end

      subroutine setmustsee(imustsee)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      k=nseed(imustsee)+1
      if (k.gt.nseedlimit) then
         print *,'Exceed nseedlimit!'
         stop
      endif
      nseed(imustsee)=k
      if (imustsee.eq.nssspinel) then
         seed(imustsee,1,k)=xcr
         seed(imustsee,2,k)=xfe
         seed(imustsee,3,k)=1d0-xcr-xfe
      else
         do j=1,ms(imustsee)
            seed(imustsee,j,k)=sx(imustsee,j)
         enddo
      endif
      mustsee(imustsee,k)=-1
      return
      end

      subroutine unsetmustsee(imustsee)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      k=nseed(imustsee)
      mustsee(imustsee,k)=0
      nseed(imustsee)=k-1
      return
      end

      subroutine calcgamma(iss,idf)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      common/ac/ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9,ac10
C     return log10(gam.) and dlog10(gam.)/dX
C        gam.: activity coeff.;   X: mole fraction
C     variable names are gamma() and dgamma(), respectively
c      print *,(wss(2)(1:8))
      do i=1,ms(iss)+mse(iss)
         gamma(i)=0d0
         do j=1,ms(iss)
            dgamma(i,j)=0d0
         enddo
      enddo
c      print *,(wss(13)(1:8))
      if ((t0.lt.feas(iss,1,3)).or.
     >             (t0.gt.feas(iss,1,4))) return
      if(wss(iss)(1:8).eq.'Olivine ') then
      a=1d0/dlog(10d0)
c      if(t0.eq.1600) then
c           print *,sx(2,1)
c          stop
c          endif
      gamma(nfa)=a*20.334*(1d0-sx(2,nfa))**2d0/(rgas1/1d3)/t0
      gamma(nfo)=a*20.334*(1d0-sx(2,nfo))**2d0/(rgas1/1d3)/t0
c       gamma(nfa)=0d0
c       gamma(nfo)=0d0
c      dgamma(nfa,nfa)=a*2d0*20.334/(rgas1/1d3)/t0*(-1d0+sx(2,nfa))
c      dgamma(nfa,nfo)=a*2d0*20.334/(rgas1/1d3)/t0*(-1d0+sx(2,nfo))
      proverka=5d0
c     Orthopyroxen
      elseif(wss(iss)(1:8).eq.'Ortho En') then
      a=1d0/dlog(10d0)
      gamma(nfs)=a*14.853*(1d0-sx(3,nfs))**2d0/(rgas1/1d3)/t0
      gamma(noen)=a*14.853*(1d0-sx(3,noen))**2d0/(rgas1/1d3)/t0
c      dgamma(nfs,nfs)=a*14.853/(rgas1/1d3)/t0*(-1d0+sx(3,nfs))
c      dgamma(noen,noen)=a*14.853/(rgas1/1d3)/t0*(-1d0+sx(3,noen))
      elseif(wss(iss)(1:8).eq.'Metal   ') then
C     Metal
C cf. L. Grossman, E. Olsen and J.M. Lattimer, Science, 206, 449-451 (1979)
c      data AC1,AC2,AC5,AC6/1.19d0,-70.7d0,-6.3d0,1.83d2/
c         print *,iss
         ac1=1.19d0
         ac2=-70.7d0
         ac5=-6.3d0
         ac6=1.83d2
         T=T0*1.d-2
         T1=1d0/T
         AC3=AC1+AC2*T1
         AC4=AC5+AC6*T1
         AC7=5.29d0-.337d0*T
         AC8=1.343d0-.08d0*T
         AC9=-.676d0+.03d0*T
         AC10=-.454d0+10.86d0*T1
          gamma(nsi)=(AC3+AC4*sx(iss,nsi))
          gamma(ncr)=(AC10*(1d0-21d0*sx(iss,ncr)))
          gamma(nni)=(AC7*(sx(iss,nni)-AC8)**2d0+AC9)
             dgamma(nsi,nsi)=AC4
             dgamma(ncr,ncr)=-21d0*ac10
             dgamma(nni,nni)=2d0*ac7*(sx(iss,nni)-ac8)
         feas(iss,ncr,2)=0.02d0
         feas(iss,nni,2)=0.5d0
      elseif(wss(iss)(1:8).eq.'Feldspar') then
C     Feldspar
C     Al-avoidance model with calorimetric calibration
C     "Thermochemistry of the high structural state plagioclases"
C     R.C. Newton, T.V. Charlu and O.J. Kleppa, GCA 44, 933-941 (1980)
C       mole fraction X
        xab=sx(iss,nab)
        xan=sx(iss,nan)
        a=1d0/dlog(10d0)
        a1=6746.1d0
        a2=2024.7d0
        rgas=rgas3
C       log10(gamma)
        gamma(nab)=dlog10(xab*(2d0-xab)
     >         *dexp(((1d0-xab)**2)/(rgas*t0)*(a1+2d0*xab*(a2-a1))))
        gamma(nan)=dlog10(((1d0+xan)**2)/4d0
     >         *dexp(((1d0-xan)**2)/(rgas*t0)*(a2+2d0*xan*(a1-a2))))
C       dg=dlog10(gamma)/dX
        if (idf.ne.0) then
        dgab=a*(1d0/xab-1d0/(2d0-xab)
     >       -2d0*(1d0-xab)/(rgas*t0)*(a1+2d0*xab*(a2-a1))
     >       +((1d0-xab)**2)/(rgas*t0)*2d0*(a2-a1))
        dgan=a*(2d0/(1d0+xan)
     >       -2d0*(1d0-xan)/(rgas*t0)*(a2+2d0*xan*(a1-a2))
     >       +((1d0-xan)**2)/(rgas*t0)*2d0*(a1-a2))
              dgamma(nab,nab)=dgab
              dgamma(nan,nan)=dgan
        endif
      elseif(wss(iss)(1:8).eq.'CMAS liq') then
C     CMAS liquid
C     Berman (1983) PhD thesis
C     cf. Berman and Brown (1984) GCA 48, 661-678 and
C         de Capitani and Brown (1987) GCA 51, 2639-2652
C     program was written by J. Beckett (Caltech)
        feas(iss,1,3)=1350d0
        if ((t0.lt.feas(iss,1,3)).or.
     >             (t0.gt.feas(iss,1,4))) return
        TK=t0
        call CMASgamma(TK,iss,idf)
      elseif(wss(iss)(1:8).eq.'AlSpinel') then
C     Aluminous Spinel
C     J. Beckett (1994) personal communication;
C     for L. Chamberlin's paper (to be submitted to Am. Mineral.)
C     These activity data are valid only at 1673 K.
        a=1d0/dlog(10d0)
C     Al2O3
        x=1d0-sx(iss,1)
c        gamma(1)=dexp(x*x*(-418.74+1477.3*x-1724.8*x*x+670.42*x*x*x))
c        dgamma(1,1)=-a*x*(-418.74*2d0+1477.3*x*3d0
c     >                     -1724.8*x*x*4d0+670.42*x*x*x*5d0)
        g=x*x*(-323.4761+1147.8456*x-1347.204*x*x
     >                +527.0512*x*x*x)
          gamma(1)=a*g
          dgamma(1,1)=-a*x*(-323.4761*2d0+1147.8456*x*3d0
     >                     -1347.204*x*x*4d0+527.0512*x*x*x*5d0)
        feas(iss,1,2)=0.25d0
C     MgAl2O4
        x=sx(iss,2)
c        gamma(2)=dexp(-87.392+837.4752*x-2634.756*x*x+3777.069*x*x*x
c     >             -2562.813*x*x*x*x+670.416*x*x*x*x*x)
c        dgamma(2,2)=a*(837.4752-2634.756*x*2d0+3777.069*x*x*3d0
c     >             -2562.813*x*x*x*4d0+670.416*x*x*x*x*5d0)
        g=-66.8588+646.9521*x-2045.24*x*x+2944.118*x*x*x
     >             -2006.02*x*x*x*x+527.0501*x*x*x*x*x
          gamma(2)=a*g
          dgamma(2,2)=a*(646.9521-2045.24*x*2d0+2944.118*x*x*3d0
     >             -2006.02*x*x*x*4d0+527.0501*x*x*x*x*5d0)
        feas(iss,2,1)=0.75d0
c        feas(iss,2,2)=1d0
        elseif(wss(iss)(1:8).eq.'Spinel  ') then
C     Spinel (Mg,Fe)(Al,Cr)2O4 solid solution
C     B.J. Wood and J. Nicholls (1978) Contrib. Mineral. Petrol. 66, 389-400
C     S.A.C. Webb and B.J. Wood (1986) Contrib. Mineral. Petrol. 92, 471-480
        deltag=ag0sp(1)+ag0sp(4)-ag0sp(2)-ag0sp(3)
        deltag=deltag/(rgas1/1d3)/t0
        w=2.7d0*4.184d0*2d0/(rgas1/1d3)/t0
        a=1d0/dlog(10d0)
        gamma(1)=a*w*(xal**2)-a*xmg*xal*deltag
        gamma(2)=a*w*(xal**2)+a*xfe*xal*deltag
        gamma(3)=a*w*(xcr**2)+a*xmg*xcr*deltag
        gamma(4)=a*w*(xcr**2)-a*xfe*xcr*deltag
        dgamma(1,1)=-2d0*a*w*xal+a*xmg*deltag
        dgamma(2,1)=-2d0*a*w*xal-a*xfe*deltag
        dgamma(3,1)=2d0*a*w*xcr+a*xmg*deltag
        dgamma(4,1)=2d0*a*w*xcr-a*xfe*deltag
        dgamma(1,2)=a*xal*deltag
        dgamma(2,2)=a*xal*deltag
        dgamma(3,2)=-a*xcr*deltag
        dgamma(4,2)=-a*xcr*deltag
      endif
c      print *,gamma(4)
C
c      print *, ms(iss)
c      print *,(gamma(i),i=1,ms(iss))
      eps2=dsqrt(eps)
      do i=1,ms(iss)
         if (sx(iss,i).ge.1d0+eps2/10d0) then
            gamma(i)=0d0
            do j=1,ms(iss)
               dgamma(i,j)=0d0
            enddo
         endif
      enddo
      do i=1,ms(iss)
         if (gamma(i).gt.100d0) then
            gamma(i)=100d0
            do j=1,ms(iss)
               dgamma(i,j)=0d0
            enddo
         endif
      enddo
      return
      end

      subroutine sxnonideal(iss,iseed)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      common/fpmin/issmin,ineg(isslimit),inflag,Ginf,ffac,xfac,G0
      dimension np(nseedlimit),x(isslimit)
      EXTERNAL my_dfunc,this_func
C     initial estimation of mole fractions for a non-ideal solid solution
C     return sx(iss,*)
c      print*,iss,iseed
      xfac=100d0
      isxold=0
      ires=0
      itmax=200
      issmin=iss
      eps2=dsqrt(eps)
      gtol=eps
      kch=ms(iss)
      ioutfeas(iseed)=0
      isamecomp(iseed)=0
      isxfail=0
      if (iseed.eq.1) then
         idf=0
         call calcgamma(iss,idf)
      endif
      if ((t0.lt.feas(iss,1,3)).or.
     >             (t0.gt.feas(iss,1,4))) then
        isxfail=1
        ioutfeas(iseed)=1
       if (iseed.eq.1) then
        if(idebug.ge.1)print 523,wss(iss),feas(iss,1,3),feas(iss,1,4)
        if(idebug.ge.2)write(8,523)wss(iss),feas(iss,1,3),feas(iss,1,4)
        if(isummary.ge.2)write(15,523)wss(iss),
     >                  feas(iss,1,3),feas(iss,1,4)
        if(iwindow.ge.1)write(16,523)wss(iss),
     >                  feas(iss,1,3),feas(iss,1,4)
       endif
        do i=1,kch
          sx(iss,i)=0d0
        enddo
        return
523     format (' Out of the feasible temperature: ',A8,
     >            F8.2,' -',F8.2,' K')
      endif
      if (iseed.eq.1) then
C     seed #1 obtained by ideal calculation
         if(wss(iss)(1:8).eq.'CMAS liq') then
c            idf=0
c            call calcgamma(iss,idf)
            xmax=0d0
            do i=1,kch
               if (sx(iss,i).gt.xmax) then
                  xmax=sx(iss,i)
                  ii=i
               endif
            enddo
            sx(iss,ii)=sx(iss,ii)/(10**gamma(ii))
         endif
         sxt=0d0
         do i=1,kch
            sxt=sxt+sx(iss,i)
         enddo
         do i=1,kch
            seed(iss,i,1)=sx(iss,i)/sxt
         enddo
      endif
C     calculation for each seed
      if (mustsee(iss,iseed).eq.0) then
            do i=1,kch
               x(i)=seed(iss,i,iseed)*xfac
            enddo
            isxfail=0
            irun=0
C 10         call frprmn(x,kch-1,gtol,icount,fret)
 10         call
     &        dfpmin(x,kch-1,gtol,icount,fret,this_func,my_dfunc,ires)
C           Exception for metal.  See "func"
            if (wss(iss)(1:8).eq.'Metal   ') fret=G0
            xl=1d0
            do i=1,kch-1
               xl=xl-x(i)/xfac
            enddo
            xl=xl*1d2
        if (ires.eq.1) then
          if (isummary.ge.2) then
           write(15,*)' Restart Minimum Search ',
     >                 'Because Slope is Positive'
           write(15,520) icount,fret,(x(i)/xfac*1d2,i=1,kch-1),xl
          endif
          if (iwindow.ge.1) then
           write(16,*)' Restart Minimum Search ',
     >                 'Because Slope is Positive'
           write(16,520) icount,fret,(x(i)/xfac*1d2,i=1,kch-1),xl
          endif
  520     format (1x,' Iter=',I3,' F=',1pe11.4,' X=',10(0pf9.4,:))
          goto 10
        endif
            irun=irun+1
            if (irun.le.1) then
             if (isummary.ge.2) then
               write(15,510)' run #1:',iseed,icount,fret,
     >                     (x(i)/xfac*1d2,i=1,kch-1),xl
             endif
             if (iwindow.ge.1) then
               write(16,510)' run #1:',iseed,icount,fret,
     >                     (x(i)/xfac*1d2,i=1,kch-1),xl
             endif
             goto 10
            endif
            if (icount.gt.itmax) then
               if (idebug.ge.1) then
                  print *,'too many iterations in dfpmin'
               endif
               if (isummary.ge.2) then
                  write(15,*)'too many iterations in dfpmin'
               endif
               if (iwindow.ge.1) then
                  write(16,*)'too many iterations in dfpmin'
               endif
               isxfail=1
            endif
c            print *,wss(2)
            sxt=0d0
            do i=1,kch-1
               sx(iss,i)=x(i)/xfac
               sxt=sxt+sx(iss,i)
            enddo
               sx(iss,kch)=1d0-sxt
            sxt=1d0
            fseed(iseed)=fret
            if (idebug.ge.1) print 510,wss(iss),iseed,icount,fret,
     >                     (sx(iss,i)*1d2,i=1,kch)
            if (idebug.ge.2) write(8,510)wss(iss),iseed,icount,fret,
     >                     (sx(iss,i)*1d2,i=1,kch)
            if (isummary.ge.2) write(15,510)wss(iss),iseed,icount,fret,
     >                     (sx(iss,i)*1d2,i=1,kch)
            if (iwindow.ge.1) write(16,510)wss(iss),iseed,icount,fret,
     >                     (sx(iss,i)*1d2,i=1,kch)
            do i=1,kch
               if ((sx(iss,i).lt.feas(iss,i,1)).or.
     >             (sx(iss,i).gt.feas(iss,i,2))) then
                  isxfail=1
                  ioutfeas(iseed)=1
                  if (idebug.ge.1) print 513,wss(iss),iseed
                  if (idebug.ge.2) write(8,513) wss(iss),iseed
                  if (isummary.ge.2) write(15,513) wss(iss),iseed
                  if (iwindow.ge.1) write(16,513) wss(iss),iseed
513     format (' Out of the feasible region: ',A8,' Seed#=',I3)
               endif
            enddo
         if (isxfail.eq.1) then
            if (isxold.eq.1) then
C**      Disabled 5/31/94
            if (sxold(iss,1,iseed).ne.0d0) then
               do i=1,kch
                  x(i)=sxold(iss,i,iseed)*xfac
               enddo
               isxfail=0
               irun=0
C 11            call frprmn(x,kch-1,gtol,icount,fret)
 11            call
     &          dfpmin(x,kch-1,gtol,icount,fret,this_func,my_dfunc,ires)
               if (wss(iss)(1:8).eq.'Metal   ') fret=G0
        if (ires.eq.1) then
          if (isummary.ge.2) then
           write(15,*)' Restart Minimum Search ',
     >                 'Because Slope is Positive'
           write(15,520) icount,fret,(x(i)/xfac*1d2,i=1,kch-1)
          endif
          if (iwindow.ge.1) then
           write(16,*)' Restart Minimum Search ',
     >                 'Because Slope is Positive'
           write(16,520) icount,fret,(x(i)/xfac*1d2,i=1,kch-1)
          endif
          goto 11
        endif
            irun=irun+1
            if (irun.le.1) then
             if (isummary.ge.2) then
               write(15,510)' run #1:',iseed,icount,fret,
     >                     (x(i)/xfac*1d2,i=1,kch-1)
             endif
             if (iwindow.ge.1) then
               write(16,510)' run #1:',iseed,icount,fret,
     >                     (x(i)/xfac*1d2,i=1,kch-1)
             endif
             goto 11
            endif
            if (icount.gt.itmax) then
               if (idebug.ge.1) then
                  print *,'too many iterations in dfpmin'
               endif
               if (isummary.ge.2) then
                  write(15,*)'too many iterations in dfpmin'
               endif
               if (iwindow.ge.1) then
                  write(16,*)'too many iterations in dfpmin'
               endif
               isxfail=1
            endif
               sxt=0d0
               do i=1,kch-1
                  sx(iss,i)=x(i)/xfac
                  sxt=sxt+sx(iss,i)
               enddo
                  sx(iss,kch)=1d0-sxt
               sxt=1d0
               fseed(iseed)=fret
               if (idebug.ge.1) print 510,wss(iss),iseed,icount,fret,
     >                     (sx(iss,i)*1d2,i=1,kch)
               if (idebug.ge.2) write(8,510)wss(iss),iseed,icount,fret,
     >                     (sx(iss,i)*1d2,i=1,kch)
               if (isummary.ge.2) write(15,510)wss(iss),iseed,
     >                     icount,fret,(sx(iss,i)*1d2,i=1,kch)
               if (iwindow.ge.1) write(16,510)wss(iss),iseed,
     >                     icount,fret,(sx(iss,i)*1d2,i=1,kch)
               do i=1,kch
                  if ((sx(iss,i).lt.feas(iss,i,1)).or.
     >                (sx(iss,i).gt.feas(iss,i,2))) then
                     isxfail=1
                     ioutfeas(iseed)=1
                  if (idebug.ge.1) print 513,wss(iss),iseed
                  if (idebug.ge.2) write(8,513) wss(iss),iseed
                  if (isummary.ge.2) write(15,513) wss(iss),iseed
                  if (iwindow.ge.1) write(16,513) wss(iss),iseed
                  endif
               enddo
            endif
            endif
         endif
         if (isxfail.eq.1) then
            do i=1,kch
               sxold(iss,i,iseed)=0d0
               sx(iss,i)=0d0
            enddo
            if (ioutfeas(iseed).eq.0) ioutfeas(iseed)=-1
         else
            do i=1,kch
               sxold(iss,i,iseed)=sx(iss,i)
            enddo
         endif
      else
         fret=-1d0
         fseed(iseed)=-1d0
            sxt=0d0
            do i=1,kch
               sxt=sxt+seed(iss,i,iseed)
            enddo
            do i=1,kch
               sx(iss,i)=seed(iss,i,iseed)/sxt
               sxold(iss,i,iseed)=seed(iss,i,iseed)/sxt
            enddo
            sxt=1d0
            if (idebug.ge.1) print 516,iseed,wss(iss),iseed,
     >                     sxt,(sx(iss,i)*1d2,i=1,kch)
            if (idebug.ge.2) write(8,516)iseed,wss(iss),iseed,
     >                     sxt,(sx(iss,i)*1d2,i=1,kch)
            if (isummary.ge.2) write(15,516)iseed,wss(iss),iseed,
     >                     sxt,(sx(iss,i)*1d2,i=1,kch)
            if (iwindow.ge.1) write(16,516)iseed,wss(iss),iseed,
     >                     sxt,(sx(iss,i)*1d2,i=1,kch)
 516     format(' Must be checked  Seed#=',I3,/,
     >          1X,A8,I3,' total=',1pe9.2,' X=',10(0pf9.4,:))
      endif
C     Check whether the composition is different from the others
      if ((iseed.ne.1).and.(isxfail.eq.0)) then
         do iseed2=1,iseed-1
            iflag=0
            do i=1,kch
             if (dabs(sx(iss,i)-sxold(iss,i,iseed2)).gt.eps2) iflag=1
            enddo
            if (iflag.eq.0) then
               do i=1,kch
                  sx(iss,i)=0d0
               enddo
               fret=1d0
               isamecomp(iseed)=iseed2
               if(idebug.ge.1) print 514,iseed2,wss(iss),iseed
               if(idebug.ge.2) write(8,514)iseed2,wss(iss),iseed
               if(isummary.ge.2) write(15,514)iseed2,wss(iss),iseed
               if(iwindow.ge.1) write(16,514)iseed2,wss(iss),iseed
               goto 515
            endif
         enddo
      endif
 514  format(' This is the same composition of #=',I3,
     >       ': ',A8,' Seed#=',I3)
 515  continue
C
      if (fret.gt.eps2) then
               do i=1,kch
                  sx(iss,i)=0d0
               enddo
      endif
C
      if (iseed.eq.nseed(iss)) then
       do iseed2=1,nseed(iss)
        if((ioutfeas(iseed2).eq.0).and.(isamecomp(iseed2).eq.0)) then
         sxt=1d0
         icount=0
         fret=fseed(iseed2)
         if (idebug.ge.1) print 510,wss(iss),iseed2,icount,fret,
     >                     (sxold(iss,i,iseed2)*1d2,i=1,kch)
         if (idebug.ge.2) write(8,510)wss(iss),iseed2,icount,fret,
     >                     (sxold(iss,i,iseed2)*1d2,i=1,kch)
         if (isummary.ge.2) write(15,510)wss(iss),iseed2,icount,fret,
     >                     (sxold(iss,i,iseed2)*1d2,i=1,kch)
         if (iwindow.ge.1) write(16,510)wss(iss),iseed2,icount,fret,
     >                     (sxold(iss,i,iseed2)*1d2,i=1,kch)
         ipcount=0
         do iseed3=iseed2,nseed(iss)
          if(isamecomp(iseed3).eq.iseed2) then
            ipcount=ipcount+1
            np(ipcount)=iseed3
          endif
         enddo
         if(ipcount.ge.1) then
          if (idebug.ge.1) print 517,(np(i),i=1,ipcount)
          if (idebug.ge.2) write(8,517)(np(i),i=1,ipcount)
          if (isummary.ge.2) write(15,517)(np(i),i=1,ipcount)
          if (iwindow.ge.1) write(16,517)(np(i),i=1,ipcount)
 517      format (' Same comp. obtained from Seed#=',12I4)
         endif
        endif
       enddo
       fmin=1d20
       iseedmin=0
       ioutflag=1
       icount=0
       sxt=1d0
       do iseed2=1,nseed(iss)
         fret=fseed(iseed2)
         if (fret.lt.fmin) then
            fmin=fret
            iseedmin=iseed2
         endif
         if (idebug.ge.2) print 510,wss(iss),iseed2,icount,fret,
     >                     (sxold(iss,i,iseed2)*1d2,i=1,kch)
         if (idebug.ge.2) write(8,510)wss(iss),iseed2,icount,fret,
     >                     (sxold(iss,i,iseed2)*1d2,i=1,kch)
 510     format(1X,A8,I3,' Iter=',I3,' F=',1pe11.4,
     >    ' X=',10(0pf9.4,:))
         if (ioutfeas(iseed2).eq.0) ioutflag=0
       enddo
       if (iseedmin.ne.0) then
         if (idebug.ge.2) print 511,wss(iss),iseedmin
         if (idebug.ge.2) write(8,511)wss(iss),iseedmin
 511     format(' Min. estimation of ',A8,' obtained from Seed#=',I3)
       elseif (ioutflag.eq.1) then
         print 512,wss(iss)
         write(8,512) wss(iss)
         if (isummary.ge.2) write(15,512) wss(iss)
         if (iwindow.ge.1) write(16,512) wss(iss)
 512     format(1X,'Failed to estimate initial mole fractions of ',A8)
       endif
      endif
      return
      end

      SUBROUTINE hnonideal(MMdim,H,DF,iss,inumerical,idf)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(MMdim),DF(MMdim,MMdim+1)
      INCLUDE 'con0621i.f'
C     non-ideal solid solution
C     add log10(gamma) and their derivatives to H and DF
      K1=MSS(ISS)-1
      a=1d0/dlog(10d0)
      call calcgamma(iss,idf)
      if (iss.eq.nssspinel) then
       do ii=1,4
        if (ispinel(ii).le.3) then
           H(k1+ispinel(ii))=H(k1+ispinel(ii))+gamma(ii)
           if (inumerical.eq.0) then
              DF(k1+ispinel(ii),k1+1)=DF(k1+ispinel(ii),k1+1)
     >                                -dgamma(ii,1)*xcr
              DF(k1+ispinel(ii),k1+2)=DF(k1+ispinel(ii),k1+2)
     >                                -dgamma(ii,2)*xfe
           endif
        endif
       enddo
      else
       do imem=1,ms(iss)
c       print *,proverka
C       mole fraction X, gamma G,  dlog10(gamma(i))/dXj  dgamma(i,j)
        x=sx(iss,imem)
        if(iss.eq.2d0) gamma(imem)=0.5*gamma(imem)
        if(iss.eq.3d0) gamma(imem)=0.5*gamma(imem)
        g=gamma(imem)
C       add +log10(gamma)*mdb to H
          H(K1+imem)=H(k1+imem)+g*dble(mdb(iss))
C       dg(i,j)=dlog10(gamma(i))/dxj                      xi=ln(Pi)
C         =Sum(dlog10(gamma(i))/dXk)(dXk/dxj)
C         (Xi=Pi/sum(Pi)  then  dXi/dxj=-XiXj   dXi/dxi=Xi(1-Xi))
C       add -dg*mbd to DF
        if (inumerical.eq.0) then
          do j=1,ms(iss)
             dg=0d0
             x=sx(iss,j)
             do k=1,ms(iss)
                if (k.eq.j) then
                   dXdx=x*(1d0-x)
                else
                   dXdx=-x*sx(iss,k)
                endif
                dg=dg+dXdx*dgamma(imem,k)
             enddo
             DF(K1+imem,K1+j)=DF(K1+imem,K1+j)-dg*dble(mdb(iss))
          enddo
        endif
       enddo
      endif
      return
      end
