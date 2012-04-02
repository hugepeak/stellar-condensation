C       "con0330e.f"  4/7/94
C
C       subroutines from "B83LIQ.FOR" written by J. Beckett (1992?)
C       and interface routines to the main condensation program
C
C       Program B83LIQ computes activities and chemical potentials
C       of oxide components CaO, MgO, Al2O3, SiO2, and TiO2 in
C       melts given temperature (deg K) and wt % of the oxides.
C       Activities and chemical potentials relative to solid
C       oxides specified by the user are also calculated. 
C       Finally, free energies of reaction between solids
C       for which Berman (1983) gave data, and CMAS melts
C       are calculated for user specified phases.
C       In all cases, oxide components are designated in the
C       order (1) CaO, (2) MgO, (3) Al2O3, (4) SiO2, (5) TiO2.
C       CMAS thermo data is from Berman's 1983 thesis.
C       Thermo data for solids is presented in file BERMSOL in
C       the order (1) alpha-quartz, (2) beta-quartz, (3) beta-
C       tridymite, (4) beta-cristobalite, (5) corundum,
C       (6) lime, (7) periclase, (8) andalusite, 
C       (9) kyanite, (10) sillimanite, (11) pseudowollastonite,
C       (12) wollastonite, (13) rankinite, (14) alpha-larnite,
C       (15) alpha'-larnite (16) gamma-larnite, (17) beta-
C       larnite, (18) tricalcium silicate, (19) tricalcium
C       aluminate, (20) calcium aluminate, (21) calcium
C       di-aluminate, (22) hibonite, (23) protoenstatite,
C       (24) forsterite, (25) spinel (26) anorthite,
C       (27) gehlenite, (28) grossularite, (29) monticellite,
C       (30) merwinite, (31) akermanite, (32) diopside,
C       (33) sapphirine.
C       See Berman and Brown (1984) GCA 48, 661-678 and de Capitani
C       and Brown (1987) GCA 51, 2639-2652 for background.
C
      subroutine CMASopen()
C       read data files for the calculation of CMAS liquid activities
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 NAMES(33,17)
      common/cmas1/CPL(5,4),HL(5),SL(5),WQS(5),WQH(5)
      common/cmas2/WBS(3,10),WBH(3,10),WTS(3,10),WTH(3,10)
      common/cmas3/SOLH(33),SOLS(33),CPS(33,4),MINOX(33,4),NAMES
C
C        INB is the run file.  BERMLIQ contains thermo data
C        for the liquid components.  Output is to OUTB.
C        BERMSOL contains thermo data for solids.
C
c      OPEN(10,FILE='INB.DAT',ACCESS='SEQUENTIAL',STATUS='OLD')
      OPEN(11,FILE='BERMLIQ.DAT',ACCESS='SEQUENTIAL',STATUS='OLD')
c      OPEN(12,FILE='OUTB.DAT',ACCESS='SEQUENTIAL',STATUS='NEW')
      OPEN(13,FILE='BERMSOL.DAT',ACCESS='SEQUENTIAL',STATUS='OLD')
C
C         Read in thermo data for liquid components. Units in Joules.
C         WBS(I,J), WBH(I,J) give binary Margules parameters of types
C         I= (1) Wiiij, (2) Wijjj, (3) Wiijj with J for each in
C         ascending order with i<j.  So WBS(1,J) has W1112, W1113,
C         W1114, W1115, W2223, W2224, W2225, W3334, W3335, W4445 for
C         J=1-10.  WTS, WTH give Margules parameters for ternary
C         subsystems in order iijk, ijjk, ijkk also in ascending order
C         with i<j<k.  WQS, WQH give Margules parameters for
C         quaternary subsystems with i<j<k<l in ascending order.
C         SL and HL give the enthalpy of formation from the 
C         elements and third law entropy, both referenced to
C         298.15K and 1 bar.  CPL gives
C         heat capacities for melt components in the same order
C         for a formulation of CP=a+cT**-2+dT**-1/2+fT**-1.
C         SOLH,SOLS are H298 from the elements and the third
C         law entropy at 298K.  CPS = heat capacity terms for solids
C         in the same form as for CPL.  IDMIN gives stoichiometric
C         coefficients of the CMAS oxides in each of the
C         phases in BERMSOL.
C
      READ(11,*)(WBS(1,I),I=1,10)
      READ(11,*)(WBS(2,I),I=1,10)
      READ(11,*)(WBS(3,I),I=1,10)
      READ(11,*)(WBH(1,I),I=1,5)
      READ(11,*)(WBH(1,I),I=6,10)
      READ(11,*)(WBH(2,I),I=1,5)
      READ(11,*)(WBH(2,I),I=6,10)
      READ(11,*)(WBH(3,I),I=1,5)
      READ(11,*)(WBH(3,I),I=6,10)
      READ(11,*)(WTS(1,I),I=1,10)
      READ(11,*)(WTS(2,I),I=1,10)
      READ(11,*)(WTS(3,I),I=1,10)
      READ(11,*)(WTH(1,I),I=1,5)
      READ(11,*)(WTH(1,I),I=6,10)
      READ(11,*)(WTH(2,I),I=1,5)
      READ(11,*)(WTH(2,I),I=6,10)
      READ(11,*)(WTH(3,I),I=1,5)
      READ(11,*)(WTH(3,I),I=6,10)
      READ(11,*)(WQS(I),I=1,5)
      READ(11,*)(WQH(I),I=1,5)
      DO 200 J=1,5
        READ(11,*)SL(J),HL(J),(CPL(J,I),I=1,4)
 200  CONTINUE
      DO 201 J=1,33
        READ(13,*)(NAMES(J,I),I=1,17)
        READ(13,*)SOLH(J),SOLS(J),(CPS(J,I),I=1,4),(MINOX(J,I),I=1,4)
 201  CONTINUE
      close(11)
      close(13)
      return
      end

      subroutine CMASgamma(TK,iss,idf)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0330p.f'
      REAL*8 X(5),ACT(5),ZMUL(5),ZMUL0(5)
C      REAL*8 COEF(5),BSTOR(3,10),TSTOR(3,10),QSTOR(5)
      REAL*8 BSTOR(3,10),TSTOR(3,10),QSTOR(5)
      REAL*8 ZMUS(4),ACTS(4)
      INTEGER IDSOL(4)
      CHARACTER*1 NAMES(33,17)
      common/ss/sx(nsslimit,isslimit),nssp(nsslimit),
     >          mse(nsslimit),msExtra
      common/nonideal/nonid(nsslimit),gamma(isslimit),
     >          dgamma(isslimit,isslimit)
      common/r/rgas1,rgas2,rgas3
      common/cmas1/CPL(5,4),HL(5),SL(5),WQS(5),WQH(5)
      common/cmas2/WBS(3,10),WBH(3,10),WTS(3,10),WTH(3,10)
      common/cmas3/SOLH(33),SOLS(33),CPS(33,4),MINOX(33,4),NAMES
c      DATA R,T0/8.314,298.15/
      R=rgas1
      T0=298.15d0
C
C         Read in number of liquids to be done and start processing.
C         If NTER=1, calculate ternary phase diagram over user
C         specified composition range.  If NPROJ=1, calculate
C         phase saturated liquidus surface.  If you want to do
C         projections, set NPROJ=1, NTER=1; NTEMP=0
C         If NPRNT=1,suppress printing of activities and just 
C         give equilibration temperatures for user specified
C         solids.  If NTEMP=1, don't calculate any temperatures.
C         TCONV gives convergence level for Teq calculations in deg.
C
C      READ(10,*)NUMLIQ,NTER,NPROJ,NPRNT,NTEMP,TCONV
C
C       Perform requested calculations for each liquid composition
C
      NPRNT=1
      NPROJ=0
C
C         Individual analyses are to be read in from INB.DAT.
C         Read in temperature, title and analysis (wt %).
C         in the order CaO, MgO, Al2O3, SiO2, TiO2.
C         IDSOL gives ID#'s of solid components to be used for
C         solid standard states of the oxides CaO, MgO, Al2O3,
C         SiO2.  In practice, the first three will always be 
C         6,7,5 as long as you restrict yourself to Berman's
C         data table.  Berman gives data for four polymorphs of
C         SiO2.  If you wish to explore alternatives, add
C         the appropriate thermo data to BERMSOL and increase
C         the size of arrays NAMES,SOLH,SOLS,CPS accordingly.
C         NMIN gives the number of minerals for which Grxn is
C         desired and IDMIN gives the ID#'s in BERMSOL of the
C         NMIN minerals.
C         Compute temperature integration of standard state
C         chemical potentials and mole fractions of liquid
C         components.  Then do activities and complete
C         chemical potentials based on liquid and solid standard
C         states.  Finally calculate free energies of reaction
C         for specified minerals.
C
C      READ(10,1007)NMIN,TITLE
C 1007  FORMAT(I2,50A1)
C      READ(10,*)TK,(WTOX(I),I=1,5),IDSOL,(IDMIN(I),I=1,NMIN)
C       WRITE(12,1000)TITLE,TK
C 1000   FORMAT(/1X,50A1,'at ',F8.2,' deg K')
c      solid standard states of the oxides CaO, MgO, Al2O3, and SiO2:
c        lime, periclase, corundum, and beta-cristobalite
       IDSOL(1)=6
       IDSOL(2)=7
       IDSOL(3)=5
       IDSOL(4)=4
C
C       IF(NPRNT.NE.1)THEN
C        WRITE(12,1009)((NAMES(IDSOL(I),J),J=1,17),I=1,4)
C 1009   FORMAT(1X,'solid stnd states =',2(1X,17A1)/1X,19X,2(1X,17A1))
C       WRITE(12,1001)
C 1001   FORMAT(19X,'CaO',9X,'MgO',8X,'Al2O3',7X,'SiO2',8X,'TiO2')
C       WRITE(12,1002)WTOX
C 1002   FORMAT(1X,'Oxide wt %',2X,5E12.5)
C       ENDIF
C
C        Convert wt%'s to mole fractions
C
C      CALL CONVRT(WTOX,X)
C
C       Go through calculations for user specified T if desired.  
C       Output all information.  Then get equilibrium liquidus T's
C       for user specified phases.
C
C      IF(NPRNT.NE.1.AND.NPROJ.NE.1)THEN
C         WRITE(12,1003)X
c 1003     FORMAT(1X,'mole frac ',2X,5e12.5)
C      ENDIF
      do i=1,4
         x(i)=sx(iss,i)
      enddo
         x(5)=0d0
c      print 1003,(x(i),i=1,4)
      DO 182 I=1,4
        CALL CPINT(HL(I),SL(I),CPL(I,1),CPL(I,2),CPL(I,3),CPL(I,4),
     #             TK,ZMUL0(I))
        CALL CPINT(SOLH(IDSOL(I)),SOLS(IDSOL(I)),CPS(IDSOL(I),1),
     #   CPS(IDSOL(I),2),CPS(IDSOL(I),3),CPS(IDSOL(I),4),TK,ZMUS(I))
 182  CONTINUE
c      DO 202 I=1,NMIN
c         CALL CPINT(SOLH(IDMIN(I)),SOLS(IDMIN(I)),CPS(IDMIN(I),1),
c     #   CPS(IDMIN(I),2),CPS(IDMIN(I),3),CPS(IDMIN(I),4),TK,ZMUM(I))
c 202  CONTINUE
      CALL ACT1(BSTOR,TSTOR,QSTOR,WBH,WBS,WTH,WTS,WQH,WQS,X,TK)
      CALL ACT2(BSTOR,TSTOR,QSTOR,X,ACT)
C      CALL ACT2(BSTOR,TSTOR,QSTOR,X,TK,ACT,COEF)
C
C       Output activities, activity coefficients and mu's.
C
c      IF(NPRNT.NE.1.AND.NPROJ.NE.1)THEN
c        WRITE(12,1004)(COEF(I),I=1,4)
c 1004    FORMAT(1X,'act coeff ',2X,5E12.5)
c      ENDIF
      DO 183 I=1,4
c        IF(X(I).LT.1.E-6)THEN
        IF(X(I).LE.0d0)THEN
           ZMUL(I)=0d0
        ELSE
c                 ACT(I)=RTln(gamma(i)liq)
           ZMUL(I)=ZMUL0(I)+R*TK*DLOG(X(I))+ACT(I)
cc           ACT(I)=COEF(I)*X(I)
        ENDIF
 183  CONTINUE
C
C        Calculate activities vs. solid standard states.
C
       DO 203 I=1,4
c         IF(X(I).LT.1.E-6)THEN
         IF(X(I).LE.0d0)THEN
           ACTS(I)=0d0
         ELSE
c           Z0=ZMUL(I)-R*TK*DLOG(ACT(I))
c           ACTS(I)=DEXP((Z0-ZMUS(I))/(R*TK))*ACT(I)
cc           ACTS(I)=DEXP((ZMUL0(I)-ZMUS(I))/(R*TK))*ACT(I)
c                 ACTS(I)=RTln(gamma(i)sol)
           ACTS(I)=(ZMUL0(I)-ZMUS(I))+ACT(I)
         ENDIF
 203   CONTINUE
c       IF(NPRNT.NE.1.AND.NPROJ.NE.1)THEN
c         WRITE(12,1005)(ACT(K),K=1,4)
c 1005     FORMAT(1X,'activity/L',2X,4E12.5)
c         WRITE(12,1006)(ZMUL(K),K=1,4)
c 1006     FORMAT(1X,'mu oxide/L',2X,4E12.5)
c         WRITE(12,1011)(ACTS(I),I=1,4)
c 1011     FORMAT(1X,'activity/S',2X,4E12.5)
c       ENDIF
          a=1d0/dlog(10d0)
c                 ACTS=RTln(g)  gamma=log10(g)
       do i=1,4
          gamma(i)=a*acts(i)/(R*TK)
       enddo
c       print 1004,(gamma(i),i=1,4)
       if (idf.ne.0) then
          do j=1,4
             call CMASdg(BSTOR,TSTOR,QSTOR,X,ACT,j)
             do i=1,4
                dgamma(i,j)=act(i)*a/(R*TK)
             enddo
          enddo
       endif
       return
       END
C
C
C      SUBROUTINE CONVRT(WTOX,X)
C
C        Convert oxide wt %'s CaO, MgO, Al2O3, SiO2, TiO2
C        to mole fractions of the oxides in the same order.
C
C      IMPLICIT REAL(A-H,O-Z)
C      REAL X(5),WTOX(5),MOLWT(5)
C      DATA MOLWT/56.079,40.304,101.961,60.085,79.879/
C      SUM=0.
C      DO 100 I=1,5
C        X(I)=WTOX(I)/MOLWT(I)
C        SUM=SUM+X(I)
C 100  CONTINUE
C      DO 101 I=1,5
C        X(I)=X(I)/SUM
C 101  CONTINUE
C      RETURN
C      END
C
C
      SUBROUTINE CPINT(H,S,A,B,C,D,T,G)
C        Perform Cp integration for phase or component.
C        H is the 298K enthalpy of formation from the elements,
C        S is the third law entropy at 298K, A-D are Cp terms
C        as Cp=A+BT**-2+CT**-1/2+DT**-1.  T is the temperature
C        in deg K; G is the apparent standard state free energy.
C
      IMPLICIT REAL*8(A-H,O-Z)
C      DATA T0/298.15/
      T0=298.15d0
C
C        First get 298K portions of integration.
C
      ZH=H-A*T0+B/T0-2*C*SQRT(T0)-D*DLOG(T0)
      ZS=S-A*DLOG(T0)+B/(2*T0**2)+2*C/SQRT(T0)+D/T0
C
C         Now get temperature portion of integration and compute G.
C
      ZHT=A*T-B/T+2*C*SQRT(T)+D*DLOG(T)
      ZST=A*DLOG(T)-B/(2*T**2)-2*C/SQRT(T)-D/T
      G=ZH+ZHT-T*(ZS+ZST)
      RETURN
      END
C
C
      SUBROUTINE ACT1(BSTOR,TSTOR,QSTOR,WBH,WBS,WTH,WTS,WQH,WQS,X,T)
C
C         Compute individual terms in activity coefficient summation,
C         each exclusive of that part (Qm/Xm-3) dependent on the 
C         particular component whose activity coefficient is
C         to be determined.  Terms involving binary Margules parameters
C         are stored in BSTOR, ternary parameters in TSTOR and
C         quaternary terms in QSTOR.  BSTOR, TSTOR and QSTOR are 
C         accessed by subroutine ACT2 in order to compute specific
C         activity coefficients.
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 BSTOR(3,10),TSTOR(3,10),QSTOR(5),WBS(3,10),WBH(3,10)
      REAL*8 WTS(3,10),WTH(3,10),WQS(5),WQH(5),X(5)
C
C       Set counters
C
      NQT=1
      NBIN=1
      NTER=1
C
C       Start summations.  Note that W=WH-T*WS.
C
      DO 100 I=1,4
      DO 101 J=I+1,5
C
C          Store terms involving binary Margules parameters
C
       BSTOR(1,NBIN)=(WBH(1,NBIN)-T*WBS(1,NBIN))*X(J)*X(I)**3
       BSTOR(2,NBIN)=(WBH(2,NBIN)-T*WBS(2,NBIN))*X(I)*X(J)**3
       BSTOR(3,NBIN)=(WBH(3,NBIN)-T*WBS(3,NBIN))*(X(I)*X(J))**2
       NBIN=NBIN+1
       IF(J.NE.5)THEN
        DO 102 K=J+1,5
C
C            Store terms involving ternary Margules parameters
C
         TSTOR(1,NTER)=(WTH(1,NTER)-T*WTS(1,NTER))*X(J)*X(K)*X(I)**2
         TSTOR(2,NTER)=(WTH(2,NTER)-T*WTS(2,NTER))*X(I)*X(K)*X(J)**2
         TSTOR(3,NTER)=(WTH(3,NTER)-T*WTS(3,NTER))*X(I)*X(J)*X(K)**2
         NTER=NTER+1
         IF(K.NE.5)THEN
          DO 103 L=K+1,5
C
C            Store terms involving quaternary Margules parameters
C
           QSTOR(NQT)=(WQH(NQT)-T*WQS(NQT))*X(I)*X(J)*X(K)*X(L)
           NQT=NQT+1
 103      CONTINUE
         ENDIF
 102    CONTINUE
       ENDIF
 101  CONTINUE
 100  CONTINUE
      RETURN
      END
C
C
C      SUBROUTINE ACT2(BSTOR,TSTOR,QSTOR,X,TK,ACT,ACTC)
      SUBROUTINE ACT2(BSTOR,TSTOR,QSTOR,X,ACT)
C
C       Calculate RTln(gamma i) and gamma(i) for each oxide.
C       BSTOR, TSTOR, AND QSTOR are terms obtained from 
C       subroutine ACT1.  X contains the mole fractions,
C       TK the temperature in degrees K, ACT(I) =
C       RTln(G=gamma(I)) and ACTC(I) = gamma(I).
C
      IMPLICIT  REAL*8(A-H,O-Z)
C      REAL*8 BSTOR(3,10),TSTOR(3,10),QSTOR(5),X(5),ACT(5),ACTC(5)
      REAL*8 BSTOR(3,10),TSTOR(3,10),QSTOR(5),X(5),ACT(5)
      common/r/rgas1,rgas2,rgas3
c      DATA R/8.314/
      R=rgas1
      DO 110 M=1,5
c      IF(X(M).GE.1.E-6)THEN
      IF(X(M).GT.0d0)THEN
C
C       X(M) is nonzero.  Set counters and calculate gamma.
C
      NQT=1
      NBIN=1
      NTER=1
      ACT(M)=0d0
C
C       If M is same as one of the indices i,j,k,l set the
C       appropriate counter MQI, MQJ, MQK or MQL TO 1 
C       (otherwise 0) so that WX terms in BSTOR, TSTOR and
C       QSTOR can be multiplied by the right coefficients.
C       Once Qm is determined for a given type of term,
C       multiply it out and add it to the running sum.
C
      DO 100 I=1,4
      IF(M.EQ.I)THEN
      MQI=1
      ELSE
      MQI=0
      ENDIF
      DO  101 J=I+1,5
       IF(M.EQ.J)THEN
       MQJ=1
       ELSE
       MQJ=0
       ENDIF
C
C        get binary Margules terms
C
      ACT(M)=ACT(M)+BSTOR(1,NBIN)*((3*MQI+MQJ)/X(M)-3d0)
      ACT(M)=ACT(M)+BSTOR(2,NBIN)*((3*MQJ+MQI)/X(M)-3d0)
      ACT(M)=ACT(M)+BSTOR(3,NBIN)*((2*(MQI+MQJ))/X(M)-3d0)
      NBIN=NBIN+1
       IF(J.NE.5)THEN
        DO  102 K=J+1,5
        IF(M.EQ.K)THEN
         MQK=1
         ELSE
         MQK=0
        ENDIF
C
C        get ternary Margules terms
C
        ACT(M)=ACT(M)+TSTOR(1,NTER)*((MQI*2+MQJ+MQK)/X(M)-3d0)
        ACT(M)=ACT(M)+TSTOR(2,NTER)*((MQJ*2+MQI+MQK)/X(M)-3d0)
        ACT(M)=ACT(M)+TSTOR(3,NTER)*((MQK*2+MQI+MQJ)/X(M)-3d0)
        NTER=NTER+1
        IF(K.NE.5)THEN
        DO  103 L=K+1,5
        IF(M.EQ.L)THEN
        MQL=1
        ELSE
        MQL=0
        ENDIF
C
C          get quaternary Margules term
C
        ACT(M)=ACT(M)+QSTOR(NQT)*((MQI+MQJ+MQK+MQL)/X(M)-3d0)
        NQT=NQT+1
 103    CONTINUE
        ENDIF
 102    CONTINUE
        ENDIF
 101    CONTINUE
 100    CONTINUE
c        ACTC(M)=DEXP(ACT(M)/(R*TK))
        ELSE
        ACT(M)=0d0
c        ACTC(M)=0d0
        ENDIF
 110    CONTINUE
        RETURN
        END
C
C
      SUBROUTINE CMASdg(BSTOR,TSTOR,QSTOR,X,ACT,N)
C
C       Calculate dRTln(gamma(m))/dX(n) for each oxide.
C       BSTOR, TSTOR, AND QSTOR are terms obtained from 
C       subroutine ACT1.  X contains the mole fractions,
C       TK the temperature in degrees K, ACT(m) =
C       dRTln(G=gamma(m))/dX(n).
C       This routine assumes that all Xs are independent
C       one another (actually sum(X)=1).
C
      IMPLICIT  REAL*8(A-H,O-Z)
      REAL*8 BSTOR(3,10),TSTOR(3,10),QSTOR(5),X(5),ACT(5)
      if (n.gt.4) then
         print *,'CMAS liquid is defined only Si, Al, Mg, and Si'
         stop
      endif
      DO 110 M=1,5
      if (m.eq.n) then
         mn=1
      else
         mn=0
      endif
c      IF(X(M).GE.1.E-6)THEN
      IF((X(M).GT.0d0).and.(X(n).gt.0d0)) THEN
C
C       X(M) is nonzero.  Set counters and calculate gamma.
C
      NQT=1
      NBIN=1
      NTER=1
      ACT(M)=0d0
C
C       If M is same as one of the indices i,j,k,l set the
C       appropriate counter MQI, MQJ, MQK or MQL TO 1 
C       (otherwise 0) so that WX terms in BSTOR, TSTOR and
C       QSTOR can be multiplied by the right coefficients.
C       Once Qm is determined for a given type of term,
C       multiply it out and add it to the running sum.
C
      DO 100 I=1,4
      IF(M.EQ.I)THEN
      MQI=1
      ELSE
      MQI=0
      ENDIF
      if (i.eq.n) then
         nqi=1
      else
         nqi=0
      endif
      DO  101 J=I+1,5
       IF(M.EQ.J)THEN
       MQJ=1
       ELSE
       MQJ=0
       ENDIF
       if (j.eq.n) then
          nqj=1
       else
          nqj=0
       endif
C
C        get binary Margules terms
C
      ACT(M)=ACT(M)+BSTOR(1,NBIN)*((3*MQI+MQJ-mn)/X(M)-3d0)
     >              *(3*nqi+nqj)/X(n)
      ACT(M)=ACT(M)+BSTOR(2,NBIN)*((3*MQJ+MQI-mn)/X(M)-3d0)
     >              *(3*nqj+nqi)/X(n)
      ACT(M)=ACT(M)+BSTOR(3,NBIN)*((2*(MQI+MQJ)-mn)/X(M)-3d0)
     >              *(2*(nqi+nqj))/X(n)
      NBIN=NBIN+1
       IF(J.NE.5)THEN
        DO  102 K=J+1,5
        IF(M.EQ.K)THEN
         MQK=1
         ELSE
         MQK=0
        ENDIF
        if (k.eq.n) then
           nqk=1
        else
           nqk=0
        endif
C
C        get ternary Margules terms
C
        ACT(M)=ACT(M)+TSTOR(1,NTER)*((MQI*2+MQJ+MQK-mn)/X(M)-3d0)
     >                *(nqi*2+nqj+nqk)/X(n)
        ACT(M)=ACT(M)+TSTOR(2,NTER)*((MQJ*2+MQI+MQK-mn)/X(M)-3d0)
     >                *(nqj*2+nqi+nqk)/X(n)
        ACT(M)=ACT(M)+TSTOR(3,NTER)*((MQK*2+MQI+MQJ-mn)/X(M)-3d0)
     >                *(nqk*2+nqi+nqj)/X(n)
        NTER=NTER+1
        IF(K.NE.5)THEN
        DO  103 L=K+1,5
        IF(M.EQ.L)THEN
        MQL=1
        ELSE
        MQL=0
        ENDIF
        if (l.eq.n) then
           nql=1
        else
           nql=0
        endif
C
C          get quaternary Margules term
C
        ACT(M)=ACT(M)+QSTOR(NQT)*((MQI+MQJ+MQK+MQL-mn)/X(M)-3d0)
     >                *(nqi+nqj+nqk+nql)/X(n)
        NQT=NQT+1
 103    CONTINUE
        ENDIF
 102    CONTINUE
        ENDIF
 101    CONTINUE
 100    CONTINUE
        ELSE
        ACT(M)=0d0
        ENDIF
 110    CONTINUE
        RETURN
        END
