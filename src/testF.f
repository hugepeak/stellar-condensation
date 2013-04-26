C     Condensation Calculation Program "cwin1123" for Windows
C     "con1123A.for" source file; part one of six (A,B,C,D,E and F)
C     Designed by          L. Grossman  GCA 36,597-619         (1972)
C     Original program by  J. Lattimer  Ap. J. 219,230-249     (1978)
C     Modified by          S. Yoneda    LPS XXV 1533-1534      (1994)
C                                       Meteoritics 29,554-555 (1994)
C                                       GCA submitted          (1994)
C
C $DEFINE mswindows
C $UNDEFINE mswindows
C $IF DEFINED (mswindows)
C      INCLUDE 'FLIB.FI'
C $ENDIF
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      CHARACTER*8  nw(2000)
      dimension pw(5)
C
      CHARACTER*80 A80
      CHARACTER*20 filename,fileele,files(nfilelimit)
C      character*32 arg
      parameter (nzmax = 100, nmax=5000, nsteps = 100)

      integer*4 i, istep, ispecies, iz(nmax), ia(nmax), i_thermo
      real*8 yel(nzmax), y(nmax),matwt,totmol,atno(nmax)

      character*50 snet, szone, sxpath
      logical do_decay

C
      data version/'cwin1123.exe 02/13/95'/
C
      data disappear/1d-11/
      data sscheck,ss2check,factor/1.0d0,0.99d0,1d-4/
      data eps/1.0D-12/
      data m2wt,nfeas/9,4/
      data rgas1/8.314510d0/
C      integer narg, i
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
C $IF DEFINED (mswindows)
      iwindow=1
C $ELSE
C      iwindow=0
C $ENDIF
C
C     Main
C
C      narg = iargc()
C      call getarg( i, arg )

      
      snet = '../data_pub/nuclear_decay_net.xml'
      sxpath = '[z <= 30]'
      szone = '../data_pub/zone.xml'

      do_decay = .true.

      print *, 'Program : ',version
      WRITE(*,'(A)') ' ENTER Control File Name : '
C     Control file name is usually "????????.CON"
      READ (*,'(A)') contfile
      filename=contfile
      i=index(filename,'.')
      j=index(filename,' ')
      if ((i.gt.32).or.
     >    ((i.eq.0).and.(j.ge.35)).or.
     >    (j-i.gt.4)) then
         write(*,*) filename,' --- name is too long !'
         stop
      endif

      OPEN (9,FILE=filename)
      read(9,'(A80)') a80
      read(a80,'(I2)') idebug
      if (idebug.lt.9) then
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
      
C      OPEN (9,FILE='super.abu')
C      do i=1, 5
C      read (9,*) pstore(i)
C      print *, pstore(i)
C      enddo
C      close(9)
C      pause

      do 5001 ifile=1,nfile
      istopflag=0
      contfile=files(ifile)
      filename=contfile
C
      call readcontfile(filename)
C
      i_thermo = 0
      call readthermofile( i_thermo )
c
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
		atno(	1	) =		1.0
		atno(	2	) =		1.0
		atno(	3	) =		2.0
		atno(	4	) =		3.0
		atno(	5	) =		4.0
		atno(	6	) =		3.0
		atno(	7	) =		4.0
		atno(	8	) =		5.0
		atno(	9	) =		6.0
		atno(	10	) =		7.0
		atno(	11	) =		8.0
		atno(	12	) =		9.0
		atno(	13	) =		4.0
		atno(	14	) =		5.0
		atno(	15	) =		6.0
		atno(	16	) =		7.0
		atno(	17	) =		8.0
		atno(	18	) =		9.0
		atno(	19	) =		10.0
		atno(	20	) =		11.0
		atno(	21	) =		6.0
		atno(	22	) =		7.0
		atno(	23	) =		8.0
		atno(	24	) =		9.0
		atno(	25	) =		10.0
		atno(	26	) =		11.0
		atno(	27	) =		12.0
		atno(	28	) =		13.0
		atno(	29	) =		14.0
		atno(	30	) =		7.0
		atno(	31	) =		8.0
		atno(	32	) =		9.0
		atno(	33	) =		10.0
		atno(	34	) =		11.0
		atno(	35	) =		12.0
		atno(	36	) =		13.0
		atno(	37	) =		14.0
		atno(	38	) =		15.0
		atno(	39	) =		16.0
		atno(	40	) =		17.0
		atno(	41	) =		18.0
		atno(	42	) =		19.0
		atno(	43	) =		8.0
		atno(	44	) =		9.0
		atno(	45	) =		10.0
		atno(	46	) =		11.0
		atno(	47	) =		12.0
		atno(	48	) =		13.0
		atno(	49	) =		14.0
		atno(	50	) =		15.0
		atno(	51	) =		16.0
		atno(	52	) =		17.0
		atno(	53	) =		18.0
		atno(	54	) =		19.0
		atno(	55	) =		20.0
		atno(	56	) =		10.0
		atno(	57	) =		11.0
		atno(	58	) =		12.0
		atno(	59	) =		13.0
		atno(	60	) =		14.0
		atno(	61	) =		15.0
		atno(	62	) =		16.0
		atno(	63	) =		17.0
		atno(	64	) =		18.0
		atno(	65	) =		19.0
		atno(	66	) =		20.0
		atno(	67	) =		21.0
		atno(	68	) =		22.0
		atno(	69	) =		12.0
		atno(	70	) =		13.0
		atno(	71	) =		14.0
		atno(	72	) =		15.0
		atno(	73	) =		16.0
		atno(	74	) =		17.0
		atno(	75	) =		18.0
		atno(	76	) =		19.0
		atno(	77	) =		20.0
		atno(	78	) =		21.0
		atno(	79	) =		22.0
		atno(	80	) =		23.0
		atno(	81	) =		24.0
		atno(	82	) =		25.0
		atno(	83	) =		26.0
		atno(	84	) =		27.0
		atno(	85	) =		28.0
		atno(	86	) =		29.0
		atno(	87	) =		30.0
		atno(	88	) =		31.0
		atno(	89	) =		32.0
		atno(	90	) =		33.0
		atno(	91	) =		34.0
		atno(	92	) =		14.0
		atno(	93	) =		15.0
		atno(	94	) =		16.0
		atno(	95	) =		17.0
		atno(	96	) =		18.0
		atno(	97	) =		19.0
		atno(	98	) =		20.0
		atno(	99	) =		21.0
		atno(	100	) =		22.0
		atno(	101	) =		23.0
		atno(	102	) =		24.0
		atno(	103	) =		25.0
		atno(	104	) =		26.0
		atno(	105	) =		27.0
		atno(	106	) =		28.0
		atno(	107	) =		29.0
		atno(	108	) =		30.0
		atno(	109	) =		31.0
		atno(	110	) =		32.0
		atno(	111	) =		33.0
		atno(	112	) =		34.0
		atno(	113	) =		35.0
		atno(	114	) =		36.0
		atno(	115	) =		37.0
		atno(	116	) =		38.0
		atno(	117	) =		16.0
		atno(	118	) =		17.0
		atno(	119	) =		18.0
		atno(	120	) =		19.0
		atno(	121	) =		20.0
		atno(	122	) =		21.0
		atno(	123	) =		22.0
		atno(	124	) =		23.0
		atno(	125	) =		24.0
		atno(	126	) =		25.0
		atno(	127	) =		26.0
		atno(	128	) =		27.0
		atno(	129	) =		28.0
		atno(	130	) =		29.0
		atno(	131	) =		30.0
		atno(	132	) =		31.0
		atno(	133	) =		32.0
		atno(	134	) =		33.0
		atno(	135	) =		34.0
		atno(	136	) =		35.0
		atno(	137	) =		36.0
		atno(	138	) =		37.0
		atno(	139	) =		38.0
		atno(	140	) =		39.0
		atno(	141	) =		40.0
		atno(	142	) =		41.0
		atno(	143	) =		17.0
		atno(	144	) =		18.0
		atno(	145	) =		19.0
		atno(	146	) =		20.0
		atno(	147	) =		21.0
		atno(	148	) =		22.0
		atno(	149	) =		23.0
		atno(	150	) =		24.0
		atno(	151	) =		25.0
		atno(	152	) =		26.0
		atno(	153	) =		27.0
		atno(	154	) =		28.0
		atno(	155	) =		29.0
		atno(	156	) =		30.0
		atno(	157	) =		31.0
		atno(	158	) =		32.0
		atno(	159	) =		33.0
		atno(	160	) =		34.0
		atno(	161	) =		35.0
		atno(	162	) =		36.0
		atno(	163	) =		37.0
		atno(	164	) =		38.0
		atno(	165	) =		39.0
		atno(	166	) =		40.0
		atno(	167	) =		41.0
		atno(	168	) =		42.0
		atno(	169	) =		43.0
		atno(	170	) =		44.0
		atno(	171	) =		19.0
		atno(	172	) =		20.0
		atno(	173	) =		21.0
		atno(	174	) =		22.0
		atno(	175	) =		23.0
		atno(	176	) =		24.0
		atno(	177	) =		25.0
		atno(	178	) =		26.0
		atno(	179	) =		27.0
		atno(	180	) =		28.0
		atno(	181	) =		29.0
		atno(	182	) =		30.0
		atno(	183	) =		31.0
		atno(	184	) =		32.0
		atno(	185	) =		33.0
		atno(	186	) =		34.0
		atno(	187	) =		35.0
		atno(	188	) =		36.0
		atno(	189	) =		37.0
		atno(	190	) =		38.0
		atno(	191	) =		39.0
		atno(	192	) =		40.0
		atno(	193	) =		41.0
		atno(	194	) =		42.0
		atno(	195	) =		43.0
		atno(	196	) =		44.0
		atno(	197	) =		45.0
		atno(	198	) =		46.0
		atno(	199	) =		47.0
		atno(	200	) =		21.0
		atno(	201	) =		22.0
		atno(	202	) =		23.0
		atno(	203	) =		24.0
		atno(	204	) =		25.0
		atno(	205	) =		26.0
		atno(	206	) =		27.0
		atno(	207	) =		28.0
		atno(	208	) =		29.0
		atno(	209	) =		30.0
		atno(	210	) =		31.0
		atno(	211	) =		32.0
		atno(	212	) =		33.0
		atno(	213	) =		34.0
		atno(	214	) =		35.0
		atno(	215	) =		36.0
		atno(	216	) =		37.0
		atno(	217	) =		38.0
		atno(	218	) =		39.0
		atno(	219	) =		40.0
		atno(	220	) =		41.0
		atno(	221	) =		42.0
		atno(	222	) =		43.0
		atno(	223	) =		44.0
		atno(	224	) =		45.0
		atno(	225	) =		46.0
		atno(	226	) =		47.0
		atno(	227	) =		48.0
		atno(	228	) =		49.0
		atno(	229	) =		50.0
		atno(	230	) =		51.0
		atno(	231	) =		22.0
		atno(	232	) =		23.0
		atno(	233	) =		24.0
		atno(	234	) =		25.0
		atno(	235	) =		26.0
		atno(	236	) =		27.0
		atno(	237	) =		28.0
		atno(	238	) =		29.0
		atno(	239	) =		30.0
		atno(	240	) =		31.0
		atno(	241	) =		32.0
		atno(	242	) =		33.0
		atno(	243	) =		34.0
		atno(	244	) =		35.0
		atno(	245	) =		36.0
		atno(	246	) =		37.0
		atno(	247	) =		38.0
		atno(	248	) =		39.0
		atno(	249	) =		40.0
		atno(	250	) =		41.0
		atno(	251	) =		42.0
		atno(	252	) =		43.0
		atno(	253	) =		44.0
		atno(	254	) =		45.0
		atno(	255	) =		46.0
		atno(	256	) =		47.0
		atno(	257	) =		48.0
		atno(	258	) =		49.0
		atno(	259	) =		50.0
		atno(	260	) =		51.0
		atno(	261	) =		52.0
		atno(	262	) =		53.0
		atno(	263	) =		54.0
		atno(	264	) =		23.0
		atno(	265	) =		24.0
		atno(	266	) =		25.0
		atno(	267	) =		26.0
		atno(	268	) =		27.0
		atno(	269	) =		28.0
		atno(	270	) =		29.0
		atno(	271	) =		30.0
		atno(	272	) =		31.0
		atno(	273	) =		32.0
		atno(	274	) =		33.0
		atno(	275	) =		34.0
		atno(	276	) =		35.0
		atno(	277	) =		36.0
		atno(	278	) =		37.0
		atno(	279	) =		38.0
		atno(	280	) =		39.0
		atno(	281	) =		40.0
		atno(	282	) =		41.0
		atno(	283	) =		42.0
		atno(	284	) =		43.0
		atno(	285	) =		44.0
		atno(	286	) =		45.0
		atno(	287	) =		46.0
		atno(	288	) =		47.0
		atno(	289	) =		48.0
		atno(	290	) =		49.0
		atno(	291	) =		50.0
		atno(	292	) =		51.0
		atno(	293	) =		52.0
		atno(	294	) =		53.0
		atno(	295	) =		54.0
		atno(	296	) =		55.0
		atno(	297	) =		56.0
		atno(	298	) =		57.0
		atno(	299	) =		24.0
		atno(	300	) =		25.0
		atno(	301	) =		26.0
		atno(	302	) =		27.0
		atno(	303	) =		28.0
		atno(	304	) =		29.0
		atno(	305	) =		30.0
		atno(	306	) =		31.0
		atno(	307	) =		32.0
		atno(	308	) =		33.0
		atno(	309	) =		34.0
		atno(	310	) =		35.0
		atno(	311	) =		36.0
		atno(	312	) =		37.0
		atno(	313	) =		38.0
		atno(	314	) =		39.0
		atno(	315	) =		40.0
		atno(	316	) =		41.0
		atno(	317	) =		42.0
		atno(	318	) =		43.0
		atno(	319	) =		44.0
		atno(	320	) =		45.0
		atno(	321	) =		46.0
		atno(	322	) =		47.0
		atno(	323	) =		48.0
		atno(	324	) =		49.0
		atno(	325	) =		50.0
		atno(	326	) =		51.0
		atno(	327	) =		52.0
		atno(	328	) =		53.0
		atno(	329	) =		54.0
		atno(	330	) =		55.0
		atno(	331	) =		56.0
		atno(	332	) =		57.0
		atno(	333	) =		58.0
		atno(	334	) =		59.0
		atno(	335	) =		60.0
		atno(	336	) =		25.0
		atno(	337	) =		26.0
		atno(	338	) =		27.0
		atno(	339	) =		28.0
		atno(	340	) =		29.0
		atno(	341	) =		30.0
		atno(	342	) =		31.0
		atno(	343	) =		32.0
		atno(	344	) =		33.0
		atno(	345	) =		34.0
		atno(	346	) =		35.0
		atno(	347	) =		36.0
		atno(	348	) =		37.0
		atno(	349	) =		38.0
		atno(	350	) =		39.0
		atno(	351	) =		40.0
		atno(	352	) =		41.0
		atno(	353	) =		42.0
		atno(	354	) =		43.0
		atno(	355	) =		44.0
		atno(	356	) =		45.0
		atno(	357	) =		46.0
		atno(	358	) =		47.0
		atno(	359	) =		48.0
		atno(	360	) =		49.0
		atno(	361	) =		50.0
		atno(	362	) =		51.0
		atno(	363	) =		52.0
		atno(	364	) =		53.0
		atno(	365	) =		54.0
		atno(	366	) =		55.0
		atno(	367	) =		56.0
		atno(	368	) =		57.0
		atno(	369	) =		58.0
		atno(	370	) =		59.0
		atno(	371	) =		60.0
		atno(	372	) =		61.0
		atno(	373	) =		62.0
		atno(	374	) =		63.0
		atno(	375	) =		27.0
		atno(	376	) =		28.0
		atno(	377	) =		29.0
		atno(	378	) =		30.0
		atno(	379	) =		31.0
		atno(	380	) =		32.0
		atno(	381	) =		33.0
		atno(	382	) =		34.0
		atno(	383	) =		35.0
		atno(	384	) =		36.0
		atno(	385	) =		37.0
		atno(	386	) =		38.0
		atno(	387	) =		39.0
		atno(	388	) =		40.0
		atno(	389	) =		41.0
		atno(	390	) =		42.0
		atno(	391	) =		43.0
		atno(	392	) =		44.0
		atno(	393	) =		45.0
		atno(	394	) =		46.0
		atno(	395	) =		47.0
		atno(	396	) =		48.0
		atno(	397	) =		49.0
		atno(	398	) =		50.0
		atno(	399	) =		51.0
		atno(	400	) =		52.0
		atno(	401	) =		53.0
		atno(	402	) =		54.0
		atno(	403	) =		55.0
		atno(	404	) =		56.0
		atno(	405	) =		57.0
		atno(	406	) =		58.0
		atno(	407	) =		59.0
		atno(	408	) =		60.0
		atno(	409	) =		61.0
		atno(	410	) =		62.0
		atno(	411	) =		63.0
		atno(	412	) =		64.0
		atno(	413	) =		65.0
		atno(	414	) =		66.0
		atno(	415	) =		67.0
		atno(	416	) =		29.0
		atno(	417	) =		30.0
		atno(	418	) =		31.0
		atno(	419	) =		32.0
		atno(	420	) =		33.0
		atno(	421	) =		34.0
		atno(	422	) =		35.0
		atno(	423	) =		36.0
		atno(	424	) =		37.0
		atno(	425	) =		38.0
		atno(	426	) =		39.0
		atno(	427	) =		40.0
		atno(	428	) =		41.0
		atno(	429	) =		42.0
		atno(	430	) =		43.0
		atno(	431	) =		44.0
		atno(	432	) =		45.0
		atno(	433	) =		46.0
		atno(	434	) =		47.0
		atno(	435	) =		48.0
		atno(	436	) =		49.0
		atno(	437	) =		50.0
		atno(	438	) =		51.0
		atno(	439	) =		52.0
		atno(	440	) =		53.0
		atno(	441	) =		54.0
		atno(	442	) =		55.0
		atno(	443	) =		56.0
		atno(	444	) =		57.0
		atno(	445	) =		58.0
		atno(	446	) =		59.0
		atno(	447	) =		60.0
		atno(	448	) =		61.0
		atno(	449	) =		62.0
		atno(	450	) =		63.0
		atno(	451	) =		64.0
		atno(	452	) =		65.0
		atno(	453	) =		66.0
		atno(	454	) =		67.0
		atno(	455	) =		68.0
		atno(	456	) =		69.0
		atno(	457	) =		70.0
		atno(	458	) =		30.0
		atno(	459	) =		31.0
		atno(	460	) =		32.0
		atno(	461	) =		33.0
		atno(	462	) =		34.0
		atno(	463	) =		35.0
		atno(	464	) =		36.0
		atno(	465	) =		37.0
		atno(	466	) =		38.0
		atno(	467	) =		39.0
		atno(	468	) =		40.0
		atno(	469	) =		41.0
		atno(	470	) =		42.0
		atno(	471	) =		43.0
		atno(	472	) =		44.0
		atno(	473	) =		45.0
		atno(	474	) =		46.0
		atno(	475	) =		47.0
		atno(	476	) =		48.0
		atno(	477	) =		49.0
		atno(	478	) =		50.0
		atno(	479	) =		51.0
		atno(	480	) =		52.0
		atno(	481	) =		53.0
		atno(	482	) =		54.0
		atno(	483	) =		55.0
		atno(	484	) =		56.0
		atno(	485	) =		57.0
		atno(	486	) =		58.0
		atno(	487	) =		59.0
		atno(	488	) =		60.0
		atno(	489	) =		61.0
		atno(	490	) =		62.0
		atno(	491	) =		63.0
		atno(	492	) =		64.0
		atno(	493	) =		65.0
		atno(	494	) =		66.0
		atno(	495	) =		67.0
		atno(	496	) =		68.0
		atno(	497	) =		69.0
		atno(	498	) =		70.0
		atno(	499	) =		71.0
		atno(	500	) =		72.0
		atno(	501	) =		73.0
		atno(	502	) =		32.0
		atno(	503	) =		33.0
		atno(	504	) =		34.0
		atno(	505	) =		35.0
		atno(	506	) =		36.0
		atno(	507	) =		37.0
		atno(	508	) =		38.0
		atno(	509	) =		39.0
		atno(	510	) =		40.0
		atno(	511	) =		41.0
		atno(	512	) =		42.0
		atno(	513	) =		43.0
		atno(	514	) =		44.0
		atno(	515	) =		45.0
		atno(	516	) =		46.0
		atno(	517	) =		47.0
		atno(	518	) =		48.0
		atno(	519	) =		49.0
		atno(	520	) =		50.0
		atno(	521	) =		51.0
		atno(	522	) =		52.0
		atno(	523	) =		53.0
		atno(	524	) =		54.0
		atno(	525	) =		55.0
		atno(	526	) =		56.0
		atno(	527	) =		57.0
		atno(	528	) =		58.0
		atno(	529	) =		59.0
		atno(	530	) =		60.0
		atno(	531	) =		61.0
		atno(	532	) =		62.0
		atno(	533	) =		63.0
		atno(	534	) =		64.0
		atno(	535	) =		65.0
		atno(	536	) =		66.0
		atno(	537	) =		67.0
		atno(	538	) =		68.0
		atno(	539	) =		69.0
		atno(	540	) =		70.0
		atno(	541	) =		71.0
		atno(	542	) =		72.0
		atno(	543	) =		73.0
		atno(	544	) =		74.0
		atno(	545	) =		75.0
		atno(	546	) =		76.0
		atno(	547	) =		34.0
		atno(	548	) =		35.0
		atno(	549	) =		36.0
		atno(	550	) =		37.0
		atno(	551	) =		38.0
		atno(	552	) =		39.0
		atno(	553	) =		40.0
		atno(	554	) =		41.0
		atno(	555	) =		42.0
		atno(	556	) =		43.0
		atno(	557	) =		44.0
		atno(	558	) =		45.0
		atno(	559	) =		46.0
		atno(	560	) =		47.0
		atno(	561	) =		48.0
		atno(	562	) =		49.0
		atno(	563	) =		50.0
		atno(	564	) =		51.0
		atno(	565	) =		52.0
		atno(	566	) =		53.0
		atno(	567	) =		54.0
		atno(	568	) =		55.0
		atno(	569	) =		56.0
		atno(	570	) =		57.0
		atno(	571	) =		58.0
		atno(	572	) =		59.0
		atno(	573	) =		60.0
		atno(	574	) =		61.0
		atno(	575	) =		62.0
		atno(	576	) =		63.0
		atno(	577	) =		64.0
		atno(	578	) =		65.0
		atno(	579	) =		66.0
		atno(	580	) =		67.0
		atno(	581	) =		68.0
		atno(	582	) =		69.0
		atno(	583	) =		70.0
		atno(	584	) =		71.0
		atno(	585	) =		72.0
		atno(	586	) =		73.0
		atno(	587	) =		74.0
		atno(	588	) =		75.0
		atno(	589	) =		76.0
		atno(	590	) =		77.0
		atno(	591	) =		78.0
		atno(	592	) =		79.0
		atno(	593	) =		80.0
		atno(	594	) =		36.0
		atno(	595	) =		37.0
		atno(	596	) =		38.0
		atno(	597	) =		39.0
		atno(	598	) =		40.0
		atno(	599	) =		41.0
		atno(	600	) =		42.0
		atno(	601	) =		43.0
		atno(	602	) =		44.0
		atno(	603	) =		45.0
		atno(	604	) =		46.0
		atno(	605	) =		47.0
		atno(	606	) =		48.0
		atno(	607	) =		49.0
		atno(	608	) =		50.0
		atno(	609	) =		51.0
		atno(	610	) =		52.0
		atno(	611	) =		53.0
		atno(	612	) =		54.0
		atno(	613	) =		55.0
		atno(	614	) =		56.0
		atno(	615	) =		57.0
		atno(	616	) =		58.0
		atno(	617	) =		59.0
		atno(	618	) =		60.0
		atno(	619	) =		61.0
		atno(	620	) =		62.0
		atno(	621	) =		63.0
		atno(	622	) =		64.0
		atno(	623	) =		65.0
		atno(	624	) =		66.0
		atno(	625	) =		67.0
		atno(	626	) =		68.0
		atno(	627	) =		69.0
		atno(	628	) =		70.0
		atno(	629	) =		71.0
		atno(	630	) =		72.0
		atno(	631	) =		73.0
		atno(	632	) =		74.0
		atno(	633	) =		75.0
		atno(	634	) =		76.0
		atno(	635	) =		77.0
		atno(	636	) =		78.0
		atno(	637	) =		79.0
		atno(	638	) =		80.0
		atno(	639	) =		81.0
		atno(	640	) =		82.0
		atno(	641	) =		83.0
		atno(	642	) =		38.0
		atno(	643	) =		39.0
		atno(	644	) =		40.0
		atno(	645	) =		41.0
		atno(	646	) =		42.0
		atno(	647	) =		43.0
		atno(	648	) =		44.0
		atno(	649	) =		45.0
		atno(	650	) =		46.0
		atno(	651	) =		47.0
		atno(	652	) =		48.0
		atno(	653	) =		49.0
		atno(	654	) =		50.0
		atno(	655	) =		51.0
		atno(	656	) =		52.0
		atno(	657	) =		53.0
		atno(	658	) =		54.0
		atno(	659	) =		55.0
		atno(	660	) =		56.0
		atno(	661	) =		57.0
		atno(	662	) =		58.0
		atno(	663	) =		59.0
		atno(	664	) =		60.0
		atno(	665	) =		61.0
		atno(	666	) =		62.0
		atno(	667	) =		63.0
		atno(	668	) =		64.0
		atno(	669	) =		65.0
		atno(	670	) =		66.0
		atno(	671	) =		67.0
		atno(	672	) =		68.0
		atno(	673	) =		69.0
		atno(	674	) =		70.0
		atno(	675	) =		71.0
		atno(	676	) =		72.0
		atno(	677	) =		73.0
		atno(	678	) =		74.0
		atno(	679	) =		75.0
		atno(	680	) =		76.0
		atno(	681	) =		77.0
		atno(	682	) =		78.0
		atno(	683	) =		79.0
		atno(	684	) =		80.0
		atno(	685	) =		81.0
		atno(	686	) =		82.0
		atno(	687	) =		83.0
		atno(	688	) =		84.0
		atno(	689	) =		85.0
		atno(	690	) =		86.0
		atno(	691	) =		40.0
		atno(	692	) =		41.0
		atno(	693	) =		42.0
		atno(	694	) =		43.0
		atno(	695	) =		44.0
		atno(	696	) =		45.0
		atno(	697	) =		46.0
		atno(	698	) =		47.0
		atno(	699	) =		48.0
		atno(	700	) =		49.0
		atno(	701	) =		50.0
		atno(	702	) =		51.0
		atno(	703	) =		52.0
		atno(	704	) =		53.0
		atno(	705	) =		54.0
		atno(	706	) =		55.0
		atno(	707	) =		56.0
		atno(	708	) =		57.0
		atno(	709	) =		58.0
		atno(	710	) =		59.0
		atno(	711	) =		60.0
		atno(	712	) =		61.0
		atno(	713	) =		62.0
		atno(	714	) =		63.0
		atno(	715	) =		64.0
		atno(	716	) =		65.0
		atno(	717	) =		66.0
		atno(	718	) =		67.0
		atno(	719	) =		68.0
		atno(	720	) =		69.0
		atno(	721	) =		70.0
		atno(	722	) =		71.0
		atno(	723	) =		72.0
		atno(	724	) =		73.0
		atno(	725	) =		74.0
		atno(	726	) =		75.0
		atno(	727	) =		76.0
		atno(	728	) =		77.0
		atno(	729	) =		78.0
		atno(	730	) =		79.0
		atno(	731	) =		80.0
		atno(	732	) =		81.0
		atno(	733	) =		82.0
		atno(	734	) =		83.0
		atno(	735	) =		84.0
		atno(	736	) =		85.0
		atno(	737	) =		86.0
		atno(	738	) =		87.0
		atno(	739	) =		88.0
		atno(	740	) =		89.0
		atno(	741	) =		42.0
		atno(	742	) =		43.0
		atno(	743	) =		44.0
		atno(	744	) =		45.0
		atno(	745	) =		46.0
		atno(	746	) =		47.0
		atno(	747	) =		48.0
		atno(	748	) =		49.0
		atno(	749	) =		50.0
		atno(	750	) =		51.0
		atno(	751	) =		52.0
		atno(	752	) =		53.0
		atno(	753	) =		54.0
		atno(	754	) =		55.0
		atno(	755	) =		56.0
		atno(	756	) =		57.0
		atno(	757	) =		58.0
		atno(	758	) =		59.0
		atno(	759	) =		60.0
		atno(	760	) =		61.0
		atno(	761	) =		62.0
		atno(	762	) =		63.0
		atno(	763	) =		64.0
		atno(	764	) =		65.0
		atno(	765	) =		66.0
		atno(	766	) =		67.0
		atno(	767	) =		68.0
		atno(	768	) =		69.0
		atno(	769	) =		70.0
		atno(	770	) =		71.0
		atno(	771	) =		72.0
		atno(	772	) =		73.0
		atno(	773	) =		74.0
		atno(	774	) =		75.0
		atno(	775	) =		76.0
		atno(	776	) =		77.0
		atno(	777	) =		78.0
		atno(	778	) =		79.0
		atno(	779	) =		80.0
		atno(	780	) =		81.0
		atno(	781	) =		82.0
		atno(	782	) =		83.0
		atno(	783	) =		84.0
		atno(	784	) =		85.0
		atno(	785	) =		86.0
		atno(	786	) =		87.0
		atno(	787	) =		88.0
		atno(	788	) =		89.0
		atno(	789	) =		90.0
		atno(	790	) =		91.0
		atno(	791	) =		92.0
		atno(	792	) =		44.0
		atno(	793	) =		45.0
		atno(	794	) =		46.0
		atno(	795	) =		47.0
		atno(	796	) =		48.0
		atno(	797	) =		49.0
		atno(	798	) =		50.0
		atno(	799	) =		51.0
		atno(	800	) =		52.0
		atno(	801	) =		53.0
		atno(	802	) =		54.0
		atno(	803	) =		55.0
		atno(	804	) =		56.0
		atno(	805	) =		57.0
		atno(	806	) =		58.0
		atno(	807	) =		59.0
		atno(	808	) =		60.0
		atno(	809	) =		61.0
		atno(	810	) =		62.0
		atno(	811	) =		63.0
		atno(	812	) =		64.0
		atno(	813	) =		65.0
		atno(	814	) =		66.0
		atno(	815	) =		67.0
		atno(	816	) =		68.0
		atno(	817	) =		69.0
		atno(	818	) =		70.0
		atno(	819	) =		71.0
		atno(	820	) =		72.0
		atno(	821	) =		73.0
		atno(	822	) =		74.0
		atno(	823	) =		75.0
		atno(	824	) =		76.0
		atno(	825	) =		77.0
		atno(	826	) =		78.0
		atno(	827	) =		79.0
		atno(	828	) =		80.0
		atno(	829	) =		81.0
		atno(	830	) =		82.0
		atno(	831	) =		83.0
		atno(	832	) =		84.0
		atno(	833	) =		85.0
		atno(	834	) =		86.0
		atno(	835	) =		87.0
		atno(	836	) =		88.0
		atno(	837	) =		89.0
		atno(	838	) =		90.0
		atno(	839	) =		91.0
		atno(	840	) =		92.0
		atno(	841	) =		93.0
		atno(	842	) =		94.0
		atno(	843	) =		95.0
		atno(	844	) =		96.0
		atno(	845	) =		46.0
		atno(	846	) =		47.0
		atno(	847	) =		48.0
		atno(	848	) =		49.0
		atno(	849	) =		50.0
		atno(	850	) =		51.0
		atno(	851	) =		52.0
		atno(	852	) =		53.0
		atno(	853	) =		54.0
		atno(	854	) =		55.0
		atno(	855	) =		56.0
		atno(	856	) =		57.0
		atno(	857	) =		58.0
		atno(	858	) =		59.0
		atno(	859	) =		60.0
		atno(	860	) =		61.0
		atno(	861	) =		62.0
		atno(	862	) =		63.0
		atno(	863	) =		64.0
		atno(	864	) =		65.0
		atno(	865	) =		66.0
		atno(	866	) =		67.0
		atno(	867	) =		68.0
		atno(	868	) =		69.0
		atno(	869	) =		70.0
		atno(	870	) =		71.0
		atno(	871	) =		72.0
		atno(	872	) =		73.0
		atno(	873	) =		74.0
		atno(	874	) =		75.0
		atno(	875	) =		76.0
		atno(	876	) =		77.0
		atno(	877	) =		78.0
		atno(	878	) =		79.0
		atno(	879	) =		80.0
		atno(	880	) =		81.0
		atno(	881	) =		82.0
		atno(	882	) =		83.0
		atno(	883	) =		84.0
		atno(	884	) =		85.0
		atno(	885	) =		86.0
		atno(	886	) =		87.0
		atno(	887	) =		88.0
		atno(	888	) =		89.0
		atno(	889	) =		90.0
		atno(	890	) =		91.0
		atno(	891	) =		92.0
		atno(	892	) =		93.0
		atno(	893	) =		94.0
		atno(	894	) =		95.0
		atno(	895	) =		96.0
		atno(	896	) =		97.0
		atno(	897	) =		98.0
		atno(	898	) =		99.0
		atno(	899	) =		48.0
		atno(	900	) =		49.0
		atno(	901	) =		50.0
		atno(	902	) =		51.0
		atno(	903	) =		52.0
		atno(	904	) =		53.0
		atno(	905	) =		54.0
		atno(	906	) =		55.0
		atno(	907	) =		56.0
		atno(	908	) =		57.0
		atno(	909	) =		58.0
		atno(	910	) =		59.0
		atno(	911	) =		60.0
		atno(	912	) =		61.0
		atno(	913	) =		62.0
		atno(	914	) =		63.0
		atno(	915	) =		64.0
		atno(	916	) =		65.0
		atno(	917	) =		66.0
		atno(	918	) =		67.0
		atno(	919	) =		68.0
		atno(	920	) =		69.0
		atno(	921	) =		70.0
		atno(	922	) =		71.0
		atno(	923	) =		72.0
		atno(	924	) =		73.0
		atno(	925	) =		74.0
		atno(	926	) =		75.0
		atno(	927	) =		76.0
		atno(	928	) =		77.0
		atno(	929	) =		78.0
		atno(	930	) =		79.0
		atno(	931	) =		80.0
		atno(	932	) =		81.0
		atno(	933	) =		82.0
		atno(	934	) =		83.0
		atno(	935	) =		84.0
		atno(	936	) =		85.0
		atno(	937	) =		86.0
		atno(	938	) =		87.0
		atno(	939	) =		88.0
		atno(	940	) =		89.0
		atno(	941	) =		90.0
		atno(	942	) =		91.0
		atno(	943	) =		92.0
		atno(	944	) =		93.0
		atno(	945	) =		94.0
		atno(	946	) =		95.0
		atno(	947	) =		96.0
		atno(	948	) =		97.0
		atno(	949	) =		98.0
		atno(	950	) =		99.0
		atno(	951	) =		100.0
		atno(	952	) =		101.0
		atno(	953	) =		102.0
		atno(	954	) =		51.0
		atno(	955	) =		52.0
		atno(	956	) =		53.0
		atno(	957	) =		54.0
		atno(	958	) =		55.0
		atno(	959	) =		56.0
		atno(	960	) =		57.0
		atno(	961	) =		58.0
		atno(	962	) =		59.0
		atno(	963	) =		60.0
		atno(	964	) =		61.0
		atno(	965	) =		62.0
		atno(	966	) =		63.0
		atno(	967	) =		64.0
		atno(	968	) =		65.0
		atno(	969	) =		66.0
		atno(	970	) =		67.0
		atno(	971	) =		68.0
		atno(	972	) =		69.0
		atno(	973	) =		70.0
		atno(	974	) =		71.0
		atno(	975	) =		72.0
		atno(	976	) =		73.0
		atno(	977	) =		74.0
		atno(	978	) =		75.0
		atno(	979	) =		76.0
		atno(	980	) =		77.0
		atno(	981	) =		78.0
		atno(	982	) =		79.0
		atno(	983	) =		80.0
		atno(	984	) =		81.0
		atno(	985	) =		82.0
		atno(	986	) =		83.0
		atno(	987	) =		84.0
		atno(	988	) =		85.0
		atno(	989	) =		86.0
		atno(	990	) =		87.0
		atno(	991	) =		88.0
		atno(	992	) =		89.0
		atno(	993	) =		90.0
		atno(	994	) =		91.0
		atno(	995	) =		92.0
		atno(	996	) =		93.0
		atno(	997	) =		94.0
		atno(	998	) =		95.0
		atno(	999	) =		96.0
		atno(	1000	) =		97.0
		atno(	1001	) =		98.0
		atno(	1002	) =		99.0
		atno(	1003	) =		100.0
		atno(	1004	) =		101.0
		atno(	1005	) =		102.0
		atno(	1006	) =		103.0
		atno(	1007	) =		104.0
		atno(	1008	) =		105.0


	 
      pind= 0
      fileele=contfile
      i=index(fileele,'.')
      j=index(fileele,' ')
      if (i.eq.0) then
         if ((j.eq.0).or.(j.gt.17)) then
            fileele(17:20)='.ele'
         else
            fileele(j:j+3)='.ele'
         endif
      else
         fileele(i:i+3)='.ele'
      endif
C               OPEN (99,FILE='time.tab')
               OPEN (99,FILE=fileele)


      nw(1)=' time'
      nw(2)=' TK'
      nw(3)=' Pbar'
      nw(4)=' Al'
      nw(5)=' C'
      nw(6)=' Ca'
      nw(7)=' Fe'
      nw(8)=' Mg'
      nw(9)=' N'
      nw(10)=' Na'
      nw(11)=' S'
      nw(12)=' Si'
      nw(13)=' Ti'
      nw(14)=' O'
      nw(15)=' H'
      nw(16)=' Ni'
      nw(17)=' Cr'
      nw(18)=' Co'
      nw(19)=' P'
      nw(20)=' Mn'
      nw(21)=' Cl'
      nw(22)=' K'
      nw(23)=' F'
      nw(24)=' He'
      nw(25)=' Ne'
      nw(26)=' Ar'

      
      call initialize( snet, sxpath )

      call getzone( szone )

      call getnumberofspecies(ispecies)

      write(*,10) ispecies
10    format( "Number of species = ", i5 )

      call printnumberofreactions

c...  Get the mass fractions before decay

      call getabundances( iz, ia, y )
      print 8
8     format( "The initial abundances:" )

      do i = 1, ispecies
C         if( y(i).gt.0.) print 13,i, iz(i), ia(i), y(i)
13       format( 3i5,   1pe12.4 )
      enddo

C            write(99,55) (nw(i),i=1,6)
 55   format(6(a8))
      write(99,33) (nw(i),char(9),i=1,26), (i,char(9),i=1,ispecies)
      write(99,33) (nw(i),char(9),i=1,26), (iz(i),char(9),i=1,ispecies)
      write(99,33)  (nw(i),char(9),i=1,26),(ia(i),char(9),i=1,ispecies)
33       format( 26(a8,a), 2000(i5,a))




c...  Get the elemental abundances before decay

      call getelementalabundances( yel )

c...  Print the initial elemental abundances

      print 9
9     format( "The initial elemental abundances:" )

      do i = 1, 30
C          print 137, i, yel(i)
      enddo    

      r(1)=yel(13)
      r(2)=yel(6)
      r(3)=yel(20)
      r(4)=yel(26)
      r(5)=yel(12)
      r(6)=yel(7)
      r(7)=yel(11)
      r(8)=yel(16)
      r(9)=yel(14)
      r(10)=yel(22)
      r(11)=yel(8)
      r(12)=yel(1)
      r(13)=yel(28)
      r(14)=yel(24)
      r(15)=yel(27)
      r(16)=yel(15)
      r(17)=yel(25)
      r(18)=yel(17)
      r(19)=yel(19)
      r(20)=yel(9)
      r(21)=yel(2)
      r(22)=yel(10)
      r(23)=yel(18)

C         pw(1)=time
C         pw(2)=t0
C         pw(3)=ptot
         iw=3
      do i=1,m
         pw(i+iw)=r(i)
       enddo
      iw=iw+m

       write(99,77) (pw(i),char(9),i=1,iw), (y(i),char(9),i=1,ispecies)

C     abundances are normalized to their sum=1
      R0=0d0
      DO I=1,M
      R0=R0+R(I)
      enddo
      DO I=1,M
      R(I)=R(I)/R0
      enddo
      WRITE(8,490)(WE(I),char(9),R(I),I=1,M)

       istep = 0            
       time = 0
       dtime = 0
       timeprev = 0

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
         if( i_thermo .eq. 1 ) time = timeprev
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
              print 110,T0,Tdone,Tnew
                print 111,PTOT,PTOT
                        
            pind= pind+1    
                   
C            t0=pstore(pind)

      if( i_thermo .eq. 0 ) then

C       time = 3.13819-8.24014E-4*t0
C     add 300K:
       time = 3.38539-8.24014E-4*t0
C     add 500K:
C       time = 3.55019-8.24014E-4*t0

C      Kozasa et al 2009 log(time) in days
C      time = 2.96388-0.000377661*t0+0.0000000694069*t0**2.0

       time = 10.0**time
C      Nozawa et al time in days
C       time = (t0/6.169e7)**(1/-1.78)*1.0


       time = time*24.0*3600.0

      endif

C    t0 = temperature ("current T")
C    timeprev = previous value of the time (corresponding to previous T)
C    time = solve for time from t0.

       timeprev = time
       if( i_thermo .eq. 1 )
     &       call gettimefromtemperature( t0, timeprev, time )

       dtime = time - timeprev
c...  Do the decay loop

C      do istep = 1, nsteps
        istep = istep + 1
        if( do_decay ) then
          call decay( dtime )
          call getelementalabundances( yel )
        endif

c...    Print the elemental abundances

        print 136, istep, time
136      format( "For step ",i5, " ,time = ", 1pe12.4 )

        do i = 1, 30
C           print 137,i,yel(i)
137         format( i5, 2(1pe12.4) )
        enddo

C      enddo

c...  Print the final abundances after decay

      call getabundances( iz, ia, y )

       print 138
138    format( "The final abundances:" )
       matwt = 0.d0
       totmol = 0.d0
      do i = 1, ispecies
      totmol = totmol + y(i)
      enddo
      do i = 1, ispecies
         if( y(i).gt.0.) print 13,i, iz(i), ia(i), y(i)
       matwt = matwt + atno(i) * y(i)/totmol
C          print 139,matwt

         
      enddo 
C      call cleanup
          print 139,matwt
139    format( " MATWT = ", 1pe16.10 )
      r(1)=yel(13)
      r(2)=yel(6)
      r(3)=yel(20)
      r(4)=yel(26)
      r(5)=yel(12)
      r(6)=yel(7)
      r(7)=yel(11)
      r(8)=yel(16)
      r(9)=yel(14)
      r(10)=yel(22)
      r(11)=yel(8)
      r(12)=yel(1)
      r(13)=yel(28)
      r(14)=yel(24)
      r(15)=yel(27)
      r(16)=yel(15)
      r(17)=yel(25)
      r(18)=yel(17)
      r(19)=yel(19)
      r(20)=yel(9)
      r(21)=yel(2)
      r(22)=yel(10)
      r(23)=yel(18)

      if( i_thermo .eq. 1 ) then
        call gettotalpressurefromtime( time, PTOT )
      else
           PTOT=-11.84210526+0.00368421*t0
           PTOT=10.0**PTOT
C      PTOT=-10.81169+0.00138*t0+5.74274E-7*t0*t0
C     from Desch add 300K
      PTOT=-11.17478+0.00104*t0+5.74274E-7*t0*t0
C     Oulette et al 2010 add 300K
      PTOT=-12.33683+0.00106*t0+0.00000056804*t0*t0
 

C     f=2e-4 interpolation from Lattimer et al.
      PTOT=-22.0+3.0*DLOG10(t0)

C      Kozasa et al 2009 log(ro)
C      PTOT=-14.83689+0.00119*t0-0.000000227153*t0**2

      PTOT=10.0**PTOT

C     density from Nozawa et al.
C      PTOT=2.0e-14*(1.0/600.0)**(-3.0)
C      PTOT=PTOT*((time/24.0/3600.0)/1.0)**(-3.0)

C     intermediate pressure:
      PTOT=PTOT/1.0
      PTOT=83.14472*t0/(matwt/PTOT)
      endif

      DLOGP=DLOG10(PTOT)



         pw(1)=time
         pw(2)=t0
         pw(3)=ptot
         iw=3
      do i=1,m
         pw(i+iw)=r(i)
       enddo
      iw=iw+m

       write(99,77) (pw(i),char(9),i=1,iw), (y(i),char(9),i=1,ispecies)
 77      format(2026(1pE12.5,a))
      
      do i=1, m
C      read (99,*) pstore(i)
C      print *, pstore(i)
      enddo
C      close(99)
C      pause
C        if(t0.lt.1000.0) then
        do k=1,m
C        r(k)=pstore(k)
C        print *, r(k)
        enddo
C        endif

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




        do k=1,m
        print *, r(k)
        enddo
 490   format(a2,a,1pe12.5)

C      WRITE(8,49)(WE(I),R(I),I=1,M)
C      WRITE(8,490)(matwt)
      WRITE(8,490)(WE(I),char(9),R(I),I=1,M)
C      WRITE(8,13)(i, iz(i), ia(i), y(i), i = 1,ispecies)


C      pause









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

      call checktemp()
C
      if (ifailconv.ne.0) then
          if ((istopflag.eq.1).or.(t0.eq.tstart)) then
             goto 1000
          elseif (dt.lt.dtmin) then
C             print *,' Failed to converge !!!'
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
      INCLUDE 'con0621i.f'
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
c          print *,I
c          print *,-xk(I+ii-1)
          do ii=1,kch
            sx(j,ii)=-xk(I+ii-1)
            DO JJ=1,M
c          print *,dnu(I+ii-1,JJ)
              sx(j,ii)=sx(j,ii)+DLP(JJ)*dnu(I+ii-1,JJ)
            end do
            sx(j,ii)=10d0**sx(j,ii)
          end do
          chek=sx(2,1)
c          print *,(sx(j,ii),ii=1,kch)
c          print *,kch
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
c        print *,chek
        sxt=0d0
        do ii=1,kch
          sxt=sxt+sx(j,ii)
        end do
       if (nonid(j).eq.0) then
        if (idebug.ge.1) print 510,wss(j),sxt*1d2,
     >   (sx(j,ii)*1d2,ii=1,kch)
        if (idebug.ge.2) write(8,510)wss(j),sxt*1d2,
     >   (sx(j,ii)*1d2,ii=1,kch)
        if (iwindow.ge.1) write(16,510)wss(j),sxt*1d2,
     >                     (sx(j,ii)*1d2,ii=1,kch)
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
C 510  format(1x,A8,' total=',1pe9.2,5(1x,0pf7.4))
 510  format(1X,A8,' total=',1pe11.4,
     >    ' X=',10(0pf9.4,:))
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
c        print * 
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
      INCLUDE 'con0621i.f'
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
      INCLUDE 'con0621i.f'
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
      IF(NITER-1000) 50,50,300
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
      INCLUDE 'con0621i.f'
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
      INCLUDE 'con0621i.f'
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
