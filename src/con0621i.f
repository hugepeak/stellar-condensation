C     "con0621i.f"  02/13/95   Include File for declarations
C     For dimension:   see  'con0330p.f'
C     m=25: # of elements             n=150: # of gaseous species
C     ns=300: total # of all species  is=40: # of condensates
C     nss=15: # of solid solutions    iss=10: # of endmembers
C
      INCLUDE 'con0330p.f'
C
      CHARACTER*2  we(mlimit)*2,wtb(mlimit)
C                  we(m)          wtb(m)
      CHARACTER*80 wee(mlimit),wc(nslimit),wss(nsslimit),wref(nslimit)
C                  wee(m),wc(ns),wss(nss),wref(ns)
      CHARACTER*7  wsp(4)
C                  wsp(# spinel endmembers)
      CHARACTER*8  wcond(ntlimit),wdisap(ntlimit)
C                  wcond(ntlimit),wdisap(ntlimit)
      CHARACTER*20 contfile,datafile
      CHARACTER*80 comdata,comcntl
      CHARACTER*50 version
      CHARACTER*1 NAMES(33,17)
      COMMON/A/ chek,t0,ptot,TG,m,n,ngp,ns,nss,nsst,IS,NIS,N1,
     >          nspecies(nslimit)
C               nspecies(ns): store species # specified in the thermodynamic
C               data file and coordinate them with internal species # (ns)
      common/b/ R(mlimit),NU(nslimit,mlimit),dnu(nslimit,mlimit),
     >          ms(nsslimit),mss(nsslimit),mdb(nsslimit)
C               R: abundaces of elements; nu&dnu: stoichiometric coeff.
C               ms: # of endmembers of s.s.
C               mss: function position of s.s. (=m+is)
      common/c/ p(mlimit+nlimit+islimit),dlp(mlimit),xk(nslimit),
     >          nq(nlimit+islimit),eps,niter
C               p(m+n+is): partial pressure of species / ptot
C               dlp(m)=log10(p(m)) for monatomic gases
C               xk(ns):-log10(Kf) (Kf: equilibrium constant of formation reaction)
C               nq(n+is): store internal species #
      common/d/ xc(iordlimit,nslimit),XK0(nslimit-nlimit),
     >          XK1(nslimit-nlimit),xc1(nslimit)
C               xc(iord,ns),XK0(ns-n = # solids),XK1(ns-n),xc1(ns)
C               logKf=xc(1,ns)+xc(2,ns)/T+xc(3,ns)/T^2+... (T=temp.(K)/100)
      common/e/ dt,gam,ti
      common/ss/sx(nsslimit,isslimit),nssp(nsslimit),
     >          mse(nsslimit),msExtra
C               sx: mole fraction of endmember of s.s.
C               nssp: internal species # of solid solution
C               mse: # of extra endmembers
c      common/nonideal/nonid(nsslimit),gamma(isslimit),
c     >          dgamma(isslimit,isslimit)
      common/nonideal/gamma(isslimit),dgamma(isslimit,isslimit),
     >          nonid(nsslimit)
C               nonid(nss): flag for non-ideality of solid solution
C               gamma(iss): activity coefficient
C               dgamma(iss,iss): derivative of gamma, cf. subroutine calcgamma
      common/ss2/sxold(nsslimit,isslimit,nseedlimit),fseed(nseedlimit),
     >           seed(nsslimit,isslimit,nseedlimit),nseed(nsslimit)
      COMMON/DIE/KPC(islimit),KPD(islimit),IC,ID,icont
      common/spin/xal,xmg,xfe,xcr,ispinel(4),xksp(4)
      common/names/we,wee,wc,wss,wref,wsp
      common/check/icps,icpd,icss,icpst,icpdt,icsst,
     > isps(islimit),ispd(islimit),isss(islimit),
     > ispst(islimit),ispdt(islimit),issst(islimit)
      common/saves/psaves(3,mlimit+nlimit+islimit),issave(3),
     >nqsave(3,nlimit+islimit),msssave(3,nsslimit),
     >sxsave(3,nsslimit,isslimit),xcrsave(3),xfesave(3),
     >xk0save(3,nslimit-nlimit),xk1save(3,nslimit-nlimit),tgsave(3)
      common/save2/isave2flag
      common/order/iord
      common/position/mfe,msi,mcr,mni,nO2,nsi,ncr,nni,nssmetal,nssspinel
     > ,mmg,mO,mal
      common/feldspar/nssfeldspar,nab,nan
      common/olivine/nssolivine,nfa,nfo
      common/ortho/nssortho,nfs,noen
      common/flag/idebug,isummary
      common/threshold/sscheck,ss2check,factor
      common/disappear/disappear
      common/failconv/ifailconv
      common/sxfail/isxfail
      common/table0/pmtable(mlimit,ntlimit),pstable(istlimit,ntlimit),
     >              nqt(istlimit),ist
C                   pmtable(m,ntlimit),pstable(ist,ntlimit),nqt(ist)
      common/table1/itable,ntable,ttable(ntlimit)
      common/table2/itbflag(mlimit,nslimit),pttb(ntlimit),
     >              pO2tb(ntlimit),istb(ntlimit)
      common/table3/wtb,wcond,wdisap
      common/table4/ptottb(ntlimit),tgtb(ntlimit)
      common/mode/imode
      common/wttbl/kind(nslimit),nu2(nslimit,m2wtlimit),
     >             aw(mlimit),awox(10),m2wt
      common/file/ contfile,datafile,comdata,comcntl
      common/version/version
      common/r/rgas1,rgas2,rgas3
      common/cmas0/nsscmas
      common/cmas1/CPL(5,4),HL(5),SL(5),WQS(5),WQH(5),pstore(100)
      common/cmas2/WBS(3,10),WBH(3,10),WTS(3,10),WTH(3,10)
      common/cmas3/SOLH(33),SOLS(33),CPS(33,4),MINOX(33,4),NAMES
      common/AlSp/nssalsp
      common/feasible/feas(nsslimit,isslimit,4),nfeas
      common/ghost/nssorg,nssghost(nsslimit)
      common/window/iwindow
      common/sxflag/ioutfeas(nseedlimit),isamecomp(nseedlimit)
      common/must/mustsee(nsslimit,nseedlimit)
      common/z/pstop,tstop,pstart,tstart,dpmin,dtmin,dtemp,dpres,dp
      common/FreeE/ag0(mlimit+nslimit),amuss(nsslimit,isslimit),amugas,
     > ptotal,Vtotal,amu(mlimit+nslimit),amusstotal(nsslimit),ag0sp(4)
      common/gas/sxgas(nlimit+mlimit)
      common/bar/ibar
      common/auto/iauto,istopflag
      real*8 dtime,time,timeprev
