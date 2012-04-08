C     "con1115f.f" 11/16/94
C     Interface routines for the use of "Numerical Recipes"

      DOUBLE PRECISION FUNCTION this_func(x)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION x(*)
      INCLUDE 'con0621i.f'
      common/fpmin/issmin,ineg(isslimit),inflag,Ginf,ffac,xfac,G0
      common/ac/ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9,ac10
      gasconst=rgas1/1d3
      a=dlog(10d0)
      Ginf=1d5
      ffac=1d0
      G=0d0
      G0=0d0
      if (wss(issmin)(1:8).eq.'Spinel  ') then
         xcr=x(1)/xfac
         xfe=x(2)/xfac
            xal=1d0-xcr
            x(3)=xal*xfac
            xmg=1d0-xfe
            x(4)=xmg*xfac
         inflag=0
         do i=1,4
            if (x(i).le.0) then
               ineg(i)=1
               inflag=1
            else
               ineg(i)=0
            endif
         enddo
        if (inflag.eq.0) then
         call ssfreeenergy(issmin,G)
         do i=1,4
            do j=1,m
               G=G-sx(issmin,i)*amu(j)*dnu(nssp(issmin)+ispinel(i)-1,j)
            enddo
         enddo
        else
         G=Ginf
         do i=1,4
            G=G-Ginf*x(i)*ineg(i)
         enddo
        endif
      else
         xt=0d0
         do i=1,ms(issmin)-1
            xt=xt+x(i)
         enddo
         x(ms(issmin))=1d0*xfac-xt
         inflag=0
         do i=1,ms(issmin)
            if (x(i).le.0) then
               ineg(i)=1
               inflag=1
            else
               ineg(i)=0
            endif
         enddo
        if (inflag.eq.0) then
         do i=1,ms(issmin)
            sx(issmin,i)=x(i)/xfac
         enddo
         call ssfreeenergy(issmin,G)
         do i=1,ms(issmin)
            do j=1,m
               G=G-sx(issmin,i)*amu(j)*dnu(nssp(issmin)+i-1,j)
            enddo
         enddo
         if (wss(issmin)(1:8).eq.'Metal   ') then
C        Gammas for metal do not comply with the Gibbs-Duhem equation !
C        Additional term <g*> is needed to obtain the same results from
C        mass action law equations by minimization of this function G.
C        <g*> is defined by simultaneous partial differential equations:
C             (dg*/dXj) = -RT*Sumi(Xi*(dlnGammai/dXj))
C                                           (i= 1 ~ k, j= 1 ~ (k-1))
C        Fortunately, since Gammas for metal are defined only by their
C        corresponding mole fractions, the summation for i vanishes.
C             (dg*/dXj) = -RT*Xj*(dlnGammaj/dXj)      (j= 1 ~ (k-1))
C        then,
C             g*=-RT*Sumj(Integ(Xj*(dlnGammaj/dXj)dXj))
C        However, evaluation of the function must be done without this
C        additional term.
           G0=G
           G=G-a*gasconst*t0*ac4/2d0*sx(issmin,nsi)**2
           G=G-a*gasconst*t0*ac10*(-21d0)/2d0*sx(issmin,ncr)**2
           G=G-a*gasconst*t0*ac7*(sx(issmin,nni)*2d0/3d0-ac8)
     >          *sx(issmin,nni)**2
         endif
        else
         G=Ginf
         do i=1,ms(issmin)
            G=G-Ginf*x(i)*ineg(i)
         enddo
            G0=G
        endif
      endif
      this_func=G/ffac
      G0=G0/ffac
      end

      SUBROUTINE my_dfunc(x,df)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'con0621i.f'
      DOUBLE PRECISION x(*),df(isslimit)
      common/fpmin/issmin,ineg(isslimit),inflag,Ginf,ffac,xfac,G0
      gasconst=rgas1/1d3
      a=dlog(10d0)
      f=this_func(x)
      if (wss(issmin)(1:8).eq.'Spinel  ') then
        if (inflag.eq.0) then
            df(1)=amuss(issmin,2)-amuss(issmin,4)
            df(2)=amuss(issmin,3)-amuss(issmin,4)
            do j=1,m
               df(1)=df(1)-amu(j)*
     >           (dnu(nssp(issmin)+ispinel(1)-1,j)*xfe
     >           +dnu(nssp(issmin)+ispinel(2)-1,j)*xmg
     >           +dnu(nssp(issmin)+ispinel(3)-1,j)*(-xfe)
     >           +dnu(nssp(issmin)+ispinel(4)-1,j)*(-xmg))
               df(2)=df(2)-amu(j)*
     >           (dnu(nssp(issmin)+ispinel(1)-1,j)*xcr
     >           +dnu(nssp(issmin)+ispinel(2)-1,j)*(-xcr)
     >           +dnu(nssp(issmin)+ispinel(3)-1,j)*xal
     >           +dnu(nssp(issmin)+ispinel(4)-1,j)*(-xal))
            enddo
        else
            df(1)=-Ginf*(ineg(1)-ineg(3))
            df(2)=-Ginf*(ineg(2)-ineg(4))
        endif
      else
        if (inflag.eq.0) then
         do i=1,ms(issmin)-1
            df(i)=amuss(issmin,i)-amuss(issmin,ms(issmin))
            do j=1,m
               df(i)=df(i)-amu(j)*dnu(nssp(issmin)+i-1,j)
     >           +amu(j)*dnu(nssp(issmin)+ms(issmin)-1,j)
            enddo
         enddo
        else
         do i=1,ms(issmin)-1
            df(i)=-Ginf*(ineg(i)-ineg(ms(issmin)))
         enddo
        endif
      endif
      do i=1,ms(issmin)
         df(i)=df(i)/ffac
      enddo
      return
      end

C     Subroutines from "Numerical Recipes in FORTRAN" 2nd Ed.
C     (IBM PC diskette) Cambridge University Press, 1992
C     Licensed to use one copy for each subroutine to S. Yoneda
C
C     Only change made on "dfpmin" is following:
C     INCLUDE 'con0330p.f'
C     PARAMETER (NMAX=isslimit,ITMAX=200,STPMX=100d0,EPS=3.d-12,
C    >           TOLX=4.d0*EPS)
C     print *,'too many iterations in dfpmin'
C     iter=itmax+1
C
C     On "lnsrch":
C     PARAMETER (ALF=1.d-4,TOLX=1.d-12)
C
C     Also lines with "ires"

      SUBROUTINE dfpmin(p,n,gtol,iter,fret,this_func,my_dfunc,ires)
      INTEGER iter,n,NMAX,ITMAX,ires
      DOUBLE PRECISION fret,gtol,p(n),this_func,EPS,STPMX,TOLX
      INCLUDE 'con0330p.f'
      PARAMETER (NMAX=isslimit,ITMAX=200,STPMX=100d0,EPS=3.d-12,
     >           TOLX=4.d0*EPS)
     *
      EXTERNAL my_dfunc,this_func
CU    USES my_dfunc,this_func,lnsrch
      INTEGER i,its,j
      LOGICAL check
      DOUBLE PRECISION den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp
     *,test,dg(NMAX),
     *g(NMAX),hdg(NMAX),hessin(NMAX,NMAX),pnew(NMAX),xi(NMAX)
      fp=this_func(p)
      call my_dfunc(p,g)
      sum=0.d0
      do 12 i=1,n
        do 11 j=1,n
          hessin(i,j)=0.d0
11      continue
        hessin(i,i)=1.d0
        xi(i)=-g(i)
        sum=sum+p(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum), dble(n))
      do 27 its=1,ITMAX
        iter=its
        call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,this_func,ires)
        fp=fret
        do 13 i=1,n
          xi(i)=pnew(i)-p(i)
          p(i)=pnew(i)
13      continue
        do 15 i=1,n
          dg(i)=g(i)
15      continue
        call my_dfunc(p,g)
        test=0.d0
        den=max(fret,1.d0)
        do 16 i=1,n
          temp=abs(g(i))*max(abs(p(i)),1.d0)/den
          if(temp.gt.test)test=temp
16      continue
        if(test.lt.gtol) then
           ires=0
           return
        endif
        if (ires.eq.1) then
c           print *,'Restart Minimum Search Because Slope is Positive'
c           print *,'Iter=',iter,' F=',fp
c           print *,'X=',(p(i),i=1,n)
c           print *,'DF=',(g(i),i=1,n)
c           print *,'DX=',(xi(i),i=1,n)
           return
        endif
        test=0.d0
        do 14 i=1,n
          temp=abs(xi(i))/max(abs(p(i)),1.d0)
          if(temp.gt.test)test=temp
14      continue
        if(test.lt.TOLX)return
        do 17 i=1,n
          dg(i)=g(i)-dg(i)
17      continue
        do 19 i=1,n
          hdg(i)=0.d0
          do 18 j=1,n
            hdg(i)=hdg(i)+hessin(i,j)*dg(j)
18        continue
19      continue
        fac=0.d0
        fae=0.d0
        sumdg=0.d0
        sumxi=0.d0
        do 21 i=1,n
          fac=fac+dg(i)*xi(i)
          fae=fae+dg(i)*hdg(i)
          sumdg=sumdg+dg(i)**2
          sumxi=sumxi+xi(i)**2
21      continue
        if(fac**2.gt.EPS*sumdg*sumxi)then
          fac=1.d0/fac
          fad=1.d0/fae
          do 22 i=1,n
            dg(i)=fac*xi(i)-fad*hdg(i)
22        continue
          do 24 i=1,n
            do 23 j=1,n
              hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+
     *
     *fae*dg(i)*dg(j)
23          continue
24        continue
        endif
        do 26 i=1,n
          xi(i)=0.d0
          do 25 j=1,n
            xi(i)=xi(i)-hessin(i,j)*g(j)
25        continue
26      continue
27    continue
C      print *,'too many iterations in dfpmin'
      iter=itmax+1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software o)'%+!:591`$!.

      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,this_func,ires)
      INTEGER n,ires
      LOGICAL check
      DOUBLE PRECISION f,fold,stpmax,g(n),p(n),x(n),xold(n),this_func,ALF
     *,TOLX
      PARAMETER (ALF=1.d-4,TOLX=1.d-12)
      EXTERNAL this_func
CU    USES this_func
      INTEGER i
      DOUBLE PRECISION a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2
     *,slope,sum,temp,
     *test,tmplam
      check=.false.
      ires=0
      sum=0.d0
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.d0
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      if (slope.ge.0d0) then
         ires=1
         return
      endif
      test=0.d0
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.d0
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=this_func(x)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.d0)then
            tmplam=-slope/(2.d0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.d0)then
              tmplam=-slope/(2.d0*b)
            else
              disc=b*b-3.d0*a*slope
              if (disc.le.0) then
                 print *,disc,a,b,f,f2,alam,alam2,(x(i),i=1,n),slope,
     >                   fold,fold2
              endif
              tmplam=(-b+sqrt(disc))/(3.d0*a)
            endif
            if(tmplam.gt..5d0*alam)tmplam=.5d0*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1d0*alam)
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software o)'%+!:591`$!.
