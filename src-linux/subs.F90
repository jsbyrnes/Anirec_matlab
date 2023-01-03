SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)

IMPLICIT NONE

INTEGER(kint) I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
REAL(kreal) A(6,6),SCALE(6)
REAL(kreal) C,F,G,R,S,B2,RADIX
REAL(kreal) DABS

LOGICAL NOCONV

do 900 I=1,6
do 901 J=1,6
if(ISNAN(A(I,J))) A(I,J)=0.0
901 continue
900 continue

!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
!     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).1

!     THIS SUBROUTINE BALANCES A REAL MATRIX AND ISOLATES
!     EIGENVALUES WHENEVER POSSIBLE.

!     ON INPUT:

!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT;

!        N IS THE ORDER OF THE MATRIX;

!        A CONTAINS THE INPUT MATRIX TO BE BALANCED.

!     ON OUTPUT:

!        A CONTAINS THE BALANCED MATRIX;

!        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)
!          IS EQUAL TO ZERO IF
!           (1) I IS GREATER THAN J AND
!           (2) J=1,...,LOW-1 OR I=IGH+1,...,N;

!        SCALE CONTAINS INFORMATION DETERMINING THE
!           PERMUTATIONS AND SCALING FACTORS USED.

!     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
!     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
!     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
!     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
!        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
!                 = D(J,J),      J = LOW,...,IGH
!                 = P(J)         J = IGH+1,...,N.
!     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
!     THEN 1 TO LOW-1.

!     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.

!     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
!     BALANC  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
!     K,L HAVE BEEN REVERSED.)

!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY

!     ------------------------------------------------------------------
!
!     :::::::::: RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.
!                RADIX = 16.0D0 FOR LONG FORM ARITHMETIC
!                ON S360 ::::::::::
DATA RADIX/2/

B2 = RADIX * RADIX
K = 1
L = N

GO TO 100
!     :::::::::: IN-LINE PROCEDURE FOR ROW AND
!                COLUMN EXCHANGE ::::::::::
20 SCALE(M) = J
IF (J .EQ. M) GO TO 50
!
DO 30 I = 1, L
F = A(I,J)
A(I,J) = A(I,M)
A(I,M) = F
30 CONTINUE

DO 40 I = K, N
F = A(J,I)
A(J,I) = A(M,I)
A(M,I) = F
40 CONTINUE

50 GO TO (80,130), IEXC
!     :::::::::: SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                AND PUSH THEM DOWN ::::::::::
80 IF (L .EQ. 1) GO TO 280
L = L - 1
!     :::::::::: FOR J=L STEP -1 UNTIL 1 DO -- ::::::::::
100 DO 120 JJ = 1, L
J = L + 1 - JJ

DO 110 I = 1, L
IF (I .EQ. J) GO TO 110
IF (A(J,I) .NE. 0.0D0) GO TO 120
110    CONTINUE

M = L
IEXC = 1
GO TO 20
120 CONTINUE

GO TO 140
!     :::::::::: SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
!                AND PUSH THEM LEFT ::::::::::
130 K = K + 1

140 DO 170 J = K, L

DO 150 I = K, L
IF (I .EQ. J) GO TO 150
IF (A(I,J) .NE. 0.0D0) GO TO 170
150    CONTINUE

M = K
IEXC = 2
GO TO 20
170 CONTINUE
!     :::::::::: NOW BALANCE THE SUBMATRIX IN ROWS K TO L ::::::::::
DO 180 I = K, L
180 SCALE(I) = 1.0D0
!     :::::::::: ITERATIVE LOOP FOR NORM REDUCTION ::::::::::

190 NOCONV = .FALSE.

DO 270 I = K, L
C = 0.0D0
R = 0.0D0

DO 200 J = K, L
!IF (J .EQ. I) GO TO 200
C = C + DABS(A(J,I))
R = R + DABS(A(I,J))
200    CONTINUE

!     :::::::::: GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ::::::::::
IF (C .EQ. 0.0D0 .OR. R .EQ. 0.0D0) GO TO 270
G = R / RADIX
F = 1.0D0
S = C + R

210    IF (C .GE. G) GO TO 220
F = F * RADIX
C = C * B2
GO TO 210
220    G = R * RADIX
230    IF (C .LT. G) GO TO 240
F = F / RADIX
C = C / B2
GO TO 230
!     :::::::::: NOW BALANCE ::::::::::
240    IF ((C + R) / F .GE. 0.95D0 * S) GO TO 270

G = 1.0D0 / F
SCALE(I) = SCALE(I) * F
NOCONV = .TRUE.

DO 250 J = K, N
250    A(I,J) = A(I,J) * G

DO 260 J = 1, L
260    A(J,I) = A(J,I) * F

270 CONTINUE

IF (NOCONV) GO TO 190

280 LOW = K
IGH = L

RETURN
!     :::::::::: LAST CARD OF BALANC ::::::::::
END



subroutine csolve(nn,a,ai,x,xi,y,yi)
!  solves the complex nxn system of equations a*x=y using gaussian elimination
!  and partial pivoting

integer(kint), intent(in)           :: nn
real(kreal), intent(inout)          :: a(6, 6), ai(6, 6), x(6), xi(6), y(6), yi(6)

integer(kint)                       :: n

    n=nn


    if(n.gt.0) then
        call clup(n,a,ai,ip)
    else
        n=-n
    endif

    call cbcktr(n,a,ai,x,xi,y,yi,ip)

    return

end



subroutine clup(n,a,ai,ip)
!  finds lu decomp of a+i*ai using partial pivoting
!  pivoting sequence returned in ip(n)

integer(kint), intent(in)          :: n
integer(kint), intent(inout)       :: ip(1000)
real(kreal), intent(inout)         :: a(6,6), ai(6,6)

real(kreal)                        :: tol, aam, aai, tem, b, bi, 
integer(kint)                      :: nm, i, j, jm, i1, ipi, ipj,

tol=1.d-14

!  initialize permutation vector
do 50 i=1,n
50 ip(i)=i

nm1=n-1

if(n.eq.1) go to 700
do 100 i=1,nm1
aam=0.d0
do 200 j=i,n
aai=a(ip(j),i)**2+ai(ip(j),i)**2
if(aam.gt.aai) go to 200
aam=aai
ipi=ip(j)
jm=j
200 continue
if(aam.lt.tol) go to 400
ip(jm)=ip(i)
ip(i)=ipi
i1=i+1
do 100 j=i1,n
ipj=ip(j)
! if victim index is already zero, dont bother to rub it out
tem=dabs(a(ipj,i))+dabs(ai(ipj,i))
if(tem.lt.tol) go to 100
b=(a(ipj,i)*a(ipi,i)+ai(ipj,i)*ai(ipi,i))/aam
bi=(ai(ipj,i)*a(ipi,i)-a(ipj,i)*ai(ipi,i))/aam
a(ipj,i)=b
ai(ipj,i)=bi
do 500 k=i1,n
a(ipj,k)=a(ipj,k)-b*a(ipi,k)+bi*ai(ipi,k)
500 ai(ipj,k)=ai(ipj,k)-b*ai(ipi,k)-bi*a(ipi,k)
100 continue
700 continue
return
!c  400 print 101,aam,i
400 continue
!c  101 format('near-zero pivot ',e12.5,'  on column',i3)
!c      stop
end


subroutine cbcktr(n1,z,zi,dr,di,er,ei,ip)

integer(kint), intent(in)           :: n1
integer(kreal), intent(inout)       :: z(6,6), zi(6,6), dr(6), di(6), er(6), ei(6)
integer(kint), intent(inout)        :: ip(1000)

integer(kint)                       :: i, j, i1, ii, ip1
real(kreal)                         :: zkk1, zkk2, uii, dri, dii

!  back transform with unit lower triangular matrix
do 300 i=1,n1
dr(i)=er(ip(i))
300 di(i)=ei(ip(i))
if(n1.eq.1) go to 400
do 310 i=2,n1
i1=i-1
do 310 j=1,i1
zkk1=z(ip(i),j)
zkk2=zi(ip(i),j)
dr(i)=dr(i)-zkk1*dr(j)+zkk2*di(j)
310 di(i)=di(i)-zkk1*di(j)-zkk2*dr(j)
400 continue
!  back transform with upper triangular matrix
do 320 ii=1,n1
i=n1+1-ii
ip1=i+1
uii=z(ip(i),i)**2+zi(ip(i),i)**2
if(i.eq.n1) go to 340
do 330 j=ip1,n1
zkk1=z(ip(i),j)
zkk2=zi(ip(i),j)
dr(i)=dr(i)-zkk1*dr(j)+zkk2*di(j)
330 di(i)=di(i)-zkk1*di(j)-zkk2*dr(j)
340 dri=dr(i)
dii=di(i)
zkk1=z(ip(i),i)
zkk2=zi(ip(i),i)
dr(i)=(dri*zkk1+dii*zkk2)/uii
320 di(i)=(dii*zkk1-dri*zkk2)/uii


return
end

subroutine solve(nn,a,x,y)
!  solves the nxn system of equations a*x=y using gaussian elimination
!  and partial pivoting
!  if n<0 the lu decomposition is already done
!  note that the matrix a is modified

implicit none

real(kreal), intent(inout) a(3,3), x(3), y(3)
integer(kint) nn, ip(1000)

n=nn
if(n.gt.0)then
call lup(n,a,ip)
else
n=-n
endif
!      print 101,((i,j,a(i,j),j=1,n),i=1,n)
!  101 format(' a(',2i2,')=',e15.5)
!      type 102,(ip(i),i=1,n)
!  102 format(10i5)

call bcktr(n,a,x,y,ip)
return


end


subroutine lup(n,a,ip)

!  finds lu decomp of a using partial pivoting
!  output in c (upper triangle/unit lower triangle) and
!  pivoting sequence returned in ip(n)

integer(kint), intent(in)     :: n
integer(kint), intent(inout)  :: ip(1000)
real(kreal), intent(inout)    :: a(3,3)

real(kreal)                   :: tol, aam, aai, tem, b
integer(kint)                 :: nm1, i, j, jm, ipi, i1, ipj

tol=1.d-14
do 50 i=1,n
50 ip(i)=i
nm1=n-1
if(n.eq.1) go to 700
do 100 i=1,nm1
aam=0.d0
do 200 j=i,n
aai=a(ip(j),i)**2
if(aam.gt.aai) go to 200
aam=aai
ipi=ip(j)
jm=j
200 continue
if(aam.lt.tol) go to 400
ip(jm)=ip(i)
ip(i)=ipi
i1=i+1
do 100 j=i1,n
ipj=ip(j)
!  if victim index is already zero, dont bother to rub it out
tem=dabs(a(ipj,i))
if(tem.lt.tol) go to 100
b=(a(ipj,i)*a(ipi,i))/aam
a(ipj,i)=b
do 500 k=i1,n
a(ipj,k)=a(ipj,k)-b*a(ipi,k)
500 continue
100 continue
700 continue
return
!  400 print 101,aam,i
400 continue
!  101 format('near-zero pivot ',e12.5,'  on column',i3)
!      stop


end


subroutine bcktr(n1,z,dr,er,ip)

!  performs backtransform on input vector er - 'y'
!  to find solution dr - 'x'

real(kreal), intent(inout) :: z(3,3), dr(3), er(3)
integer(kint), intent(inout) :: ip(1000), n1

integer(kint) i, i1, ii, ip1, j

!  back transform with unit lower triangular matrix

do 300 i=1,n1
300 dr(i)=er(ip(i))
if(n1.eq.1) go to 400
do 310 i=2,n1
i1=i-1
do 310 j=1,i1
310 dr(i)=dr(i)-z(ip(i),j)*dr(j)
400 continue
!  back transform with upper triangular matrix
do 320 ii=1,n1
i=n1+1-ii
ip1=i+1
if(i.eq.n1) go to 320
do 330 j=ip1,n1
330 dr(i)=dr(i)-z(ip(i),j)*dr(j)
320 dr(i)=dr(i)/z(ip(i),i)
return

end