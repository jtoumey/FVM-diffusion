!*************************************************************************!
!                                                                         !
!  Module:       DIFSRC1D.F90                                             !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         27 January 2015                                          !
!                                                                         !
!  Language:     FORTRAN90                                                !
!                                                                         !
!  Description:  This code solves heat conduction in a 1-D plate with a   !
!                constant source. The method follows example 4.2 in       !
!                Versteeg and Malalasekera, 2nd ed.                       !
!                                                                         !
!                ** Modified to include different B.C. formulation **     !
!                ** from Dr. Blaisdell's AAE412 class notes        **     !
!                                                                         !
!*************************************************************************!
PROGRAM DIFSRC1D
IMPLICIT NONE
!...n-2 is the number of cells -- added two extra to account for B.C.
integer n,ii,jj
parameter (n = 7)
real k,q,area,Ta,Tb,T(n)
real xmax,x(n),dx
real aw,ap,ae,Su,Sp,aa,ab
real a(n),b(n),c(n),d(n)
real alphaA,betaA,gammaA,alphaB,betaB,gammaB
!
!...thermal conductivity [W/m.K] and uniform heat generation [kW/m^3]
!   area [m^2]
!
k    = 0.5
q    = 1000000.
area = 1.
!
!...boundary temperatures: fixed [*C]
!
Ta = 100.
Tb = 200.
!
!...set up the grid
!   xmax is the domain length [m]
!
xmax = 0.02
dx = xmax/float(n-2)
x(1) = 0.0
do jj = 2,n-1
  x(jj) = (jj-1-0.5)*dx
end do
x(n) = xmax
!
!...set up system of equations
!                                  [alphaA betaA gammaA]
!   Dirichlet: T_A = const         [1 0 Ta]
!   Neumann  : (dT/dx)_A = -q      [0 1 q ] 
!   Robin    : hT_A-k(dT/dx) = hT  [h k hT]
!...Left side: 
alphaA = 1.
betaA  = 0.
gammaA = Ta
!...Right side
alphaB = 0.
betaB  = 1.
gammaB = 0 
!
!...Left boundary 
a(1) = 0.
b(1) = alphaA - 2.*betaA/dx
c(1) = 2.*betaA/dx
d(1) = gammaA
!...Left Cell
aw = 0.0 
aa = 2.*k*area/dx
ae = k*area/dx
Su = q*area*dx
Sp = 0.0
ap = ae + aa - Sp
!
a(2) = -aa
b(2) =  ap
c(2) = -ae
d(2) =  Su
!...Interior cells:
do ii = 3,n-2
   aw = k*area/dx
   ae = k*area/dx
   Su = q*area*dx
   Sp = 0.0
   ap = ae + aw - Sp
   !
   a(ii) = -aw
   b(ii) =  ap
   c(ii) = -ae
   d(ii) =  Su
end do 
!...Right cell
aw = k*area/dx
ae = 0.
ab = 2.*k*area/dx
Su = q*area*dx
Sp = 0
ap = aw + ab - Sp
!
a(n-1) = -aw
b(n-1) =  ap
c(n-1) = -ab
d(n-1) =  Su
!
!...Right Boundary
a(n) = -2.*betaB/dx
b(n) = alphaB + 2.*betaB/dx
c(n) = 0.
d(n) = gammaB
!
!...Solve the system using the Thomas Algorithm 
!
call thomas(n,a,b,c,d,T)
!
!...Write results
!
101 format(5x,'____x(j)___',3x,'____T(j)__')
201 format(2x,f12.5,2x,f12.5)
write(6,101)
open(unit=7,file='temp_distr.dat')
do jj = 1,n
   write(6,201)x(jj),T(jj)
   write(7,201)x(jj),T(jj)
end do
END
