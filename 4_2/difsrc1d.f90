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
!                constant source. The method follows example 4.3 in       !
!                Versteeg and Malalasekera, 2nd ed.                       !
!                                                                         !
!*************************************************************************!
PROGRAM DIFSRC1D
IMPLICIT NONE

integer n,ii,jj
parameter (n = 5)
real k,q,area,Ta,Tb,T(n)
real xmax,x(n),dx
real aw,ap,ae,Su,Sp
real a(n),b(n),c(n),d(n)
!
!...thermal conductivity [W/m.K] and uniform heat generation [kW/m^3]
!   area [m^2]
!
k    = 0.05
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
dx = xmax/float(n)
!write(6,201)dx
do jj = 1,n
  x(jj) = (jj-0.5)*dx
end do
!
!...set up system of equations
!   Left Boundary:
!
aw = 0.
ae = k*area/dx
Su =  2.*k*area*Ta/dx + q*area*dx
Sp = -2.*k*area/dx
ap =  ae + aw - Sp
!
a(1) = -aw
b(1) =  ap
c(1) = -ae
d(1) =  Su
!...Interior cells:
do ii = 2,n-1
   aw = 0.
   ae = k*area/dx
   Su =  2.*k*area*Ta/dx + q*area*dx
   Sp = -2.*k*area/dx
   ap =  ae + aw - Sp
   !
   a(1) = -aw
   b(1) =  ap
   c(1) = -ae
   d(1) =  Su
end do 
!...Right Boundary:
aw = k*area/dx
ae = 0.
Su =  2.*k*area*Tb/dx + q*area*dx
Sp = -2.*k*area/dx
ap =  ae + aw - Sp
!
a(n) = -aw
b(n) =  ap
c(n) = -ae
d(n) =  Su
























201      format(3x,f12.5)!,3x,f12.5)

END