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
dx = xmax/float(n)
do jj = 1,n
  x(jj) = (jj-0.5)*dx
end do
!
!...set up system of equations
!   Left Boundary:
!
aw =  0.
ae =  k*area/dx
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
!...Right Boundary:
aw =  k*area/dx
ae =  0.
Su =  2.*k*area*Tb/dx + q*area*dx
Sp = -2.*k*area/dx
ap =  ae + aw - Sp
!
a(n) = -aw
b(n) =  ap
c(n) = -ae
d(n) =  Su
!
!...Solve the system using the Thomas Algorithm 
!
call thomas(n,a,b,c,d,T)
!
!...Write results
!
101   format(5x,'____x(j)___',3x,'____T(j)__')
201      format(2x,f12.5,2x,f12.5)
write(6,101)
open(unit=7,file='temp_distr.dat')
write(6,201)0.,Ta
write(7,201)0.,Ta
do jj = 1,n
   write(6,201)x(jj),T(jj)
   write(7,201)x(jj),T(jj)
end do
write(6,201)xmax,Tb
write(7,201)xmax,Tb



END