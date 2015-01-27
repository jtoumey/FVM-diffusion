!*************************************************************************!
!                                                                         !
!  Module:       DIFFUSION1D.F90                                          !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         26 January 2015                                          !
!                                                                         !
!  Language:     FORTRAN90                                                !
!                                                                         !
!  Description:  This code solves source-free heat conduction in a 1-D    !
!                bar. The method follows example 4.3 in Versteeg and      !
!                Malalasekera, 2nd ed.                                    !
!                                                                         !
!*************************************************************************!
PROGRAM DIFFUSION1D
IMPLICIT NONE
!
!...number of nodes/CVs is `n'
!
integer n,ii,jj
parameter (n = 5)
real k,area,Ta,Tb,T(n)
real xmax,x(n),dx
real aw,ap,ae,Su,Sp
real a(n),b(n),c(n),d(n)
!
!...thermal conductivity [W/m.K] and cross-sectional area [m^2]
!
k    = 1000. 
area = 0.01
!
!...boundary condition temperatures [*C]
!
Ta = 100.
Tb = 500.  
!
!...set up the grid
!   xmax = domain length
!
xmax = 0.5
dx = xmax/float(n)
do jj = 1,n
  x(jj) = (jj-0.5)*dx
end do
!
!...set up system of equations
!   Left Boundary:
aw =  0.
ae =  k*area/dx
Su =  2.*k*area*Ta/dx
Sp = -2.*k*area/dx
ap =  aw + ae - Sp
!
a(1) = -aw
b(1) =  ap
c(1) = -ae
d(1) =  Su 
!...Interior cells
do ii = 2,n-1
   aw = k*area/dx
   ae = k*area/dx
   Su = 0.
   ap = ae + aw - Sp
   !
   a(ii) = -aw
   b(ii) =  ap
   c(ii) = -ae
   d(ii) =  Su
end do
!...Right Boundary
aw =  k*area/dx
ae =  0.
Su =  2.*k*area*Tb/dx
Sp = -2.*k*area/dx
ap =  aw + ae - Sp
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































END

