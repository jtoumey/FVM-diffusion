!*************************************************************************!
!                                                                         !
!  Module:       DIFFUSION2D.F90                                          !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         January 2015                                             !
!                                                                         !
!  Language:     FORTRAN90                                                !
!                                                                         !
!  Description:  This code solves heat conduction in a 2-D plate. with a  !
!                constant source. The method follows example 4.2 in       !
!                Versteeg and Malalasekera, 2nd ed.                       !
!                                                                         !
!*************************************************************************!
PROGRAM DIFFUSION2D
IMPLICIT NONE
! n = number of CVs + 2 (addt'l index on each side for the node on the boundary. This way, we can add an 
! additional equation for the Temperature at each boundary and modify some coefficients in order to 
! quickly specify any type of BC without re-formulating the whole problem and re-coding the B.C. 
! implementation)

integer IL,JL,ii,jj,kk
parameter (IL=3,JL=4)
real dx,dy,xmax,ymax,x(IL),y(JL)
real aw,ae,as,an,Su,Sp,ap
real qw,k,area,Tn
real a(JL
real, dimension(JL,1) :: T


T = 0.
!
!...set up the grid
!
xmax = 0.3
dx = xmax/float(IL)
!x(1) = 0.
do ii = 1,IL
   x(ii) = (ii-0.5)*dx
end do
!x(IL) = xmax
!
ymax = 0.4
dy = ymax/float(JL)
!y(1) = 0
do jj = 2,JL
   y(jj) = (jj-0.5)*dy
end do
!y(JL) = ymax
!
!...physical properties
!
qw = 500000.
k = 1000.
Tn = 100.
!
!...Set up system 
!   SW Corner
aw = 0.
ae = k*area/dx
as = 0.
an = k*area/dx
Sp = 0
Su = area*qw
ap = aw + ae + as + an + Sp
!
a(1) =  0.
b(1) =  ap
c(1) = -an
d(1) =  Su + ae*T(2,1)
!
do jj = 2,JL-1
   aw = 0
   ae = k*area/dx
   as = k*area/dx
   an = k*area/dx
   Sp = 0
   Su = area*qw
   ap = aw + ae + as + an + Sp
   !
   a(1) = -as
   b(1) =  ap
   c(1) = -an
   d(1) =  Su + ae*T(2,JL)
end do
!...NW corner
aw = 0.
ae = k*area/dx
as = 0.
an = k*area/dx
Sp = 0
Su = area*qw
ap = aw + ae + as + an + Sp
!
a(1) =  0.
b(1) =  ap
c(1) = -an
d(1) =  Su + ae*T(2,1)
!



101 format(3x,f12.5,3x,f12.5)

!do kk = 1,IL 
!   write(6,101)x(kk),y(kk)
!end do


END
