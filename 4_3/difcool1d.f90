!*************************************************************************!
!                                                                         !
!  Module:       DIFCOOL1D.F90                                            !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         28 January 2015                                          !
!                                                                         !
!  Language:     FORTRAN90                                                !
!                                                                         !
!  Description:  This code solves cooling of a circular fin by convective !
!                heat transfer along its length. Base (left boundary) has !
!                a fixed temperature. Right boundary is insulated (zero   !
!                normal gradient). The method follows example 4.3 in      !
!                Versteeg and Malalasekera, 2nd ed.                       !
!                                                                         !
!*************************************************************************!
PROGRAM DIFCOOL1D
IMPLICIT NONE 

integer n,ii,jj
parameter (n = 5)
real n2,Ta,Tinf,T(n)
real xmax,x(n),dx
real aw,ap,ae,Su,Sp
real a(n),b(n),c(n),d(n)
!
!...n^2 = hP/kA [1/m^2]
!
n2 = 25.
!
!...boundary temperature: fixed at left [*C]
!   freestream temp [*C]
!
Ta   = 100.
Tinf = 20.
!
!...set up the grid
!   xmax is the domain length [m]
!
xmax = 1.
dx = xmax/float(n)
do jj = 1,n
  x(jj) = (jj-0.5)*dx
end do
!
!...set up system of equations
!   Left Boundary:
!
aw =  0.
ae =  1./dx
Su =  2.*Ta/dx + n2*Tinf*dx
Sp = -n2*dx - 2./dx
ap =  aw + ae - Sp 
!
a(1) = -aw
b(1) =  ap
c(1) = -ae
d(1) =  Su
!...Interior cells
do ii = 2,n-1
   aw =  1./dx
   ae =  1./dx
   Su =  n2*Tinf*dx
   Sp = -n2*dx
   ap =  aw + ae - Sp 
   !
   a(ii) = -aw
   b(ii) =  ap
   c(ii) = -ae
   d(ii) =  Su
end do
!....Right boundary
aw =  1./dx
ae =  0.
Su =  n2*Tinf*dx
Sp = -n2*dx
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
101 format(5x,'____x(j)___',3x,'____T(j)__')
201 format(2x,f12.5,2x,f12.5)
write(6,101)
open(unit=7,file='temp_distr.dat')
write(6,201)0.,Ta
write(7,201)0.,Ta
do jj = 1,n
   write(6,201)x(jj),T(jj)
   write(7,201)x(jj),T(jj)
end do























END