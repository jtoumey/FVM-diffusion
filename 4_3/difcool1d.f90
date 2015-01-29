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
Ta = 100.
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
Su =  2.*dx*Ta + n2*Tinf*dx
Sp = -n2*dx - 2./dx
ap = aw + ae - Sp 
!
a(1) = -aw
b(1) =  ap
c(1) =  ae
d(1) =  Su
!...Interior cells














END