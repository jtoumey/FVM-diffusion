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
parameter (IL=6,JL=6)
real dx,dy,xmax,ymax,x(IL),y(JL)


real, dimension(IL,JL) :: T
!
!...set up the grid
!
xmax = 0.3
dx = xmax/float(IL-2)
x(1) = 0.
do ii = 2,IL-1
   x(ii) = (ii-1-0.5)*dx
end do
x(IL) = xmax
!
ymax = 0.4
dy = ymax/float(JL-2)
y(1) = 0
do jj = 2,JL-1
   y(jj) = (jj-1-0.5)*dy
end do
y(JL) = ymax
!
!...Set up system 
!






101 format(3x,f12.5,3x,f12.5)

do kk = 1,IL 
   write(6,101)x(kk),y(kk)
end do


END
