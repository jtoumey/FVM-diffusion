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

integer IL,JL,ii,jj,kk,iter
parameter (IL=3,JL=4)
real dx,dy,xmax,ymax,x(IL),y(JL)
real aw,ae,as,an,Su,Sp,ap
real qw,k,area,Tn
real a(JL),b(JL),c(JL),d(JL)
real Tsol(JL),Tprev,resid
real, dimension(JL,IL) :: T
!
!...Initial temperature distribution: 0 [*C] everywhere
!
T = 0.
resid = 1000.
iter = 1
!
!...set up the grid
!
xmax = 0.3
dx = xmax/float(IL)
do ii = 1,IL
   x(ii) = (ii-0.5)*dx
end do
!
ymax = 0.4
dy = ymax/float(JL)
do jj = 2,JL
   y(jj) = (jj-0.5)*dy
end do
!
!...physical properties
!
qw = 500000.
k  = 1000.
Tn = 100.
area = 0.001

do while (resid >= 1.)
   Tprev = maxval(T)
   !
   !...West Cells
   !
   !...Set up system 
   !   SW Corner
   aw = 0.
   ae = k*area/dx
   as = 0.
   an = k*area/dx
   Sp = 0.
   Su = area*qw
   ap = aw + ae + as + an - Sp
   !write(6,301)an,as,aw,ae,ap,Su
   !
   a(1) = -as
   b(1) =  ap
   c(1) = -an
   d(1) =  Su + ae*T(1,2)
   !
   do jj = 2,JL-1
      aw = 0.
      ae = k*area/dx
      as = k*area/dx
      an = k*area/dx
      Sp = 0.
      Su = area*qw
      ap = aw + ae + as + an - Sp
      !
      !write(6,301)an,as,aw,ae,ap,Su
      !
      a(jj) = -as
      b(jj) =  ap
      c(jj) = -an
      d(jj) =  Su + ae*T(jj,2)
   end do
   !...NW corner
   aw = 0.
   ae = k*area/dx
   as = k*area/dx
   an = 0.
   Sp = -2.*k*area/dx
   Su = area*qw + 2.*k*area*Tn/dx
   ap = aw + ae + as + an - Sp

   !write(6,301)an,as,aw,ae,ap,Su
   !
   a(JL) = -as
   b(JL) =  ap
   c(JL) = an
   d(JL) =  Su + ae*T(JL,2)

   !
   !...solve system with the TDMA
   !
   call thomas(JL,a,b,c,d,Tsol)
   !...Store temperature solution
T(:,1) = Tsol


do ii = 2,IL-1
   !-----------------------------------------------------------------------!
   !
   !...Interior cells
   !
   !...Set up system 
   !   South boundary
   aw = k*area/dx
   ae = k*area/dx
   as = 0.
   an = k*area/dx
   Sp = 0
   Su = 0
   ap = aw + ae + as + an - Sp
   !write(6,301)an,as,aw,ae,ap,Su
   !
   a(1) = -as
   b(1) =  ap
   c(1) = -an
   d(1) =  Su + ae*T(1,ii+1) + aw*T(1,ii-1)
   !
   do jj = 2,JL-1
      aw = k*area/dx
      ae = k*area/dx
      as = k*area/dx
      an = k*area/dx
      Sp = 0.
      Su = 0.
      ap = aw + ae + as + an - Sp
      !
      !write(6,301)an,as,aw,ae,ap,Su
      !
      a(jj) = -as
      b(jj) =  ap
      c(jj) = -an
      d(jj) =  Su + ae*T(jj,ii+1) + aw*T(jj,ii-1)
   end do
   !...North boundary
   aw = k*area/dx
   ae = k*area/dx
   as = k*area/dx
   an = 0.
   Sp = -2.*k*area/dx
   Su =  2.*k*area*Tn/dx
   ap = aw + ae + as + an - Sp
   !write(6,301)an,as,aw,ae,ap,Su
   !
   a(JL) = -as
   b(JL) =  ap
   c(JL) = -an
   d(JL) =  Su + ae*T(jj,ii+1) + aw*T(jj,ii-1)
   !
   !...solve system with the TDMA
   !
   call thomas(JL,a,b,c,d,Tsol)
   !...Store temperature solution
   T(:,ii) = Tsol

end do

   !
   !...East Cells
   !
   !...Set up system 
   !   SE Corner
   aw = k*area/dx
   ae = 0.
   as = 0.
   an = k*area/dx
   Sp = 0.
   Su = 0.
   ap = aw + ae + as + an - Sp
   !write(6,301)an,as,aw,ae,ap,Su
   !
   a(1) = -as
   b(1) =  ap
   c(1) = -an
   d(1) =  Su + aw*T(1,IL-1)
   !
   do jj = 2,JL-1
      aw = k*area/dx
      ae = 0.
      as = k*area/dx
      an = k*area/dx
      Sp = 0.
      Su = 0.
      ap = aw + ae + as + an - Sp
      !
      !write(6,301)an,as,aw,ae,ap,Su
      !
      a(jj) = -as
      b(jj) =  ap
      c(jj) = -an
      d(jj) =  Su + aw*T(jj,IL-1)
   end do
   !...NE corner
   aw = k*area/dx
   ae = 0.
   as = k*area/dx
   an = 0.
   Sp = -2.*k*area/dx
   Su = 2.*k*area*Tn/dx
   ap = aw + ae + as + an - Sp

   !write(6,301)an,as,aw,ae,ap,Su
   !
   a(JL) = -as
   b(JL) =  ap
   c(JL) = -an
   d(JL) =  Su + aw*T(JL,IL-1)
   !
   !...solve system with the TDMA
   !
   call thomas(JL,a,b,c,d,Tsol)
   !...Store temperature solution
   T(:,IL) = Tsol


   resid = maxval(T) - Tprev
   write(6,*)
   write(6,401)
   write(6,*)
   
end do

write(6,201)T




101 format(3x,f12.5,3x,f12.5)
201 format(3x,f12.5)
301 format(2x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5)
401 format(3x,'*** Iteration ***')
501 format(2x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5)
!do kk = 1,IL 
!   (6,101)x(kk),y(kk)
!end do


END