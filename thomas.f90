!*************************************************************************!
!                                                                         !
!  Module:       THOMAS.F90                                               !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         26 January 2015                                          !
!                                                                         !
!  Language:     FORTRAN90                                                !
!                                                                         !
!  Description:  Solves a tri-diagonal system of linear equations using   !
!                the Thomas Algorithm (TDMA). The algorithm in this       !
!                function comes from Versteeg and Malalasekera, 2nd ed.,  !
!                Section 7.2                                              !
!                                                                         !
!                  a     Sub diagonal                                     !
!                  b     Main diagonal                                    !
!                  c     Super diagonal                                   !
!                  d     Right side of Linear System                      !
!                  phi   Solution vector                                  !
!                                                                         !
!*************************************************************************!
FUNCTION THOMAS(n,a,b,c,d,phi)
!
integer n,ii
real a(n),b(n),c(n),d(n)
real phi(n)
!
!...Forward sweep
!
do ii = 2,n
   d(ii) = d(ii) - (b(ii)*a(ii-1))/d(ii-1)
   c(ii) = c(ii) - (b(ii)*c(ii-1))/d(ii-1)
end do
!
!...Backward sweep
!
phi(n) = c(n)/d(n)
!
do ii = n-1,1,-1
   phi(ii) = (c(ii) - a(ii)*phi(ii+1))/d(ii)
end do
END FUNCTION THOMAS