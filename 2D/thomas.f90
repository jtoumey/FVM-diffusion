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
real cp(n),dp(n),m
!
!...Forward sweep
! initialize c-prime and d-prime
cp(1) = c(1)/b(1)
dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
do ii = 2,n
   m = b(ii)-cp(ii-1)*a(ii)
   cp(ii) = c(ii)/m
   dp(ii) = (d(ii)-dp(ii-1)*a(ii))/m
enddo
! initialize x
phi(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
do ii = n-1,1,-1
   phi(ii) = dp(ii)-cp(ii)*phi(ii+1)
end do
END FUNCTION THOMAS