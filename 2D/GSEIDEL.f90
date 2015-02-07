SUBROUTINE GSEIDEL(IL,JL,phi)
integer ii,jj
integer IL,JL
real, dimension :: phi




do jj = 2,JL-1
   do ii = 2,IL-1
   
      aw
      ae
      as
      phi(ii,jj) = (aw*phi(ii  ,jj-1) + ae*phi(ii,jj+1) + as*phi(ii-1,jj) +
                    an*phi(ii+1,jj  )/ap
   end do
end do

END SUBROUTINE GSEIDEL