!======================================================================
      Real(8) FUNCTION WKy (i1,j1,i2,j2,k)
!======================================================================
!                 k                          k
!     Evaluates  W (i1, j1; i2, j2) =   2   V   (i1, j1; i2, j2)
!                                            k+1
!                                    -(k+1) N   (i1, j1; i2, j2) 
!                                            k-1
!                                    +(k+2) N   (j1, i1, j2, i2) 
!----------------------------------------------------------------------
      IMPLICIT NONE
    
      Integer, intent(in) :: i1,j1,i2,j2,k
    
      Real(8), external :: VKy, NKy
    
      Wky = 2 * VKy(i1,j1,i2,j2,k) - (k+1) * NKy(i1,j1,i2,j2,k+1)  &
                                   + (k+2) * NKy(j1,i1,j2,i2,k-1)
    
      END FUNCTION WKy
    
    
 