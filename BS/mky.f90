!=======================================================================
      Real(8) Function MKy (I1,J1,I2,J2,K)
!=======================================================================
!                 k
!     Evaluates  M (i1, j1; i2, j2)  = 
!
!                 k                      k
!                N (i1, j1; i2, j2)  +  N (j1, i1; j2, i2)    *  1/2  ?? 
!
!------------------------------------------------------------------------


      IMPLICIT NONE

      INTEGER, INTENT(in) :: I1,J1,I2,J2,K
      REAL(8), EXTERNAL :: NKy

      MKy = NKy(i1,j1,i2,j2,k) + Nky(j1,i1,j2,i2,k)
      
      MKy = MKy / 2.d0   ! ???
      
      END Function MKy
      

